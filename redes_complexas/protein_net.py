"""
protein_net.py
==============
Pipeline de detecção de domínios funcionais em proteínas via redes complexas.

Fluxo: PDB -> rede de contatos de resíduos -> topologia (graus, centralidades)
-> detecção de comunidades (Louvain) -> validação (ARI/NMI vs anotação).

Autor: Pedro (PPGCO/FACOM-UFU) | Disciplina: Redes Complexas (Prof. Gilberto Nakamura)
Uso de IA declarado em declaracao_ia.md.

Modelo de rede (justificativa resumida; detalhes no relatório):
- Nó  = um resíduo de aminoácido (identificado por cadeia + número do resíduo).
- Aresta = dois resíduos "em contato", i.e. com pelo menos um par de átomos
  pesados a distância <= cutoff (modo 'heavy', padrão 5.0 A), ou cujos Cα estão
  a <= cutoff (modo 'ca', padrão 8.0 A).
- Peso  = número de pares de átomos em contato (força do contato) no modo 'heavy';
  1/d_CC no modo 'ca'. Resíduos consecutivos (ligação peptídica) são mantidos.

A escolha de granularidade por RESÍDUO (e não por átomo) é o que torna o método
viável para alvos grandes como 6B1T: a rede passa de ~10^6 átomos para ~10^4-10^5
resíduos, preservando a topologia de domínios.
"""

from __future__ import annotations
import os
import warnings
import urllib.request
from collections import defaultdict

import numpy as np
import networkx as nx

warnings.filterwarnings("ignore")

# ----------------------------------------------------------------------------- #
# 1. Obtenção e leitura do PDB
# ----------------------------------------------------------------------------- #
def fetch_pdb(pdb_id: str, out_dir: str = "data", url: str | None = None) -> str:
    """Baixa um PDB. Use `url` para fontes alternativas (ex. GitHub raw) quando o
    RCSB não estiver acessível. Para o RCSB padrão, deixe url=None.

    Retorna o caminho local do arquivo .pdb (ou .cif)."""
    os.makedirs(out_dir, exist_ok=True)
    pdb_id = pdb_id.upper()
    dst = os.path.join(out_dir, f"{pdb_id}.pdb")
    if os.path.exists(dst):
        return dst
    if url is None:
        # RCSB oficial (funciona em máquina com acesso a rcsb.org)
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    urllib.request.urlretrieve(url, dst)
    return dst


def load_structure(path: str, pdb_id: str | None = None):
    """Lê o arquivo com o Biopython. Aceita .pdb e .cif."""
    from Bio.PDB import PDBParser, MMCIFParser
    pid = pdb_id or os.path.splitext(os.path.basename(path))[0]
    parser = MMCIFParser(QUIET=True) if path.lower().endswith(".cif") else PDBParser(QUIET=True)
    return parser.get_structure(pid, path)


def _residue_key(residue):
    """Chave estável e legível para um resíduo: 'CHAIN:RESNAME:RESSEQ'."""
    chain = residue.get_parent().id
    hetflag, resseq, icode = residue.id
    return f"{chain}:{residue.resname}:{resseq}{icode.strip()}"


# ----------------------------------------------------------------------------- #
# 2. Construção da rede de contatos
# ----------------------------------------------------------------------------- #
def build_residue_network(structure, mode: str = "heavy",
                          cutoff: float | None = None,
                          model_index: int = 0,
                          weighted: bool = True) -> nx.Graph:
    """Constrói a rede de contatos de resíduos.

    Parameters
    ----------
    mode : 'heavy' (contato por átomos pesados, padrão) ou 'ca' (Cα-Cα).
    cutoff : limiar em Angstrom. Padrão 5.0 ('heavy') ou 8.0 ('ca').
    model_index : qual modelo usar (NMR pode ter vários; usamos o primeiro).
    weighted : se True, atribui peso às arestas.

    A busca de vizinhos usa KD-tree (Bio.PDB.NeighborSearch) -> O(n) na prática,
    essencial para escalar até proteínas grandes.
    """
    from Bio.PDB import NeighborSearch
    from Bio.PDB.Polypeptide import is_aa

    if cutoff is None:
        cutoff = 5.0 if mode == "heavy" else 8.0

    model = list(structure)[model_index]

    # Coleta resíduos de aminoácido (descarta água/ligantes/heteroátomos)
    residues = [r for r in model.get_residues() if is_aa(r, standard=False)]
    G = nx.Graph()
    for r in residues:
        key = _residue_key(r)
        G.add_node(key,
                   chain=r.get_parent().id,
                   resname=r.resname,
                   resseq=r.id[1])
        # coordenada do Cα (ou centroide) p/ visualização 3D
        if "CA" in r:
            ca = r["CA"].coord
        else:
            coords = np.array([a.coord for a in r.get_atoms()])
            ca = coords.mean(axis=0)
        G.nodes[key]["x"], G.nodes[key]["y"], G.nodes[key]["z"] = map(float, ca)

    res_of_atom = {}
    atoms = []
    if mode == "ca":
        for r in residues:
            if "CA" in r:
                a = r["CA"]; atoms.append(a); res_of_atom[id(a)] = _residue_key(r)
        ns = NeighborSearch(atoms)
        for a, b in ns.search_all(cutoff):
            ka, kb = res_of_atom[id(a)], res_of_atom[id(b)]
            if ka == kb:
                continue
            d = float(np.linalg.norm(a.coord - b.coord))
            w = 1.0 / d if weighted else 1.0
            G.add_edge(ka, kb, weight=w, dist=d)
    else:  # heavy
        for r in residues:
            rk = _residue_key(r)
            for a in r.get_atoms():
                if a.element == "H":
                    continue
                atoms.append(a); res_of_atom[id(a)] = rk
        ns = NeighborSearch(atoms)
        contacts = defaultdict(int)        # (ka,kb) -> nº de pares atômicos
        mindist = {}
        for a, b in ns.search_all(cutoff):
            ka, kb = res_of_atom[id(a)], res_of_atom[id(b)]
            if ka == kb:
                continue
            e = (ka, kb) if ka < kb else (kb, ka)
            contacts[e] += 1
            d = float(np.linalg.norm(a.coord - b.coord))
            mindist[e] = min(mindist.get(e, 1e9), d)
        for (ka, kb), n in contacts.items():
            w = float(n) if weighted else 1.0
            G.add_edge(ka, kb, weight=w, dist=mindist[(ka, kb)], n_contacts=n)
    return G


def giant_component(G: nx.Graph) -> nx.Graph:
    """Retorna o maior componente conexo (subgrafo)."""
    if G.number_of_nodes() == 0:
        return G
    comp = max(nx.connected_components(G), key=len)
    return G.subgraph(comp).copy()


# ----------------------------------------------------------------------------- #
# 3. Análise topológica
# ----------------------------------------------------------------------------- #
def topology_summary(G: nx.Graph) -> dict:
    """Métricas globais da rede."""
    n, m = G.number_of_nodes(), G.number_of_edges()
    degs = [d for _, d in G.degree()]
    comps = list(nx.connected_components(G))
    s = {
        "n_nodes": n,
        "n_edges": m,
        "density": nx.density(G),
        "avg_degree": float(np.mean(degs)) if degs else 0.0,
        "max_degree": int(np.max(degs)) if degs else 0,
        "n_components": len(comps),
        "giant_frac": (max(len(c) for c in comps) / n) if n else 0.0,
        "avg_clustering": nx.average_clustering(G),
    }
    Gc = giant_component(G)
    # diâmetro/caminho médio só no componente gigante e se não for grande demais
    if Gc.number_of_nodes() <= 5000 and Gc.number_of_nodes() > 1:
        s["avg_shortest_path"] = nx.average_shortest_path_length(Gc)
        s["diameter"] = nx.diameter(Gc)
    else:
        s["avg_shortest_path"] = None
        s["diameter"] = None
    return s


def centralities(G: nx.Graph, betweenness_k: int | None = None,
                 seed: int = 42) -> "pd.DataFrame":
    """Calcula centralidades por nó. Para grafos grandes, `betweenness_k` ativa
    a estimativa por amostragem (k pivôs) — o exato é O(nm)."""
    import pandas as pd
    deg = dict(G.degree())
    deg_c = nx.degree_centrality(G)
    n = G.number_of_nodes()
    if betweenness_k is None:
        betweenness_k = None if n <= 1500 else min(500, n)
    btw = nx.betweenness_centrality(G, k=betweenness_k, weight="dist",
                                    normalized=True, seed=seed)
    Gc = giant_component(G)
    clo = nx.closeness_centrality(Gc)
    clo = {x: clo.get(x, 0.0) for x in G.nodes()}
    try:
        eig = nx.eigenvector_centrality_numpy(G, weight="weight")
    except Exception:
        eig = {x: np.nan for x in G.nodes()}
    df = pd.DataFrame({
        "node": list(G.nodes()),
        "chain": [G.nodes[x].get("chain") for x in G.nodes()],
        "degree": [deg[x] for x in G.nodes()],
        "degree_centrality": [deg_c[x] for x in G.nodes()],
        "betweenness": [btw[x] for x in G.nodes()],
        "closeness": [clo[x] for x in G.nodes()],
        "eigenvector": [eig.get(x, np.nan) for x in G.nodes()],
    }).set_index("node")
    return df


# ----------------------------------------------------------------------------- #
# 4. Detecção de comunidades
# ----------------------------------------------------------------------------- #
def detect_louvain(G: nx.Graph, resolution: float = 1.0,
                   seed: int = 42) -> tuple[dict, float]:
    """Louvain (implementação do networkx, ciente de pesos).
    Retorna (partição {nó: id_comunidade}, modularidade)."""
    comms = nx.community.louvain_communities(
        G, weight="weight", resolution=resolution, seed=seed)
    comms = sorted(comms, key=len, reverse=True)
    part = {node: cid for cid, com in enumerate(comms) for node in com}
    Q = nx.community.modularity(G, comms, weight="weight")
    return part, Q


def louvain_stability(G: nx.Graph, n_seeds: int = 20,
                      resolution: float = 1.0) -> dict:
    """Roda Louvain com várias sementes e mede estabilidade (nº de comunidades
    e modularidade). Útil porque Louvain é estocástico."""
    ks, Qs = [], []
    for s in range(n_seeds):
        part, Q = detect_louvain(G, resolution=resolution, seed=s)
        ks.append(len(set(part.values()))); Qs.append(Q)
    return {"k_mean": float(np.mean(ks)), "k_std": float(np.std(ks)),
            "k_mode": int(np.bincount(ks).argmax()),
            "Q_mean": float(np.mean(Qs)), "Q_std": float(np.std(Qs)),
            "ks": ks, "Qs": Qs}


# ----------------------------------------------------------------------------- #
# 5. Validação contra anotação (cadeias ou domínios conhecidos)
# ----------------------------------------------------------------------------- #
def validate_against(part: dict, truth: dict) -> dict:
    """Compara a partição detectada com um rótulo verdadeiro (ex.: cadeia).
    Retorna ARI, NMI, homogeneidade, completude e a matriz de confusão."""
    from sklearn.metrics import (adjusted_rand_score,
                                  normalized_mutual_info_score,
                                  homogeneity_completeness_v_measure)
    import pandas as pd
    nodes = [x for x in part if x in truth]
    y_pred = [part[x] for x in nodes]
    y_true = [truth[x] for x in nodes]
    ari = adjusted_rand_score(y_true, y_pred)
    nmi = normalized_mutual_info_score(y_true, y_pred)
    hom, comp, v = homogeneity_completeness_v_measure(y_true, y_pred)
    cm = pd.crosstab(pd.Series(y_true, name="anotação"),
                     pd.Series(y_pred, name="comunidade"))
    return {"ARI": ari, "NMI": nmi, "homogeneity": hom,
            "completeness": comp, "v_measure": v, "confusion": cm}


def chain_truth(G: nx.Graph) -> dict:
    """Rótulo verdadeiro = cadeia do PDB (para multímeros)."""
    return {x: G.nodes[x]["chain"] for x in G.nodes()}


def load_6b1t_graph(cif_path: str = "data/6B1T.cif",
                    pkl_path: str = "data/6B1T_graph.pkl",
                    mode: str = "heavy", cutoff: float = 5.0) -> nx.Graph:
    """Carrega o grafo de contatos do 6B1T do cache (pickle); se ausente,
    constrói a partir do mmCIF e salva. Torna os scripts do alvo reprodutíveis
    a partir apenas de data/6B1T.cif."""
    import os, pickle
    if os.path.exists(pkl_path):
        return pickle.load(open(pkl_path, "rb"))
    if not os.path.exists(cif_path):
        raise FileNotFoundError(
            f"Nem {pkl_path} nem {cif_path} encontrados. Baixe 6B1T.cif do RCSB "
            f"(https://files.rcsb.org/download/6B1T.cif) e coloque em data/.")
    structure = load_structure(cif_path, "6B1T")
    G = giant_component(build_residue_network(structure, mode=mode, cutoff=cutoff))
    os.makedirs(os.path.dirname(pkl_path) or ".", exist_ok=True)
    pickle.dump(G, open(pkl_path, "wb"))
    return G
