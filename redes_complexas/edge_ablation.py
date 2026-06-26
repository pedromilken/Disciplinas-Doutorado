"""
edge_ablation.py
================
Ablação dos CRITÉRIOS DE ARESTA (resposta ao apontamento do orientador:
"a análise dos critérios de formação de aresta merece mais profundidade").

Em vez de só *argumentar* que a distância basta, aqui a gente MEDE. Fixados o
mesmo cutoff e o mesmo conjunto de nós, comparamos quatro esquemas de aresta/peso
sobre a(s) molécula(s)-TESTE e medimos quanto a partição muda (ARI vs anotação e
ARI cruzado entre esquemas). Se mudar pouco, fica provado que a distância já
recupera a topologia; se mudar muito, achamos algo digno de nota.

Esquemas comparados (cutoff e nós idênticos, só o peso/critério muda):
  1. topo      : aresta unitária (presença de contato), peso = 1
  2. contacts  : peso = nº de pares de átomos pesados em contato      [BASELINE atual]
  3. physchem  : peso = contacts * (1 + alpha * afinidade(i,j))
                 afinidade = atração hidrofóbica (KD) + complementaridade de carga
  4. ca8       : conjunto de ARESTAS diferente (modo Cα a 8 A), peso = 1/d
                 -> testa o efeito da DEFINIÇÃO de aresta, não só do peso

Uso:
  python edge_ablation.py --pdb data/2HHB.pdb --truth chain --gamma 0.2
  python edge_ablation.py --pdb data/SYNTH.pdb --truth data/SYNTH_truth.json \
                          --base-mode ca --base-cutoff 10 --gamma 1.0

Observação metodológica: gamma e cutoff são FIXADOS (não varridos no alvo).
Os valores default seguem a calibração nas moléculas-teste (ver relatório, Seção 2.x).
"""
from __future__ import annotations
import argparse, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import networkx as nx
import protein_net as pn

# --- propriedades físico-químicas por resíduo (3 letras) -------------------- #
# Hidrofobicidade de Kyte-Doolittle (positivo = hidrofóbico).
KD = {"ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9,
      "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8, "TRP": -0.9,
      "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5, "GLN": -3.5,
      "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5,
      "MSE": 1.9}  # MSE (selenometionina) ~ MET
# Carga aproximada em pH ~7.
CHG = {"ASP": -1.0, "GLU": -1.0, "LYS": 1.0, "ARG": 1.0, "HIS": 0.1}


def affinity(ri: str, rj: str) -> float:
    """Score de afinidade físico-química em ~[-0.5, 1].
    Positivo = interação favorável (hidrofóbico-hidrofóbico ou cargas opostas)."""
    hi, hj = KD.get(ri, 0.0), KD.get(rj, 0.0)
    qi, qj = CHG.get(ri, 0.0), CHG.get(rj, 0.0)
    hydro = (min(hi, hj) / 4.5) if (hi > 0 and hj > 0) else 0.0   # 0..1
    if qi * qj < 0:      salt = 1.0     # cargas opostas: ponte salina
    elif qi * qj > 0:    salt = -0.5    # cargas iguais: repulsão
    else:                salt = 0.0
    return 0.6 * hydro + 0.4 * salt


def reweight(G: nx.Graph, scheme: str, alpha: float = 1.0) -> nx.Graph:
    """Devolve uma CÓPIA de G com 'weight' segundo o esquema (mesmos nós/arestas)."""
    H = G.copy()
    for u, v, d in H.edges(data=True):
        nc = d.get("n_contacts", 1)
        if scheme == "topo":
            w = 1.0
        elif scheme == "contacts":
            w = float(nc)
        elif scheme == "physchem":
            ri = H.nodes[u].get("resname", ""); rj = H.nodes[v].get("resname", "")
            w = float(nc) * (1.0 + alpha * affinity(ri, rj))
            w = max(w, 0.05)               # evita peso <= 0
        else:
            raise ValueError(scheme)
        d["weight"] = w
    return H


def load_truth(spec: str, G: nx.Graph) -> dict:
    """truth = 'chain' (rótulo = cadeia) ou caminho de um JSON {nó: rótulo}."""
    if spec == "chain":
        return pn.chain_truth(G)
    truth = json.load(open(spec))
    return {x: truth[x] for x in G.nodes() if x in truth}


def evaluate(G: nx.Graph, truth: dict, gamma: float, seed: int = 42) -> dict:
    part, Q = pn.detect_louvain(G, resolution=gamma, seed=seed)
    val = pn.validate_against(part, truth)
    return {"k": len(set(part.values())), "Q": Q,
            "ARI": val["ARI"], "NMI": val["NMI"], "part": part}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--truth", default="chain",
                    help="'chain' ou caminho de JSON {nó: rótulo}")
    ap.add_argument("--gamma", type=float, default=0.2,
                    help="resolução FIXADA na calibração (não varrer no alvo)")
    ap.add_argument("--base-mode", default="heavy", choices=["heavy", "ca"])
    ap.add_argument("--base-cutoff", type=float, default=5.0)
    ap.add_argument("--alpha", type=float, default=1.0)
    args = ap.parse_args()

    st = pn.load_structure(args.pdb)
    base = pn.giant_component(
        pn.build_residue_network(st, mode=args.base_mode, cutoff=args.base_cutoff))
    truth = load_truth(args.truth, base)

    # conjunto de aresta alternativo: Cα a 8 A (testa DEFINIÇÃO de aresta)
    ca = pn.giant_component(pn.build_residue_network(st, mode="ca", cutoff=8.0))

    print(f"\n=== Ablação de arestas | {args.pdb} | gamma={args.gamma} (fixo) ===")
    print(f"rede base: {args.base_mode}@{args.base_cutoff}A  "
          f"n={base.number_of_nodes()} m={base.number_of_edges()}\n")
    header = f"{'esquema':<10} {'k':>3} {'Q':>7} {'ARI':>7} {'NMI':>7}  obs"
    print(header); print("-" * len(header))

    results = {}
    for scheme in ("topo", "contacts", "physchem"):
        G = reweight(base, scheme, alpha=args.alpha)
        r = evaluate(G, truth, args.gamma)
        results[scheme] = r
        obs = "baseline atual" if scheme == "contacts" else ""
        print(f"{scheme:<10} {r['k']:>3} {r['Q']:>7.4f} {r['ARI']:>7.4f} "
              f"{r['NMI']:>7.4f}  {obs}")

    # esquema com conjunto de aresta diferente (ca@8)
    truth_ca = load_truth(args.truth, ca)
    rca = evaluate(reweight(ca, "topo"), truth_ca, args.gamma)
    results["ca8"] = rca
    print(f"{'ca8':<10} {rca['k']:>3} {rca['Q']:>7.4f} {rca['ARI']:>7.4f} "
          f"{rca['NMI']:>7.4f}  conjunto de aresta alternativo")

    # ARI cruzado: contacts vs physchem (mudou a partição?)
    from sklearn.metrics import adjusted_rand_score
    nodes = list(base.nodes())
    pc = [results["contacts"]["part"][n] for n in nodes]
    pp = [results["physchem"]["part"][n] for n in nodes]
    cross = adjusted_rand_score(pc, pp)
    print(f"\nARI(contacts, physchem) = {cross:.4f}  "
          f"-> {'partição praticamente idêntica' if cross > 0.9 else 'partição mudou'}")
    if cross > 0.9:
        print("Leitura: a topologia da distância já fixa as comunidades; o peso "
              "físico-químico não altera o resultado.\n")
    else:
        print("Leitura: o peso físico-químico altera a partição. Compare os ARIs "
              "vs anotação acima para ver se ajuda ou atrapalha a recuperação.\n")

    out = {k: {kk: vv for kk, vv in v.items() if kk != "part"}
           for k, v in results.items()}
    out["cross_ARI_contacts_physchem"] = cross
    json.dump(out, open("ablation_results.json", "w"), indent=2)
    print("-> ablation_results.json")


if __name__ == "__main__":
    main()
