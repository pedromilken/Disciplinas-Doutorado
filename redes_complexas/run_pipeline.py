"""
run_experiments.py
==================
Executa o pipeline completo nos casos do relatório e salva figuras + um JSON
de resultados (results.json) que alimenta o relatório.

Casos:
  (A) SYNTH  - benchmark sintético com ground truth plantado (4 domínios)
  (B) 1A8O   - proteína-teste real, pequena, 1 cadeia (validação do pipeline)
  (C) 2BEG   - multímero real (5 cadeias) -> validação biológica (comunidade=cadeia)
"""
import os, json, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
import protein_net as pn
import viz
from make_synthetic import make_synthetic_pdb

os.makedirs("figs", exist_ok=True); os.makedirs("data", exist_ok=True)
RESULTS = {}


def run_case(name, pdb_path, mode, cutoff, truth=None, n_seeds=20):
    t0 = time.time()
    structure = pn.load_structure(pdb_path, name)
    G = pn.build_residue_network(structure, mode=mode, cutoff=cutoff)
    G = pn.giant_component(G)        # foca no componente principal
    topo = pn.topology_summary(G)
    df = pn.centralities(G)
    part, Q = pn.detect_louvain(G)
    stab = pn.louvain_stability(G, n_seeds=n_seeds)
    k = len(set(part.values()))
    elapsed = time.time() - t0

    # figuras
    viz.plot_degree_distribution(G, f"figs/{name}_degree.png", f" — {name}")
    viz.plot_centrality_hists(df, f"figs/{name}_centrality.png")
    viz.plot_network_2d(G, part, f"figs/{name}_net2d.png",
                        f"{name}: rede colorida por comunidade (k={k})")
    viz.plot_network_3d(G, part, f"figs/{name}_net3d.png",
                        f"{name}: comunidades sobre a estrutura 3D")
    viz.plot_stability(stab, f"figs/{name}_stability.png")

    # top hubs por betweenness (resíduos "ponte" entre domínios)
    top_btw = df.sort_values("betweenness", ascending=False).head(5)
    top_btw_list = [(idx, round(r.betweenness, 4), int(r.degree))
                    for idx, r in top_btw.iterrows()]

    res = {"name": name, "mode": mode, "cutoff": cutoff,
           "topology": {kk: (round(v, 4) if isinstance(v, float) else v)
                        for kk, v in topo.items()},
           "modularity": round(Q, 4), "n_communities": k,
           "stability": {kk: (round(v, 4) if isinstance(v, float) else v)
                         for kk, v in stab.items() if kk not in ("ks", "Qs")},
           "elapsed_s": round(elapsed, 2),
           "top_betweenness": top_btw_list}

    if truth is not None:
        val = pn.validate_against(part, truth)
        viz.plot_confusion(val["confusion"], f"figs/{name}_confusion.png")
        res["validation"] = {kk: round(v, 4) for kk, v in val.items()
                             if kk != "confusion"}
        res["confusion"] = val["confusion"].to_dict()
    RESULTS[name] = res
    print(f"[{name}] n={topo['n_nodes']} m={topo['n_edges']} "
          f"k={k} Q={Q:.3f} t={elapsed:.1f}s "
          f"{'ARI=%.3f' % res['validation']['ARI'] if truth else ''}")
    return res


# (A) sintético -------------------------------------------------------------
truth = make_synthetic_pdb("data/SYNTH.pdb", n_domains=4, per_domain=60)
json.dump(truth, open("data/SYNTH_truth.json", "w"))
run_case("SYNTH", "data/SYNTH.pdb", mode="ca", cutoff=10.0, truth=truth)

# (B) 1A8O ------------------------------------------------------------------
run_case("1A8O", "data/1A8O.pdb", mode="heavy", cutoff=5.0)

# (C) 2BEG (multímero) ------------------------------------------------------
st = pn.load_structure("data/2BEG.pdb", "2BEG")
G2 = pn.giant_component(pn.build_residue_network(st, mode="heavy", cutoff=5.0))
chain_t = pn.chain_truth(G2)
run_case("2BEG", "data/2BEG.pdb", mode="heavy", cutoff=5.0, truth=chain_t)

# (D) 2HHB (hemoglobina, multímero globular 2alfa+2beta) -------------------
st = pn.load_structure("data/2HHB.pdb", "2HHB")
G4 = pn.giant_component(pn.build_residue_network(st, mode="heavy", cutoff=5.0))
chain_t4 = pn.chain_truth(G4)
run_case("2HHB", "data/2HHB.pdb", mode="heavy", cutoff=5.0, truth=chain_t4)

json.dump(RESULTS, open("results.json", "w"), indent=2)
print("\nOK -> results.json e figs/*.png")
