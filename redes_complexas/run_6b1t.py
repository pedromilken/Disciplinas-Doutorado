"""run_6b1t.py — análise completa do alvo (capsídeo do adenovírus 5, PDB 6B1T)."""
import warnings; warnings.filterwarnings("ignore")
import time, json, pickle
import numpy as np
import networkx as nx
import protein_net as pn
import viz

G = pickle.load(open("data/6B1T_graph.pkl", "rb"))
N = G.number_of_nodes()
print(f"grafo: n={N} m={G.number_of_edges()}")

# ---- anotações ----
# nível 1: cadeia (25 rótulos)
chain_truth = {x: G.nodes[x]["chain"] for x in G.nodes()}
# nível 2: macro família estrutural
HEXON = set("ABCDEFGHIJKL")       # 12 cadeias ~900 res -> hexon (4 trímeros)
PENTON = {"M"}                     # base do penton
def macro(ch):
    if ch in HEXON: return "hexon"
    if ch in PENTON: return "penton"
    return "minor"
macro_truth = {x: macro(G.nodes[x]["chain"]) for x in G.nodes()}

# ---- topologia ----
t = time.time()
degs = np.array([d for _, d in G.degree()])
topo = {"n_nodes": N, "n_edges": G.number_of_edges(),
        "density": nx.density(G), "avg_degree": float(degs.mean()),
        "max_degree": int(degs.max()), "avg_clustering": nx.average_clustering(G),
        "n_components": nx.number_connected_components(G)}
print("topologia:", round(time.time() - t, 1), "s ->", {k: round(v,4) if isinstance(v,float) else v for k,v in topo.items()})

# ---- centralidades (eficiente p/ grafo grande) ----
t = time.time()
deg_c = nx.degree_centrality(G)
btw = nx.betweenness_centrality(G, k=600, weight="dist", normalized=True, seed=42)
try:
    eig = nx.eigenvector_centrality_numpy(G, weight="weight")
except Exception:
    eig = {x: np.nan for x in G.nodes()}
import pandas as pd
dfc = pd.DataFrame({
    "node": list(G.nodes()),
    "chain": [G.nodes[x]["chain"] for x in G.nodes()],
    "degree": [dict(G.degree())[x] for x in G.nodes()],
    "degree_centrality": [deg_c[x] for x in G.nodes()],
    "betweenness": [btw[x] for x in G.nodes()],
    "eigenvector": [eig.get(x, np.nan) for x in G.nodes()],
}).set_index("node")
dfc.to_csv("out_6b1t/6B1T_centralities.csv") if __import__("os").path.isdir("out_6b1t") else dfc.to_csv("6B1T_centralities.csv")
print("centralidades:", round(time.time() - t, 1), "s")

# ---- Louvain default (gamma=1) ----
t = time.time()
part1, Q1 = pn.detect_louvain(G, resolution=1.0)
k1 = len(set(part1.values()))
print(f"Louvain g=1: k={k1} Q={Q1:.3f} ({round(time.time()-t,1)}s)")

# ---- varredura de gamma vs anotações ----
rows = []
parts = {}
for g in [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]:
    p, Q = pn.detect_louvain(G, resolution=g)
    parts[g] = p
    vm = pn.validate_against(p, macro_truth)
    vc = pn.validate_against(p, chain_truth)
    k = len(set(p.values()))
    rows.append([g, k, round(Q,3),
                 round(vm["ARI"],3), round(vm["NMI"],3),
                 round(vc["ARI"],3), round(vc["NMI"],3)])
    print(f"  g={g:<4} k={k:<3} Q={Q:.3f} | macro ARI={vm['ARI']:.3f} NMI={vm['NMI']:.3f} | chain ARI={vc['ARI']:.3f} NMI={vc['NMI']:.3f}")

# escolhe gamma cujo k melhor casa com nº de cadeias (25) p/ validação de cadeia
best_chain = min(rows, key=lambda r: (abs(r[1]-25), -r[5]))
g_chain = best_chain[0]
# e gamma que melhor separa macro (3 famílias) -> ARI macro máx
best_macro = max(rows, key=lambda r: r[3])
g_macro = best_macro[0]
print("melhor g p/ cadeia:", g_chain, "| melhor g p/ macro:", g_macro)

# ---- estabilidade (em g=1) ----
t = time.time()
stab = pn.louvain_stability(G, n_seeds=15, resolution=1.0)
print("estabilidade:", round(time.time()-t,1),"s", {k:round(v,3) for k,v in stab.items() if k not in('ks','Qs')})

# ---- validações finais nos gammas escolhidos ----
val_chain = pn.validate_against(parts[g_chain], chain_truth)
val_macro = pn.validate_against(parts[g_macro], macro_truth)

# ---- figuras ----
import os; os.makedirs("figs", exist_ok=True)
viz.plot_degree_distribution(G, "figs/6B1T_degree.png", " — 6B1T")
viz.plot_centrality_hists(dfc.assign(closeness=np.nan), "figs/6B1T_centrality.png")
viz.plot_stability(stab, "figs/6B1T_stability.png")
viz.plot_network_3d(G, part1, "figs/6B1T_net3d_g1.png",
                    f"6B1T: comunidades (gamma=1, k={k1})")
viz.plot_network_3d(G, parts[g_chain], "figs/6B1T_net3d_chain.png",
                    f"6B1T: comunidades (gamma={g_chain}, k={len(set(parts[g_chain].values()))})")
viz.plot_confusion(val_macro["confusion"], "figs/6B1T_confusion_macro.png",
                   f"6B1T (gamma={g_macro}): comunidade x família")

# matriz comunidade x cadeia (heatmap grande)
viz.plot_confusion(val_chain["confusion"], "figs/6B1T_confusion_chain.png",
                   f"6B1T (gamma={g_chain}): comunidade x cadeia")

# top hubs por betweenness
top_btw = dfc.sort_values("betweenness", ascending=False).head(8)
top_btw_list = [(i, round(r.betweenness,4), int(r.degree), r.chain) for i,r in top_btw.iterrows()]

out = {
    "topology": {k:(round(v,4) if isinstance(v,float) else v) for k,v in topo.items()},
    "louvain_g1": {"k": k1, "Q": round(Q1,4)},
    "sweep": rows,                      # [g,k,Q,macroARI,macroNMI,chainARI,chainNMI]
    "g_chain": g_chain, "g_macro": g_macro,
    "val_chain": {k:round(v,4) for k,v in val_chain.items() if k!="confusion"},
    "val_macro": {k:round(v,4) for k,v in val_macro.items() if k!="confusion"},
    "stability": {k:round(v,4) for k,v in stab.items() if k not in('ks','Qs')},
    "confusion_macro": val_macro["confusion"].to_dict(),
    "top_betweenness": top_btw_list,
    "annotation_counts": {"hexon": sum(1 for v in macro_truth.values() if v=='hexon'),
                          "penton": sum(1 for v in macro_truth.values() if v=='penton'),
                          "minor": sum(1 for v in macro_truth.values() if v=='minor')},
}
json.dump(out, open("data/6B1T_results.json","w"), indent=2, default=str)
json.dump(part1, open("data/6B1T_partition_g1.json","w"))
print("\nOK -> data/6B1T_results.json")
print("val_macro:", out["val_macro"])
print("val_chain:", out["val_chain"])
print("confusion macro:\n", val_macro["confusion"])
