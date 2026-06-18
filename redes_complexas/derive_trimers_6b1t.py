"""derive_trimers_6b1t.py — define a anotação por trímero do héxon a partir dos
contatos ENTRE cadeias (independente do Louvain) e valida a partição do alvo."""
import warnings; warnings.filterwarnings("ignore")
import pickle, json, networkx as nx
import protein_net as pn, viz

G = pickle.load(open("data/6B1T_graph.pkl", "rb"))
HEXON = set("ABCDEFGHIJKL")

# grafo de contatos entre cadeias (peso = nº de contatos de resíduo)
cc = {}
for u, v in G.edges():
    a, b = G.nodes[u]["chain"], G.nodes[v]["chain"]
    if a == b: continue
    e = tuple(sorted((a, b))); cc[e] = cc.get(e, 0) + 1
Cg = nx.Graph()
for (a, b), w in cc.items(): Cg.add_edge(a, b, weight=w)

# trímeros = comunidades no subgrafo das cadeias de héxon
Hsub = Cg.subgraph([c for c in Cg if c in HEXON]).copy()
tri = sorted(nx.community.louvain_communities(Hsub, weight="weight", seed=1),
             key=lambda s: sorted(s))
lab = {}
for i, s in enumerate(tri):
    for c in s: lab[c] = f"hexon_T{i}"
lab["M"] = "penton"
truth = {x: lab.get(G.nodes[x]["chain"], "minor") for x in G.nodes()}

part = pn.detect_louvain(G, resolution=0.1, seed=42)[0]
val = pn.validate_against(part, truth)
print("trímeros:", [sorted(s) for s in tri])
print(f"ARI={val['ARI']:.3f} NMI={val['NMI']:.3f} "
      f"hom={val['homogeneity']:.3f} comp={val['completeness']:.3f}")
print(val["confusion"])
viz.plot_confusion(val["confusion"], "figs/6B1T_confusion_trimer.png",
                   "6B1T (gamma=0.1): comunidade x (trimero/penton/menores)")
