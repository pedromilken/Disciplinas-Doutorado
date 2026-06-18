"""validate_official_6b1t.py — valida a partição do 6B1T contra a anotação
oficial de proteínas do relatório de validação wwPDB (8 tipos de molécula)."""
import warnings; warnings.filterwarnings("ignore")
import pickle, json
import protein_net as pn, viz

G = pickle.load(open("data/6B1T_graph.pkl", "rb"))
# mapeamento cadeia -> proteína (relatório wwPDB de 6B1T / EMD-7034)
official = {}
for c in "ABCDEFGHIJKL": official[c] = "Hexon"          # Molécula 1
official["M"] = "Penton"                                  # Molécula 2
official["N"] = "IIIa"                                    # Molécula 3
for c in "OP": official[c] = "VIII"                       # Molécula 4
for c in "QRST": official[c] = "IX"                       # Molécula 5 (interlacing)
for c in "UVY": official[c] = "VI"                        # Molécula 6
official["W"] = "VII"                                     # Molécula 7
official["X"] = "VI"                                      # Molécula 8
truth = {x: official[G.nodes[x]["chain"]] for x in G.nodes()}

part = pn.detect_louvain(G, resolution=0.1, seed=42)[0]
val = pn.validate_against(part, truth)
print(f"vs anotação oficial (8 famílias): ARI={val['ARI']:.3f} NMI={val['NMI']:.3f} "
      f"hom={val['homogeneity']:.3f} comp={val['completeness']:.3f}")
print(val["confusion"])
viz.plot_confusion(val["confusion"], "figs/6B1T_confusion_official.png",
                   "6B1T (gamma=0.1): comunidade x proteina (anotacao wwPDB)")
