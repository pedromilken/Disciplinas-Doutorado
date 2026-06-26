"""
derive_trimers_g02.py
=====================
Re-derivação da validação por TRÍMERO do 6B1T com gamma FIXADO na calibração
das moléculas-teste, em resposta ao apontamento do orientador:
"os parâmetros devem ser fixados nos testes com as outras moléculas".

Mudança em relação a derive_trimers_6b1t.py:
  - Antes: gamma escolhido como argmax do ARI macro NO PRÓPRIO 6B1T (vazamento).
  - Agora: gamma = 0.2, FIXADO a partir da hemoglobina (2HHB), onde gamma=0.2
           recupera as cadeias montadas (subunidades) com ARI maximo (0.839).
           A escala "subunidade montada" da hemoglobina transfere para a escala
           "capsômero" do 6B1T. Aplicamos gamma=0.2 ao alvo SEM olhar a anotação
           do alvo.

Roda gamma=0.2 (resultado-manchete) e gamma=0.1 (sensibilidade) lado a lado.
A anotação por trímero é a mesma de antes; o trímero permanece definido pelos
contatos entre cadeias de héxon (rótulo grosseiro e essencialmente determinístico;
ver discussão de circularidade no relatório).

Requer data/6B1T.cif (ou o cache data/6B1T_graph.pkl). Não roda no sandbox
restrito; rode na sua máquina.
"""
import warnings; warnings.filterwarnings("ignore")
import json, networkx as nx
import protein_net as pn

GAMMA_FIX = 0.2            # fixado na hemoglobina (subunidade montada)
GAMMA_SENS = 0.1          # sensibilidade (valor antigo)
HEXON = set("ABCDEFGHIJKL")


def trimer_truth(G: nx.Graph) -> dict:
    """Rótulo verdadeiro: 4 trímeros de héxon + pénton (M) + menores.
    Trímeros = comunidades no grafo de contatos ENTRE cadeias de héxon
    (independente da partição de resíduos que vamos validar)."""
    cc = {}
    for u, v in G.edges():
        a, b = G.nodes[u]["chain"], G.nodes[v]["chain"]
        if a == b:
            continue
        e = tuple(sorted((a, b)))
        cc[e] = cc.get(e, 0) + 1
    Cg = nx.Graph()
    for (a, b), w in cc.items():
        Cg.add_edge(a, b, weight=w)
    Hsub = Cg.subgraph([c for c in Cg if c in HEXON]).copy()
    tri = sorted(nx.community.louvain_communities(Hsub, weight="weight", seed=1),
                 key=lambda s: sorted(s))
    lab = {}
    for i, s in enumerate(tri):
        for c in s:
            lab[c] = f"hexon_T{i}"
    lab["M"] = "penton"
    truth = {x: lab.get(G.nodes[x]["chain"], "minor") for x in G.nodes()}
    return truth, [sorted(s) for s in tri]


def run(G, truth, gamma):
    part = pn.detect_louvain(G, resolution=gamma, seed=42)[0]
    val = pn.validate_against(part, truth)
    return {"gamma": gamma, "k": len(set(part.values())),
            "ARI": round(val["ARI"], 4), "NMI": round(val["NMI"], 4),
            "homogeneity": round(val["homogeneity"], 4),
            "completeness": round(val["completeness"], 4),
            "confusion": val["confusion"].to_dict()}


def main():
    G = pn.load_6b1t_graph()
    truth, tri = trimer_truth(G)
    print("trímeros de héxon (definidos por contato entre cadeias):", tri)

    res = {}
    print(f"\n{'gamma':>6} {'k':>3} {'ARI':>7} {'NMI':>7} {'homog':>7} {'compl':>7}  papel")
    for g, role in [(GAMMA_FIX, "MANCHETE (fixado na hemoglobina)"),
                    (GAMMA_SENS, "sensibilidade (valor antigo)")]:
        r = run(G, truth, g)
        res[f"g{g}"] = r
        print(f"{g:>6} {r['k']:>3} {r['ARI']:>7.4f} {r['NMI']:>7.4f} "
              f"{r['homogeneity']:>7.4f} {r['completeness']:>7.4f}  {role}")

    json.dump(res, open("6B1T_trimer_g02.json", "w"), indent=2)
    print("\n-> 6B1T_trimer_g02.json")
    print("Use o resultado de gamma=0.2 como número-manchete; gamma=0.1 vira "
          "linha de sensibilidade. Se os dois baterem (ambos ~0.88-0.89), o "
          "resultado é robusto à escolha e a proveniência fica limpa.")


if __name__ == "__main__":
    main()
