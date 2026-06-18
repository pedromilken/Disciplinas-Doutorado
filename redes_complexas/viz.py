"""
viz.py
======
Funções de visualização para o relatório. Salva PNGs em alta resolução.
Paleta alinhada à identidade visual: Ocean Deep, Royal Purple, gold/burnt-orange.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

OCEAN = "#1F4E79"; PURPLE = "#5B2A86"; GOLD = "#B35900"
plt.rcParams.update({"figure.dpi": 140, "font.size": 11,
                     "axes.titlesize": 12, "axes.titleweight": "bold",
                     "savefig.bbox": "tight"})


def _palette(k):
    base = [OCEAN, PURPLE, GOLD, "#2E8B57", "#C71585", "#1f77b4",
            "#ff7f0e", "#9467bd", "#8c564b", "#17becf"]
    if k <= len(base):
        return base[:k]
    return [cm.tab20(i / k) for i in range(k)]


def plot_degree_distribution(G, path, title=""):
    degs = np.array([d for _, d in G.degree()])
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].hist(degs, bins=min(30, max(5, degs.max())), color=OCEAN,
               edgecolor="white")
    ax[0].set(xlabel="grau k", ylabel="nº de resíduos",
              title=f"Distribuição de graus{title}")
    ax[0].axvline(degs.mean(), color=GOLD, ls="--", lw=2,
                  label=f"média = {degs.mean():.1f}")
    ax[0].legend()
    vals, counts = np.unique(degs, return_counts=True)
    ax[1].loglog(vals, counts, "o", color=PURPLE, ms=6)
    ax[1].set(xlabel="grau k (log)", ylabel="P(k) contagem (log)",
              title="Escala log–log")
    ax[1].grid(True, which="both", alpha=0.3)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def plot_centrality_hists(df, path):
    cols = ["degree_centrality", "betweenness", "closeness", "eigenvector"]
    colors = [OCEAN, PURPLE, GOLD, "#2E8B57"]
    fig, axes = plt.subplots(2, 2, figsize=(9, 6))
    for ax, c, col in zip(axes.ravel(), cols, colors):
        vals = df[c].dropna().values
        ax.hist(vals, bins=25, color=col, edgecolor="white")
        ax.set(title=c, ylabel="nº resíduos")
    fig.suptitle("Distribuições de centralidade", fontweight="bold")
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def plot_network_2d(G, partition, path, title="", layout_seed=42):
    import networkx as nx
    k = len(set(partition.values()))
    cols = _palette(k)
    node_color = [cols[partition[n]] for n in G.nodes()]
    pos = nx.spring_layout(G, seed=layout_seed, weight="weight",
                           k=1.5 / np.sqrt(max(1, G.number_of_nodes())))
    fig, ax = plt.subplots(figsize=(7.5, 7))
    nx.draw_networkx_edges(G, pos, alpha=0.12, width=0.5, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color=node_color, node_size=45,
                           linewidths=0.3, edgecolors="white", ax=ax)
    ax.set_title(title); ax.axis("off")
    fig.savefig(path); plt.close(fig)


def plot_network_3d(G, partition, path, title=""):
    """Projeção 3D usando coordenadas reais dos Cα, coloridas por comunidade."""
    k = len(set(partition.values()))
    cols = _palette(k)
    xs = np.array([G.nodes[n]["x"] for n in G.nodes()])
    ys = np.array([G.nodes[n]["y"] for n in G.nodes()])
    zs = np.array([G.nodes[n]["z"] for n in G.nodes()])
    c = [cols[partition[n]] for n in G.nodes()]
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(xs, ys, zs, c=c, s=30, depthshade=True, edgecolors="white",
               linewidths=0.2)
    ax.set_title(title)
    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)"); ax.set_zlabel("z (Å)")
    fig.savefig(path); plt.close(fig)


def plot_stability(stab, path):
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    ax[0].hist(stab["ks"], bins=range(min(stab["ks"]), max(stab["ks"]) + 2),
               color=OCEAN, edgecolor="white", align="left")
    ax[0].set(xlabel="nº de comunidades k", ylabel="frequência",
              title=f"Estabilidade de k ({len(stab['ks'])} sementes)")
    ax[1].hist(stab["Qs"], bins=15, color=PURPLE, edgecolor="white")
    ax[1].set(xlabel="modularidade Q", ylabel="frequência",
              title=f"Q = {stab['Q_mean']:.3f} ± {stab['Q_std']:.3f}")
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def plot_confusion(cm, path, title="Comunidade × anotação"):
    fig, ax = plt.subplots(figsize=(1.2 + 0.6 * cm.shape[1],
                                    1.2 + 0.5 * cm.shape[0]))
    im = ax.imshow(cm.values, cmap="BuPu", aspect="auto")
    ax.set_xticks(range(cm.shape[1])); ax.set_xticklabels(cm.columns)
    ax.set_yticks(range(cm.shape[0])); ax.set_yticklabels(cm.index)
    ax.set_xlabel("comunidade detectada"); ax.set_ylabel("anotação (cadeia)")
    ax.set_title(title)
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            v = cm.values[i, j]
            ax.text(j, i, int(v), ha="center", va="center",
                    color="white" if v > cm.values.max() / 2 else "black",
                    fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.8); fig.tight_layout()
    fig.savefig(path); plt.close(fig)
