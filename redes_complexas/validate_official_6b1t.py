#!/usr/bin/env python3
"""
run_pipeline.py — CLI geral do pipeline de detecção de domínios funcionais.

Exemplos
--------
# Proteína-teste pequena (baixa do RCSB):
python run_pipeline.py --pdb 1A8O --mode heavy --cutoff 5.0

# Multímero, validando contra as cadeias do PDB, com resolução ajustada:
python run_pipeline.py --pdb 2HHB --mode heavy --cutoff 5.0 \
       --resolution 0.2 --validate-chains

# ALVO 6B1T (rode na sua máquina, com acesso ao rcsb.org):
#  - use mode='ca' e cutoff 8.0 para caber em memória (rede por resíduo, não por átomo)
#  - resolução baixa para separar grandes unidades (penton/hexon)
python run_pipeline.py --pdb 6B1T --mode ca --cutoff 8.0 \
       --resolution 0.5 --validate-chains --outdir out_6b1t

# Validar contra uma anotação própria (CSV: node,label  -> node = 'CHAIN:RESNAME:RESSEQ'):
python run_pipeline.py --pdb 6B1T --annotation familias_6b1t.csv
"""
import argparse, os, json, time, warnings
warnings.filterwarnings("ignore")
import protein_net as pn
import viz


def main():
    ap = argparse.ArgumentParser(description="Detecção de domínios funcionais via redes de contato.")
    ap.add_argument("--pdb", required=True, help="PDB ID (baixa do RCSB) ou caminho .pdb/.cif")
    ap.add_argument("--url", default=None, help="URL alternativa do PDB (ex.: GitHub raw)")
    ap.add_argument("--mode", choices=["heavy", "ca"], default="heavy")
    ap.add_argument("--cutoff", type=float, default=None, help="limiar (Å). Padrão 5.0(heavy)/8.0(ca)")
    ap.add_argument("--resolution", type=float, default=1.0, help="resolução do Louvain (γ)")
    ap.add_argument("--seeds", type=int, default=20, help="nº de sementes p/ estabilidade")
    ap.add_argument("--validate-chains", action="store_true",
                    help="usa a cadeia do PDB como anotação verdadeira")
    ap.add_argument("--annotation", default=None,
                    help="CSV node,label com anotação verdadeira customizada")
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    figdir = os.path.join(args.outdir, "figs"); os.makedirs(figdir, exist_ok=True)
    name = os.path.splitext(os.path.basename(args.pdb))[0].upper()

    t0 = time.time()
    path = args.pdb if os.path.exists(args.pdb) else pn.fetch_pdb(args.pdb, args.outdir, url=args.url)
    print(f"[1/5] estrutura: {path}")
    structure = pn.load_structure(path, name)

    print(f"[2/5] construindo rede ({args.mode}, cutoff={args.cutoff or 'padrão'}) ...")
    G = pn.build_residue_network(structure, mode=args.mode, cutoff=args.cutoff)
    G = pn.giant_component(G)

    print("[3/5] topologia ...")
    topo = pn.topology_summary(G)
    df = pn.centralities(G)
    df.to_csv(os.path.join(args.outdir, f"{name}_centralities.csv"))

    print(f"[4/5] Louvain (γ={args.resolution}) + estabilidade ...")
    part, Q = pn.detect_louvain(G, resolution=args.resolution)
    stab = pn.louvain_stability(G, n_seeds=args.seeds, resolution=args.resolution)
    k = len(set(part.values()))

    print("[5/5] figuras ...")
    viz.plot_degree_distribution(G, f"{figdir}/{name}_degree.png", f" — {name}")
    viz.plot_centrality_hists(df, f"{figdir}/{name}_centrality.png")
    viz.plot_network_2d(G, part, f"{figdir}/{name}_net2d.png", f"{name}: comunidades (k={k})")
    viz.plot_network_3d(G, part, f"{figdir}/{name}_net3d.png", f"{name}: comunidades 3D")
    viz.plot_stability(stab, f"{figdir}/{name}_stability.png")

    summary = {"name": name, "params": vars(args),
               "topology": topo, "modularity": Q, "n_communities": k,
               "stability": {kk: v for kk, v in stab.items() if kk not in ("ks", "Qs")},
               "elapsed_s": round(time.time() - t0, 2)}

    truth = None
    if args.annotation:
        import pandas as pd
        ann = pd.read_csv(args.annotation)
        truth = dict(zip(ann.iloc[:, 0], ann.iloc[:, 1]))
    elif args.validate_chains:
        truth = pn.chain_truth(G)
    if truth:
        val = pn.validate_against(part, truth)
        viz.plot_confusion(val["confusion"], f"{figdir}/{name}_confusion.png")
        summary["validation"] = {kk: v for kk, v in val.items() if kk != "confusion"}
        print(f"   ARI={val['ARI']:.3f}  NMI={val['NMI']:.3f}")

    # exporta partição (para colorir no visualizador 3D do RCSB, p.ex.)
    json.dump(part, open(os.path.join(args.outdir, f"{name}_partition.json"), "w"))
    json.dump({kk: (str(v) if not isinstance(v, (int, float, list, dict, type(None))) else v)
               for kk, v in summary.items()},
              open(os.path.join(args.outdir, f"{name}_summary.json"), "w"),
              indent=2, default=str)
    print(f"\nOK: n={topo['n_nodes']} m={topo['n_edges']} k={k} Q={Q:.3f} "
          f"em {summary['elapsed_s']}s -> {args.outdir}/")


if __name__ == "__main__":
    main()
