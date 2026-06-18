"""
make_synthetic.py
=================
Gera uma "proteína" sintética com K domínios espaciais plantados (ground truth
conhecido), escrita como PDB válido (apenas Cα). Serve de benchmark controlado:
sabemos a resposta certa, então medimos ARI/NMI do Louvain contra a verdade.

Analogia: é como desenhar de propósito 4 ilhas ligadas por pontes finas e depois
pedir ao algoritmo que redescubra as ilhas. Se ele acerta aqui, confiamos nele
no caso real, onde não conhecemos as ilhas de antemão.
"""
import numpy as np

AA = ["ALA", "GLY", "SER", "LEU", "VAL", "THR", "LYS", "GLU", "ASP", "ILE"]


def make_synthetic_pdb(path="data/SYNTH.pdb", n_domains=4, per_domain=60,
                       spread=6.0, sep=26.0, linker=5, seed=7):
    """Cria K blobs gaussianos (domínios) dispostos em torno de um anel, ligados
    por pequenos linkers sequenciais. Retorna o dict {residuo: dominio_verdadeiro}."""
    rng = np.random.default_rng(seed)
    centers = []
    for k in range(n_domains):
        ang = 2 * np.pi * k / n_domains
        centers.append(np.array([sep * np.cos(ang), sep * np.sin(ang), 0.0]))
    coords, truth, chain_of = [], {}, {}
    serial = 0
    resseq = 0
    for k in range(n_domains):
        for _ in range(per_domain):
            resseq += 1; serial += 1
            xyz = centers[k] + rng.normal(0, spread, size=3)
            coords.append((serial, resseq, xyz))
            truth[f"A:{AA[serial % len(AA)]}:{resseq}"] = k
        # linker para o próximo domínio (residuos intermediários)
        if k < n_domains - 1:
            a, b = centers[k], centers[k + 1]
            for t in range(1, linker + 1):
                resseq += 1; serial += 1
                xyz = a + (b - a) * t / (linker + 1) + rng.normal(0, 1.5, size=3)
                coords.append((serial, resseq, xyz))
                truth[f"A:{AA[serial % len(AA)]}:{resseq}"] = k  # linker conta p/ dom. k
    # escreve PDB
    with open(path, "w") as f:
        for serial, resseq, xyz in coords:
            res = AA[serial % len(AA)]
            f.write(
                f"ATOM  {serial:5d}  CA  {res} A{resseq:4d}    "
                f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00           C\n")
        f.write("END\n")
    return truth


if __name__ == "__main__":
    import os, json
    os.makedirs("data", exist_ok=True)
    truth = make_synthetic_pdb()
    json.dump(truth, open("data/SYNTH_truth.json", "w"))
    print(f"SYNTH.pdb gerado: {len(truth)} resíduos, "
          f"{len(set(truth.values()))} domínios verdadeiros.")
