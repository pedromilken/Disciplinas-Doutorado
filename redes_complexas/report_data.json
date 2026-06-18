# Detecção de domínios funcionais em proteínas via redes complexas

Pipeline reprodutível que modela uma proteína como **rede de contatos de
resíduos**, analisa sua **topologia** (graus e centralidades) e detecta
**comunidades** (Louvain), validando-as contra anotações estruturais conhecidas.

Projeto da disciplina de Redes Complexas — Prof. Gilberto Nakamura (UFU).
Uso de IA declarado em [`declaracao_ia.md`](declaracao_ia.md).

## Ideia em uma frase

Uma proteína é uma "rede social de aminoácidos": cada resíduo é uma pessoa, e há
um laço entre dois que ficam fisicamente próximos. Os **domínios funcionais** são
as **comunidades** dessa rede — grupos densamente conectados por dentro e
esparsamente ligados entre si.

## Instalação

```bash
python -m venv .venv && source .venv/bin/activate   # opcional
pip install -r requirements.txt
```

## Estrutura

| Arquivo | Papel |
|---|---|
| `protein_net.py` | Biblioteca: fetch/parse PDB, rede de contatos, topologia, Louvain, validação |
| `viz.py` | Figuras (graus, centralidades, rede 2D/3D, estabilidade, matriz de confusão) |
| `make_synthetic.py` | Gera proteína sintética com domínios plantados (*ground truth*) |
| `run_experiments.py` | Reproduz todos os casos do relatório (SYNTH, 1A8O, 2BEG, 2HHB) |
| `run_pipeline.py` | CLI geral para rodar qualquer proteína (inclui o alvo 6B1T) |

## Reproduzir os resultados do relatório

```bash
python run_experiments.py        # gera figs/ e results.json
```

## Rodar em uma proteína qualquer

```bash
# proteína-teste pequena
python run_pipeline.py --pdb 1A8O --mode heavy --cutoff 5.0

# multímero, validando contra as cadeias, resolução ajustada
python run_pipeline.py --pdb 2HHB --mode heavy --cutoff 5.0 \
       --resolution 0.2 --validate-chains
```

## Rodar o ALVO 6B1T (na sua máquina, com acesso ao rcsb.org)

O 6B1T é um capsídeo grande (~10 GB se carregado por átomo). A rede é construída
**por resíduo** com KD-tree, o que reduz drasticamente o custo. Recomenda-se
`--mode ca --cutoff 8.0` e resolução baixa para separar penton/hexon:

```bash
python run_pipeline.py --pdb 6B1T --mode ca --cutoff 8.0 \
       --resolution 0.5 --validate-chains --outdir out_6b1t
```

Para validar contra as famílias anotadas (em vez das cadeias), passe um CSV
`node,label` (com `node = "CADEIA:RESNAME:RESSEQ"`):

```bash
python run_pipeline.py --pdb 6B1T --annotation familias_6b1t.csv --outdir out_6b1t
```

A partição é exportada em `out_6b1t/6B1T_partition.json`; use-a para colorir os
clusters no visualizador 3D do RCSB (https://www.rcsb.org/3d-view/6B1T).

### Alvo 6B1T — já executado neste trabalho

O alvo foi rodado sobre a unidade assimétrica (`6B1T.cif`, ~100 mil átomos →
12.544 resíduos). Para reproduzir exatamente os resultados e figuras do relatório:

```bash
# coloque 6B1T.cif em data/, então:
python run_6b1t.py              # topologia, varredura de γ, Louvain, figuras
python derive_trimers_6b1t.py   # anotação por trímero + validação (ARI ≈ 0,89)
```

Resultado: Q = 0,90 (γ=1); em γ=0,1 as comunidades reproduzem os 4 trímeros de
héxon + base do pénton + proteínas menores, com **ARI = 0,89** (NMI = 0,86).
Saídas em `data/6B1T_results.json` e `data/6B1T_trimer.json`.

### Anotação oficial via API do RCSB (opcional)

`fetch_rcsb_metadata.py` consulta a API de dados do RCSB (GraphQL, em
`data.rcsb.org`) e gera `data/6B1T_annotation.csv` (cadeia → nome da proteína via
UniProt) e a simetria do assembly. A query usada está em `rcsb_query_full.graphql`.
Requer acesso a `data.rcsb.org` (rode na sua máquina):

```bash
python fetch_rcsb_metadata.py 6B1T
python run_pipeline.py --pdb 6B1T --annotation data/6B1T_annotation.csv
```

A validação contra a anotação oficial wwPDB (8 famílias) é reproduzida por
`validate_official_6b1t.py`. O `.cif` e o PDB bundle produzem redes idênticas
(12.544 nós, 73.276 arestas).

## Notas de desempenho

- Rede **por resíduo** (não por átomo) → de ~10⁶ nós para ~10⁴–10⁵.
- Busca de vizinhos com **KD-tree** (`Bio.PDB.NeighborSearch`).
- Betweenness usa **amostragem por pivôs** automaticamente em grafos grandes.
- Caminho médio/diâmetro só são calculados se o componente gigante for ≤ 5000 nós.
