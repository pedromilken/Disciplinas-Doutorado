# Migração do repositório — versão revisada

Resumo das mudanças após a devolutiva do orientador (parâmetros fixados nos
testes, ablação de critérios de aresta, verificação contra o RCSB). O número
de manchete passou de **ARI 0,89 (gamma escolhido no alvo)** para **ARI 0,812
(gamma 0,2 herdado da hemoglobina)**, removendo o vazamento de parâmetro.

---

## 1. ADICIONAR (arquivos novos)

| Arquivo | O que é |
|---|---|
| `edge_ablation.py` | Ablação dos critérios de aresta (topo / contacts / physchem / ca). |
| `derive_trimers_g02.py` | Re-derivação dos trímeros com gamma fixo em 0,2 (e 0,1 como sensibilidade). |
| `6B1T_trimer_g02.json` | Resultado da validação por capsômero em gamma 0,2 e 0,1. |
| `ablation_results.json` | Saída da ablação de arestas. |
| `figs/6B1T_net3d_g02.png` | Comunidades em gamma 0,2 sobre a estrutura 3D. |
| `figs/6B1T_confusion_trimer_g02.png` | Confusão comunidade x trímero/pénton/menores em gamma 0,2. |
| `figs/6B1T_confusion_macro_g02.png` | Confusão macro (3 classes) em gamma 0,2 (secundária). |

---

## 2. SUBSTITUIR (corrigidos)

| Arquivo | Mudança |
|---|---|
| `run_6b1t.py` | Removida a escolha de gamma por argmax da métrica no próprio alvo (linhas 74 a 80). Agora `g_macro = g_chain = 0.2`, fixado na calibração da hemoglobina. A varredura vira só sensibilidade. |

Substituir também os entregáveis de documento:

| Antigo | Novo |
|---|---|
| `relatorio_redes_complexas_proteinas.pdf` | `redes_complexas_relatorio_corrigido.docx` |
| `slides_redes_complexas_proteinas (1).pdf` | `redes_complexas_slides_corrigido.pptx` |

---

## 3. DESCARTAR (obsoletos)

| Arquivo | Motivo |
|---|---|
| `derive_trimers_6b1t.py` | Substituído por `derive_trimers_g02.py` (que roda 0,2 e 0,1). Se preferir manter, troque a linha 29 de `resolution=0.1` para `resolution=0.2`. |
| `figs/6B1T_net3d_g01.png` | Figura em gamma 0,1; superada por `6B1T_net3d_g02.png`. |
| `figs/6B1T_confusion_trimer.png` | Confusão de trímero em gamma 0,1; superada pela versão g02. |
| `figs/6B1T_confusion_macro.png` | Macro em gamma 0,1; superada por `6B1T_confusion_macro_g02.png`. |
| `6B1T_trimer.json` | Resultado em gamma 0,1 como manchete; rebaixado a sensibilidade dentro de `6B1T_trimer_g02.json`. Pode remover. |

> Não descartar `6B1T_net3d_g1.png` (gamma 1, escala fina, ainda citada na topologia) nem `figs/6B1T_confusion_official.png` (verificação RCSB).

---

## 4. ATUALIZAR (referências cruzadas)

- `README.md`: trocar a chamada `python derive_trimers_6b1t.py` por
  `python derive_trimers_g02.py`; acrescentar `python edge_ablation.py ...`.
- Qualquer pipeline (`run_pipeline.py`) que importe `derive_trimers_6b1t`
  deve apontar para o novo script.

---

## 5. Ordem de execução (reprodução)

```bash
python run_experiments.py          # testes (sintético, 1A8O, 2HHB, 2BEG)
python run_6b1t.py                 # topologia + figuras do alvo (gamma fixo 0,2)
python derive_trimers_g02.py       # ARI de capsômero: 0,812 (g=0,2) | 0,890 (g=0,1)
python edge_ablation.py --pdb data/2HHB.pdb --truth chain --gamma 0.2
python edge_ablation.py --pdb data/6B1T.cif --truth data/6B1T_trimer_truth.json --gamma 0.2
python fetch_rcsb_metadata.py 6B1T # anotação oficial (verificação)
```

## 6. Nota sobre a validação "macro" (3 classes)

Em gamma 0,2 a validação macro (héxon / pénton / menor) tem ARI baixo (~0,05),
porque os héxons se separam corretamente em 4 trímeros nessa escala. O número de
manchete é a validação por **trímero/capsômero** (`derive_trimers_g02.py`), não
a macro. A figura macro entra apenas como material secundário.
