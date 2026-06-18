# Declaração de uso de Inteligência Artificial

Conforme exigido no enunciado do projeto.

- **Ferramenta/Modelo:** Claude Opus 4.8 (Anthropic), acessado via interface de
  chat (claude.ai).
- **Temperatura:** a interface de chat não expõe nem permite ajustar a
  temperatura de amostragem (valor padrão do produto). Não foi feita chamada de
  API com temperatura/seed definidos pelo usuário.
- **Escopo do uso:**
  - Geração e organização do código do pipeline (`protein_net.py`, `viz.py`,
    `run_pipeline.py`, `run_experiments.py`, `make_synthetic.py`, `run_6b1t.py`,
    `derive_trimers_6b1t.py`).
  - Redação e formatação do relatório.
  - Sugestão de modelagem (rede de contatos por resíduo, esquema de pesos),
    métricas e estrutura experimental.
- **Verificação humana:** todos os resultados numéricos e figuras foram
  produzidos pela execução real do código sobre as estruturas reais, incluindo o
  arquivo `6B1T.cif` fornecido pelo autor. Não há valores inventados. O autor
  revisou, executou e validou o pipeline.
- **Prompt (resumo):** "Elaborar relatório completo e código reproduzível para
  o projeto de redes complexas de detecção de domínios funcionais em proteínas
  (6B1T), com modelagem de rede justificada, análise topológica, detecção de
  comunidades (Louvain) e validação contra anotações." Transcrição completa
  disponível mediante solicitação.
