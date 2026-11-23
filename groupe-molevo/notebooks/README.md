# Notebook de analise

`analise.ipynb` deve ler `resultados/registro.jsonl` e `resultados/resumo.json`
e produzir as figuras. Ele **nao** importa nada de `src/ga.py` nem roda a
busca: refazer uma figura nao pode custar uma execucao inteira.

## Figuras que o relatorio precisa

1. **Validade sob mutacao** -- SMILES contra SELFIES, medido por voces.
   E a justificativa do projeto inteiro, em um grafico de barras.
2. **Convergencia** -- GA, aleatoria e hill climbing. Eixo x em
   **avaliacoes**, faixa de desvio entre sementes.
3. **Aptidao e diversidade** -- duas curvas, dois eixos y.
4. **Frente de Pareto** -- QED contra sa_norm, dominados em cinza.
5. **As melhores moleculas** -- grade desenhada com
   `rdkit.Chem.Draw.MolsToGridImage`, com QED e SA em cada legenda.
   Inclua tambem as moleculas *estranhas* de QED alto: sao a evidencia da
   secao de limites.

## Esqueleto

```python
import json, pandas as pd
registro = pd.read_json("../resultados/registro.jsonl", lines=True)
resumo   = json.load(open("../resultados/resumo.json", encoding="utf-8"))

for metodo, c in resumo["curvas"].items():
    plt.plot(c["avaliacoes"], c["melhor_media"], label=metodo)
plt.xlabel("avaliacoes de aptidao"); plt.ylabel("melhor aptidao")
plt.legend()
```
