# Dados

## Procedencia

Este projeto **nao depende de nenhum conjunto de dados baixado**. As moleculas
sao geradas pelo proprio algoritmo, a partir do alfabeto robusto da biblioteca
`selfies`.

As duas funcoes objetivo trazem seus proprios dados embutidos:

| Fonte | O que e | Licenca |
|---|---|---|
| `rdkit.Chem.QED` | parametros das funcoes de desejabilidade, ajustados por Bickerton et al. (2012) a farmacos orais aprovados | BSD-3-Clause, junto com o RDKit |
| `rdkit Contrib/SA_Score` | tabela de frequencia de fragmentos derivada do PubChem, de Ertl e Schuffenhauer (2009) | BSD-3-Clause, junto com o RDKit |

Nenhum dado pessoal e nenhum dado de saude real e usado em qualquer etapa.

## Se forem usar sementes de moleculas conhecidas

A funcao `representation.populacao_inicial` aceita uma lista de SMILES para
comecar a busca a partir de moleculas reais em vez de cadeias aleatorias.

Se usarem um conjunto publico para isso (ChEMBL, ZINC, o conjunto de referencia
do GuacaMol):

1. Escreva aqui a fonte exata, a versao e a licenca.
2. Adicione um script de download em `src/`.
3. **Nao comite os arquivos baixados** -- o `.gitignore` ja cobre `data/`.
4. Verifique se a licenca permite o uso; varias proibem redistribuicao.

## Comparar os dois pontos de partida e um bom experimento

Partir de moleculas reais costuma convergir mais rapido, mas enviesa o
resultado para a vizinhanca delas. Rodar as duas configuracoes e comparar a
diversidade final custa pouco e rende uma secao inteira do relatorio.
