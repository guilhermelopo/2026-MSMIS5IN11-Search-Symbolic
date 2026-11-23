"""Operadores geneticos sobre cadeias SELFIES.

Todos recebem um `random.Random` explicito em vez de usar o modulo
`random` global. Sem isso nao ha reprodutibilidade: duas execucoes com a
mesma semente dariam resultados diferentes se qualquer outra parte do
programa sorteasse um numero no meio.

A validade dos filhos nao e verificada em lugar nenhum -- e garantida
pela representacao. E esse o ponto do projeto.
"""

from __future__ import annotations

import random

from .representation import alfabeto, fichas


def mutar(cadeia: str, rng: random.Random, p: float = 0.15) -> str:
    """Aplica trocas, insercoes e remocoes ficha a ficha.

    Cada ficha sofre, com probabilidade p, uma das tres operacoes
    (igualmente provaveis entre si):

      trocar  -> muda um atomo ou uma ligacao
      inserir -> a molecula cresce
      remover -> a molecula encolhe

    Inserir e remover sao indispensaveis. Sem elas, o tamanho da molecula
    fica travado no tamanho inicial, e o tamanho e um dos descritores que
    mais pesam no QED.
    """
    alfa = alfabeto()
    saida: list[str] = []
    for ficha in fichas(cadeia):
        r = rng.random()
        if r < p / 3:
            saida.append(rng.choice(alfa))
        elif r < 2 * p / 3:
            saida.append(rng.choice(alfa))
            saida.append(ficha)
        elif r < p:
            continue
        else:
            saida.append(ficha)
    # Remover tudo e possivel quando p e alto; devolver "" quebraria o
    # decodificador, entao caimos numa ficha unica.
    return "".join(saida) or rng.choice(alfa)


def cruzar(a: str, b: str, rng: random.Random) -> tuple[str, str]:
    """Cruzamento de um ponto: corta as duas listas e troca as caudas."""
    ta, tb = fichas(a), fichas(b)
    i = rng.randrange(1, max(2, len(ta)))
    j = rng.randrange(1, max(2, len(tb)))
    f1 = "".join(ta[:i] + tb[j:])
    f2 = "".join(tb[:j] + ta[i:])
    return (f1 or a), (f2 or b)


def torneio(
    pop: list[str],
    aptidoes: list[float],
    rng: random.Random,
    k: int = 3,
) -> str:
    """Sorteia k individuos e devolve o melhor deles.

    Preferido a roleta por tres motivos praticos: nao se importa com a
    escala absoluta das notas, a pressao seletiva se ajusta so no k, e nao
    ha divisao por zero quando todas as aptidoes empatam.
    """
    k = min(k, len(pop))
    indices = rng.sample(range(len(pop)), k)
    melhor = max(indices, key=lambda i: aptidoes[i])
    return pop[melhor]
