"""Dominancia, frente de Pareto e hipervolume.

Este modulo e puro: trabalha com tuplas de numeros e nao importa RDKit
nem selfies. Isso o torna o mais facil de testar com respostas conhecidas
calculadas a mao -- e e por isso que os testes dele sao os mais valiosos
do projeto.

Convencao em todo o modulo: TODOS os objetivos sao para MAXIMIZAR.
O SA score, que e para minimizar, ja entra convertido em `sa_norm` pelo
modulo objectives. Misturar as direcoes e o erro classico aqui.
"""

from __future__ import annotations

from typing import Sequence

Ponto = Sequence[float]


def domina(a: Ponto, b: Ponto) -> bool:
    """A domina B: melhor ou igual em tudo, estritamente melhor em algo."""
    melhor_ou_igual = all(x >= y for x, y in zip(a, b))
    estritamente_melhor = any(x > y for x, y in zip(a, b))
    return melhor_ou_igual and estritamente_melhor


def indices_nao_dominados(pontos: Sequence[Ponto]) -> list[int]:
    """Indices da frente de Pareto.

    Complexidade quadratica. Para as centenas ou milhares de pontos deste
    projeto isso e irrelevante; se um dia passar de ~1e5, troque por um
    algoritmo de varredura.
    """
    return [
        i
        for i, p in enumerate(pontos)
        if not any(domina(q, p) for j, q in enumerate(pontos) if j != i)
    ]


def frente(pontos: Sequence[Ponto]) -> list[Ponto]:
    return [pontos[i] for i in indices_nao_dominados(pontos)]


def camadas(pontos: Sequence[Ponto]) -> list[list[int]]:
    """Ordenacao por camadas de dominancia, como no NSGA-II.

    Camada 0 e a frente de Pareto; camada 1 e a frente do que sobra; e
    assim por diante. Camada mais baixa, melhor.
    """
    restantes = list(range(len(pontos)))
    saida: list[list[int]] = []
    while restantes:
        sub = [pontos[i] for i in restantes]
        locais = set(indices_nao_dominados(sub))
        camada = [restantes[i] for i in sorted(locais)]
        saida.append(camada)
        restantes = [r for k, r in enumerate(restantes) if k not in locais]
    return saida


def distancia_aglomeracao(pontos: Sequence[Ponto]) -> list[float]:
    """Distancia de aglomeracao do NSGA-II, para desempatar dentro da camada.

    Pontos em regioes vazias da frente recebem distancia maior e sao
    preferidos, de modo que a frente se espalhe em vez de se amontoar.
    Os extremos recebem infinito e sao sempre preservados.
    """
    n = len(pontos)
    if n == 0:
        return []
    if n <= 2:
        return [float("inf")] * n

    n_obj = len(pontos[0])
    dist = [0.0] * n
    for m in range(n_obj):
        ordem = sorted(range(n), key=lambda i: pontos[i][m])
        dist[ordem[0]] = float("inf")
        dist[ordem[-1]] = float("inf")
        menor = pontos[ordem[0]][m]
        maior = pontos[ordem[-1]][m]
        faixa = maior - menor
        if faixa == 0:
            continue
        for k in range(1, n - 1):
            if dist[ordem[k]] == float("inf"):
                continue
            anterior = pontos[ordem[k - 1]][m]
            seguinte = pontos[ordem[k + 1]][m]
            dist[ordem[k]] += (seguinte - anterior) / faixa
    return dist


def hipervolume_2d(pontos: Sequence[Ponto], referencia: Ponto = (0.0, 0.0)) -> float:
    """Area dominada pela frente em relacao a um ponto de referencia.

    Um numero so para comparar duas frentes inteiras: maior e melhor.
    Vale apenas para dois objetivos, que e o caso deste projeto.

    O ponto de referencia tem de ser o MESMO em todas as comparacoes do
    relatorio, senao os numeros nao querem dizer nada. (0, 0) e a escolha
    natural aqui, porque os dois objetivos ja vivem em 0..1.
    """
    uteis = [
        p for p in frente(pontos)
        if p[0] > referencia[0] and p[1] > referencia[1]
    ]
    if not uteis:
        return 0.0

    # Varredura da direita para a esquerda, somando faixas retangulares.
    uteis.sort(key=lambda p: p[0], reverse=True)
    area = 0.0
    y_anterior = referencia[1]
    for x, y in uteis:
        if y > y_anterior:
            area += (x - referencia[0]) * (y - y_anterior)
            y_anterior = y
    return area
