"""Linhas de base: busca aleatoria e hill climbing.

Sem elas, uma curva de aptidao subindo nao prova nada -- qualquer metodo
que guarde o melhor resultado produz uma curva subindo.

As duas devolvem o mesmo `Resultado` do modulo ga e gastam o mesmo
`Orcamento`, de modo que os tres metodos sao comparaveis sem ajuste.
O eixo horizontal dos graficos deve ser `historico_avaliacoes`, nunca o
numero de iteracoes: iteracoes de metodos diferentes custam coisas
diferentes.
"""

from __future__ import annotations

import random

from . import objectives as obj
from . import representation as rep
from .ga import Achado, Arquivo, Resultado


def _registrar(
    arquivo: Arquivo,
    cadeia: str,
    n: obj.Notas,
    apt: float,
    passo: int,
    metodo: str,
    semente: int,
) -> None:
    if not n.valida:
        return
    smi = rep.canonico(cadeia)
    if smi is None:
        return
    arquivo.registrar(
        Achado(
            smiles=smi, selfies=cadeia, qed=n.qed, sa=n.sa, sa_norm=n.sa_norm,
            aptidao=apt, geracao=passo, metodo=metodo, semente=semente,
        )
    )


def busca_aleatoria(
    semente: int,
    orcamento: obj.Orcamento,
    peso_qed: float = 0.5,
    passo_registro: int = 100,
) -> Resultado:
    """O piso. Sorteia cadeias ao acaso e guarda a melhor.

    Se o algoritmo genetico empatar com isto, ele nao esta evoluindo nada.
    A causa habitual e um bug: aptidao calculada mas nao usada na selecao,
    ou taxa de mutacao tao alta que destroi qualquer heranca.
    """
    rng = random.Random(semente)
    arquivo = Arquivo()
    hist_melhor: list[float] = []
    hist_media: list[float] = []
    hist_aval: list[int] = []
    melhor = 0.0
    soma, conta = 0.0, 0
    passo = 0

    while not orcamento.esgotado:
        cadeia = rep.aleatoria(rng)
        n = obj.notas(cadeia, orcamento)
        apt = peso_qed * n.qed + (1 - peso_qed) * n.sa_norm
        _registrar(arquivo, cadeia, n, apt, passo, "aleatoria", semente)
        melhor = max(melhor, apt)
        soma += apt
        conta += 1
        passo += 1
        if passo % passo_registro == 0:
            hist_melhor.append(melhor)
            hist_media.append(soma / conta)
            hist_aval.append(orcamento.usadas)

    if not hist_melhor or hist_aval[-1] != orcamento.usadas:
        hist_melhor.append(melhor)
        hist_media.append(soma / max(1, conta))
        hist_aval.append(orcamento.usadas)

    return Resultado(
        metodo="aleatoria", semente=semente, arquivo=arquivo,
        historico_melhor=hist_melhor, historico_media=hist_media,
        historico_avaliacoes=hist_aval, populacao_final=[],
    )


def hill_climbing(
    semente: int,
    orcamento: obj.Orcamento,
    peso_qed: float = 0.5,
    p_mutacao: float = 0.15,
    reinicios: bool = False,
    paciencia: int = 500,
    passo_registro: int = 100,
) -> Resultado:
    """Busca local gananciosa: muta e so aceita melhorias.

    Sobe rapido no comeco e depois empaca. Ao chegar num pico local, toda
    mutacao piora e ele fica preso ali -- e esse o contraste com a
    populacao do GA, que cobre varios picos ao mesmo tempo.

    `reinicios=True` reinicia de uma cadeia nova depois de `paciencia`
    tentativas sem melhora. E uma linha de base mais forte e mais honesta;
    rode as duas e mostre a diferenca.
    """
    from . import operators as ops

    rng = random.Random(semente)
    arquivo = Arquivo()
    hist_melhor: list[float] = []
    hist_media: list[float] = []
    hist_aval: list[int] = []

    atual = rep.aleatoria(rng)
    n_atual = obj.notas(atual, orcamento)
    apt_atual = peso_qed * n_atual.qed + (1 - peso_qed) * n_atual.sa_norm
    _registrar(arquivo, atual, n_atual, apt_atual, 0, "hill", semente)

    melhor = apt_atual
    soma, conta = apt_atual, 1
    passo, sem_melhora = 0, 0

    while not orcamento.esgotado:
        candidato = ops.mutar(atual, rng, p_mutacao)
        n = obj.notas(candidato, orcamento)
        apt = peso_qed * n.qed + (1 - peso_qed) * n.sa_norm
        _registrar(arquivo, candidato, n, apt, passo, "hill", semente)

        if apt > apt_atual:
            atual, apt_atual = candidato, apt
            sem_melhora = 0
        else:
            sem_melhora += 1
            if reinicios and sem_melhora >= paciencia:
                atual = rep.aleatoria(rng)
                n_atual = obj.notas(atual, orcamento)
                apt_atual = peso_qed * n_atual.qed + (1 - peso_qed) * n_atual.sa_norm
                sem_melhora = 0

        melhor = max(melhor, apt)
        soma += apt
        conta += 1
        passo += 1
        if passo % passo_registro == 0:
            hist_melhor.append(melhor)
            hist_media.append(soma / conta)
            hist_aval.append(orcamento.usadas)

    if not hist_melhor or hist_aval[-1] != orcamento.usadas:
        hist_melhor.append(melhor)
        hist_media.append(soma / max(1, conta))
        hist_aval.append(orcamento.usadas)

    return Resultado(
        metodo="hill" + ("_reinicios" if reinicios else ""),
        semente=semente, arquivo=arquivo,
        historico_melhor=hist_melhor, historico_media=hist_media,
        historico_avaliacoes=hist_aval, populacao_final=[atual],
    )
