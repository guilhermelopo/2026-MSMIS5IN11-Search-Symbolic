"""Testes dos operadores geneticos.

O que precisa ficar provado: os operadores nunca produzem lixo, e o
tamanho da molecula pode variar (sem isso o GA so rearranja o que ja
tinha). A validade nao e verificada pelos operadores -- e garantida pela
representacao -- entao e aqui que isso vira teste.
"""

import random

from src import operators as ops
from src import representation as rep


def cadeias_de_teste(n=40, semente=0):
    rng = random.Random(semente)
    return [rep.aleatoria(rng) for _ in range(n)]


# --------------------------------------------------------------------------
# Mutacao
# --------------------------------------------------------------------------

def test_mutacao_sempre_devolve_molecula_valida():
    rng = random.Random(0)
    for cadeia in cadeias_de_teste():
        for _ in range(10):
            filho = ops.mutar(cadeia, rng, p=0.3)
            assert rep.valida(filho), f"{cadeia} -> {filho}"


def test_mutacao_nunca_devolve_cadeia_vazia():
    rng = random.Random(0)
    # p=1.0 forca uma operacao em toda ficha, inclusive remocoes em serie.
    for cadeia in cadeias_de_teste(n=60):
        assert ops.mutar(cadeia, rng, p=1.0) != ""


def test_mutacao_com_p_zero_nao_muda_nada():
    rng = random.Random(0)
    for cadeia in cadeias_de_teste(n=10):
        assert ops.mutar(cadeia, rng, p=0.0) == cadeia


def test_mutacao_muda_alguma_coisa_com_p_alto():
    rng = random.Random(0)
    originais = cadeias_de_teste(n=30)
    mudou = sum(1 for c in originais if ops.mutar(c, rng, p=0.5) != c)
    assert mudou > 20


def test_mutacao_altera_o_tamanho_nos_dois_sentidos():
    """Inserir e remover precisam acontecer de verdade.

    Se so houvesse troca, o numero de fichas nunca mudaria e o GA jamais
    exploraria moleculas maiores ou menores que as iniciais.
    """
    rng = random.Random(7)
    base = rep.aleatoria(random.Random(3), n_min=12, n_max=12)
    tamanhos = {len(rep.fichas(ops.mutar(base, rng, p=0.4))) for _ in range(200)}
    assert any(t > 12 for t in tamanhos), "nenhuma insercao em 200 tentativas"
    assert any(t < 12 for t in tamanhos), "nenhuma remocao em 200 tentativas"


def test_mutacao_respeita_a_semente():
    base = rep.aleatoria(random.Random(1))
    a = ops.mutar(base, random.Random(99), p=0.3)
    b = ops.mutar(base, random.Random(99), p=0.3)
    assert a == b


# --------------------------------------------------------------------------
# Cruzamento
# --------------------------------------------------------------------------

def test_cruzamento_devolve_duas_moleculas_validas():
    rng = random.Random(0)
    cadeias = cadeias_de_teste()
    for a, b in zip(cadeias, cadeias[1:]):
        f1, f2 = ops.cruzar(a, b, rng)
        assert rep.valida(f1) and rep.valida(f2)


def test_cruzamento_nunca_devolve_cadeia_vazia():
    rng = random.Random(0)
    cadeias = cadeias_de_teste()
    for a, b in zip(cadeias, cadeias[1:]):
        f1, f2 = ops.cruzar(a, b, rng)
        assert f1 and f2


def test_cruzamento_funciona_com_cadeia_de_uma_ficha():
    rng = random.Random(0)
    uma = rep.alfabeto()[0]
    f1, f2 = ops.cruzar(uma, rep.aleatoria(random.Random(2)), rng)
    assert rep.valida(f1) and rep.valida(f2)


def test_cruzamento_mistura_os_dois_pais():
    rng = random.Random(5)
    a = rep.aleatoria(random.Random(10), n_min=15, n_max=15)
    b = rep.aleatoria(random.Random(20), n_min=15, n_max=15)
    diferentes = sum(1 for _ in range(50) if ops.cruzar(a, b, rng)[0] not in (a, b))
    assert diferentes > 40


# --------------------------------------------------------------------------
# Torneio
# --------------------------------------------------------------------------

def test_torneio_com_k_igual_a_populacao_escolhe_o_melhor():
    pop = ["a", "b", "c", "d"]
    apt = [0.1, 0.9, 0.3, 0.2]
    assert ops.torneio(pop, apt, random.Random(0), k=4) == "b"


def test_torneio_nunca_escolhe_fora_da_populacao():
    pop = ["a", "b", "c", "d", "e"]
    apt = [0.1, 0.9, 0.3, 0.2, 0.5]
    rng = random.Random(0)
    assert all(ops.torneio(pop, apt, rng, k=3) in pop for _ in range(200))


def test_torneio_nao_e_sempre_o_melhor():
    """A folga do torneio e o que impede o colapso na primeira geracao."""
    pop = ["a", "b", "c", "d", "e"]
    apt = [0.1, 0.9, 0.3, 0.2, 0.5]
    rng = random.Random(0)
    escolhidos = {ops.torneio(pop, apt, rng, k=2) for _ in range(200)}
    assert len(escolhidos) > 1


def test_torneio_com_k_maior_que_a_populacao():
    pop = ["a", "b"]
    apt = [0.1, 0.9]
    assert ops.torneio(pop, apt, random.Random(0), k=10) == "b"


def test_torneio_com_aptidoes_empatadas():
    pop = ["a", "b", "c"]
    apt = [0.5, 0.5, 0.5]
    assert ops.torneio(pop, apt, random.Random(0), k=3) in pop
