"""Testes de dominancia, frente de Pareto e hipervolume.

Este arquivo nao importa RDKit nem selfies: pareto.py e Python puro.
Por isso e o unico que roda em qualquer maquina, e tambem o que tem as
respostas mais faceis de conferir a mao -- cada caso abaixo foi calculado
no papel antes de virar codigo.
"""

from src.pareto import (
    camadas,
    distancia_aglomeracao,
    domina,
    frente,
    hipervolume_2d,
    indices_nao_dominados,
)


# --------------------------------------------------------------------------
# Dominancia
# --------------------------------------------------------------------------

def test_domina_estritamente_melhor_nos_dois():
    assert domina((0.8, 0.8), (0.5, 0.5))


def test_domina_melhor_em_um_e_igual_no_outro():
    assert domina((0.8, 0.5), (0.5, 0.5))


def test_nao_domina_quando_e_igual():
    assert not domina((0.5, 0.5), (0.5, 0.5))


def test_nao_domina_quando_troca_um_pelo_outro():
    # O caso que da sentido ao projeto: melhor no QED, pior na sintese.
    assert not domina((0.8, 0.2), (0.2, 0.8))
    assert not domina((0.2, 0.8), (0.8, 0.2))


def test_domina_nao_e_simetrico():
    a, b = (0.9, 0.9), (0.1, 0.1)
    assert domina(a, b)
    assert not domina(b, a)


# --------------------------------------------------------------------------
# Frente
# --------------------------------------------------------------------------

def test_frente_de_um_ponto_so():
    assert indices_nao_dominados([(0.5, 0.5)]) == [0]


def test_frente_com_um_dominado_no_meio():
    pontos = [
        (0.9, 0.1),   # 0  frente
        (0.5, 0.5),   # 1  frente
        (0.1, 0.9),   # 2  frente
        (0.4, 0.4),   # 3  dominado por 1
        (0.2, 0.2),   # 4  dominado por todos
    ]
    assert indices_nao_dominados(pontos) == [0, 1, 2]


def test_frente_quando_todos_sao_iguais():
    # Pontos identicos nao se dominam: todos ficam na frente.
    pontos = [(0.5, 0.5)] * 4
    assert indices_nao_dominados(pontos) == [0, 1, 2, 3]


def test_frente_quando_um_domina_todo_mundo():
    pontos = [(0.1, 0.1), (0.2, 0.2), (1.0, 1.0)]
    assert indices_nao_dominados(pontos) == [2]


def test_frente_devolve_os_pontos_certos():
    # (0.4, 0.05) e dominado por (0.9, 0.1): pior nos dois objetivos.
    # Cuidado: (0.4, 0.4) NAO seria dominado por nenhum dos extremos --
    # um ponto mediano so sai da frente se perder em TODOS os objetivos.
    pontos = [(0.9, 0.1), (0.4, 0.05), (0.1, 0.9)]
    assert frente(pontos) == [(0.9, 0.1), (0.1, 0.9)]


# --------------------------------------------------------------------------
# Camadas (NSGA-II)
# --------------------------------------------------------------------------

def test_camadas_separa_em_duas():
    pontos = [
        (0.9, 0.1),   # 0  camada 0
        (0.1, 0.9),   # 1  camada 0
        (0.5, 0.05),  # 2  camada 1 (dominado por 0)
        (0.05, 0.5),  # 3  camada 1 (dominado por 1)
    ]
    assert camadas(pontos) == [[0, 1], [2, 3]]


def test_camadas_cobre_todos_os_pontos_uma_vez_so():
    pontos = [(0.9, 0.1), (0.5, 0.5), (0.1, 0.9), (0.4, 0.4), (0.2, 0.2)]
    todos = [i for camada in camadas(pontos) for i in camada]
    assert sorted(todos) == list(range(len(pontos)))


# --------------------------------------------------------------------------
# Distancia de aglomeracao
# --------------------------------------------------------------------------

def test_aglomeracao_extremos_sao_infinitos():
    pontos = [(0.1, 0.9), (0.5, 0.5), (0.9, 0.1)]
    d = distancia_aglomeracao(pontos)
    assert d[0] == float("inf")
    assert d[2] == float("inf")
    assert 0 < d[1] < float("inf")


def test_aglomeracao_prefere_o_ponto_mais_isolado():
    # 1 fica num vao grande; 2 fica colado em 3. 1 deve pontuar mais.
    pontos = [(0.0, 1.0), (0.5, 0.5), (0.9, 0.12), (0.92, 0.1)]
    d = distancia_aglomeracao(pontos)
    assert d[1] > d[2]


# --------------------------------------------------------------------------
# Hipervolume
# --------------------------------------------------------------------------

def test_hipervolume_de_um_ponto_e_o_retangulo():
    assert hipervolume_2d([(0.5, 0.4)]) == 0.5 * 0.4


def test_hipervolume_escada_de_tres_pontos():
    # Calculo a mao:
    #   y de 0,0 a 0,3 -> x ate 0,8  ->  0,8 * 0,3 = 0,24
    #   y de 0,3 a 0,6 -> x ate 0,5  ->  0,5 * 0,3 = 0,15
    #   y de 0,6 a 0,9 -> x ate 0,2  ->  0,2 * 0,3 = 0,06
    #                                    total     = 0,45
    pontos = [(0.8, 0.3), (0.5, 0.6), (0.2, 0.9)]
    assert abs(hipervolume_2d(pontos) - 0.45) < 1e-9


def test_hipervolume_ignora_os_dominados():
    frente_pura = [(0.8, 0.3), (0.5, 0.6), (0.2, 0.9)]
    com_lixo = frente_pura + [(0.3, 0.2), (0.1, 0.1), (0.4, 0.5)]
    assert abs(hipervolume_2d(com_lixo) - hipervolume_2d(frente_pura)) < 1e-9


def test_hipervolume_cresce_quando_a_frente_melhora():
    pior = [(0.5, 0.5)]
    melhor = [(0.6, 0.5), (0.5, 0.6)]
    assert hipervolume_2d(melhor) > hipervolume_2d(pior)


def test_hipervolume_e_zero_no_ponto_de_referencia():
    assert hipervolume_2d([(0.0, 0.0)]) == 0.0


def test_hipervolume_de_conjunto_vazio():
    assert hipervolume_2d([]) == 0.0
