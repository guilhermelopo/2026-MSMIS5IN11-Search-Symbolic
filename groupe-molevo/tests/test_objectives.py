"""Testes das funcoes objetivo e do orcamento.

O erro de sinal na normalizacao do SA e o bug mais provavel deste
projeto: as duas notas tem direcoes opostas, e trocar uma pela outra
produz um programa que roda, converge e otimiza exatamente o contrario
do pretendido, sem avisar. Daí os testes de direcao abaixo.
"""

import pytest

from src import objectives as obj
from src import representation as rep

ASPIRINA = "CC(=O)Oc1ccccc1C(=O)O"
CAFEINA = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"


# --------------------------------------------------------------------------
# Normalizacao do SA
# --------------------------------------------------------------------------

def test_sa_norm_nos_extremos():
    assert obj.normalizar_sa(1.0) == pytest.approx(1.0)
    assert obj.normalizar_sa(10.0) == pytest.approx(0.0)


def test_sa_norm_inverte_a_direcao():
    # SA menor = mais facil de fazer = sa_norm MAIOR.
    assert obj.normalizar_sa(2.0) > obj.normalizar_sa(8.0)


def test_sa_norm_fica_preso_no_intervalo():
    assert obj.normalizar_sa(-5.0) == 1.0
    assert obj.normalizar_sa(99.0) == 0.0


# --------------------------------------------------------------------------
# Notas
# --------------------------------------------------------------------------

def test_notas_de_uma_molecula_conhecida():
    n = obj.notas(rep.de_smiles(ASPIRINA))
    assert n.valida
    assert 0.0 <= n.qed <= 1.0
    assert 1.0 <= n.sa <= 10.0
    assert 0.0 <= n.sa_norm <= 1.0


def test_aspirina_e_facil_de_sintetizar():
    """Teste de sanidade com resposta conhecida.

    A aspirina e uma das moleculas mais simples da farmacopeia. Se o SA
    dela sair alto, o sascorer foi carregado errado.
    """
    n = obj.notas(rep.de_smiles(ASPIRINA))
    assert n.sa < 4.0, f"SA da aspirina = {n.sa}, esperado baixo"
    assert n.sa_norm > 0.6


def test_objetivos_sao_o_par_para_maximizar():
    n = obj.notas(rep.de_smiles(CAFEINA))
    assert n.objetivos == (n.qed, n.sa_norm)


def test_notas_invalidas_nao_quebram_nada():
    assert obj.NOTAS_INVALIDAS.valida is False
    assert obj.NOTAS_INVALIDAS.qed == 0.0
    assert obj.NOTAS_INVALIDAS.sa_norm == 0.0


def test_cadeia_degenerada_nao_levanta():
    """O bug que parou uma execucao no Colab.

    `[Ring1]` decodifica para um SMILES vazio. `MolFromSmiles("")` devolve
    uma molecula sem atomos, nao None; o sascorer devolvia None para ela e
    `float(None)` levantava TypeError no meio da terceira geracao.
    """
    n = obj.notas("[Ring1]")
    assert n.valida is False
    assert n.qed == 0.0
    assert n.sa_norm == 0.0


def test_aptidao_de_cadeia_degenerada_e_zero():
    """A aptidao tem de ser um numero: o torneio compara com `max`."""
    a = obj.aptidao("[Ring1]")
    assert isinstance(a, float)
    assert a == 0.0


@pytest.mark.parametrize("cadeia", ["[Ring1]", "[Ring2]", "[Branch1]", "[#Branch2]"])
def test_fichas_sozinhas_nunca_quebram(cadeia):
    """Qualquer ficha isolada tem de ser pontuavel sem excecao."""
    n = obj.notas(cadeia)
    assert isinstance(n.qed, float)
    assert isinstance(n.sa, float)
    assert isinstance(n.sa_norm, float)


# --------------------------------------------------------------------------
# Aptidao
# --------------------------------------------------------------------------

def test_aptidao_fica_entre_zero_e_um():
    cadeia = rep.de_smiles(ASPIRINA)
    for w in (0.0, 0.25, 0.5, 0.75, 1.0):
        assert 0.0 <= obj.aptidao(cadeia, peso_qed=w) <= 1.0


def test_peso_um_usa_so_o_qed():
    cadeia = rep.de_smiles(ASPIRINA)
    n = obj.notas(cadeia)
    assert obj.aptidao(cadeia, peso_qed=1.0) == pytest.approx(n.qed)


def test_peso_zero_usa_so_a_sintese():
    cadeia = rep.de_smiles(ASPIRINA)
    n = obj.notas(cadeia)
    assert obj.aptidao(cadeia, peso_qed=0.0) == pytest.approx(n.sa_norm)


# --------------------------------------------------------------------------
# Orcamento
# --------------------------------------------------------------------------

def test_orcamento_conta_cada_avaliacao():
    o = obj.Orcamento(3)
    cadeia = rep.de_smiles("CCO")
    for _ in range(3):
        obj.notas(cadeia, o)
    assert o.usadas == 3
    assert o.esgotado
    assert o.restantes == 0


def test_orcamento_comeca_vazio():
    o = obj.Orcamento(10)
    assert not o.esgotado
    assert o.restantes == 10


def test_orcamento_nao_fica_negativo():
    o = obj.Orcamento(2)
    o.gastar(5)
    assert o.restantes == 0


def test_cache_nao_altera_a_contagem_do_orcamento():
    """O cache economiza tempo, nao orcamento.

    Se ele economizasse orcamento, dois metodos com o mesmo limite fariam
    numeros diferentes de avaliacoes reais e a comparacao seria injusta.
    """
    obj.limpar_cache()
    o = obj.Orcamento(100)
    cadeia = rep.de_smiles(ASPIRINA)
    for _ in range(10):
        obj.notas(cadeia, o)
    assert o.usadas == 10
