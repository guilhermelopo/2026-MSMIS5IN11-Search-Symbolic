"""Testes da representacao SELFIES.

O teste `test_robustez_do_alfabeto` e o mais importante do projeto: ele
verifica a propriedade em que todo o resto se apoia. Se ele falhar, o
algoritmo genetico inteiro perde a justificativa.
"""

import random

import pytest
from rdkit import Chem

from src import representation as rep

# Alguns farmacos conhecidos, para a ida e volta.
CONHECIDOS = [
    "CC(=O)Oc1ccccc1C(=O)O",            # aspirina
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",       # ibuprofeno
    "CC(=O)Nc1ccc(O)cc1",               # paracetamol
    "Cn1cnc2c1c(=O)n(C)c(=O)n2C",       # cafeina
    "CCO",                              # etanol
    "c1ccccc1",                         # benzeno
]


def canonizar(smiles: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


# --------------------------------------------------------------------------
# Ida e volta
# --------------------------------------------------------------------------

@pytest.mark.parametrize("smiles", CONHECIDOS)
def test_ida_e_volta_preserva_a_molecula(smiles):
    cadeia = rep.de_smiles(smiles)
    volta = rep.para_smiles(cadeia)
    assert canonizar(volta) == canonizar(smiles)


def test_canonico_concorda_com_o_rdkit():
    cadeia = rep.de_smiles("CC(=O)Oc1ccccc1C(=O)O")
    assert rep.canonico(cadeia) == canonizar("CC(=O)Oc1ccccc1C(=O)O")


def test_cadeias_diferentes_podem_dar_a_mesma_molecula():
    # Por isso o cache e o arquivo de achados indexam pelo SMILES canonico
    # e nao pela cadeia SELFIES.
    a = rep.de_smiles("CCO")
    b = rep.de_smiles("OCC")
    assert rep.canonico(a) == rep.canonico(b)


# --------------------------------------------------------------------------
# A propriedade central
# --------------------------------------------------------------------------

def test_robustez_do_alfabeto():
    """Toda sequencia de fichas do alfabeto robusto e uma molecula valida."""
    rng = random.Random(0)
    invalidas = []
    for _ in range(2000):
        cadeia = rep.aleatoria(rng, n_min=1, n_max=30)
        if not rep.valida(cadeia):
            invalidas.append(cadeia)
    assert not invalidas, f"{len(invalidas)} cadeias invalidas, ex.: {invalidas[:3]}"


def test_alfabeto_nao_e_vazio_e_e_estavel():
    a, b = rep.alfabeto(), rep.alfabeto()
    assert len(a) > 10
    assert a == b          # a ordem fixa e o que garante reprodutibilidade


def test_fichas_recompoem_a_cadeia():
    rng = random.Random(1)
    cadeia = rep.aleatoria(rng)
    assert "".join(rep.fichas(cadeia)) == cadeia


# --------------------------------------------------------------------------
# Geracao
# --------------------------------------------------------------------------

def test_aleatoria_respeita_a_semente():
    assert rep.aleatoria(random.Random(42)) == rep.aleatoria(random.Random(42))


def test_populacao_inicial_tem_o_tamanho_pedido():
    pop = rep.populacao_inicial(random.Random(0), 25)
    assert len(pop) == 25
    assert all(rep.valida(s) for s in pop)


def test_populacao_inicial_usa_as_sementes():
    pop = rep.populacao_inicial(random.Random(0), 10, sementes=CONHECIDOS[:3])
    canonicos = {rep.canonico(s) for s in pop}
    assert canonizar(CONHECIDOS[0]) in canonicos
    assert len(pop) == 10


def test_populacao_inicial_ignora_semente_invalida():
    pop = rep.populacao_inicial(random.Random(0), 5, sementes=["isto nao e um smiles"])
    assert len(pop) == 5
