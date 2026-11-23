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
    """A promessa do SELFIES: toda sequencia de fichas DECODIFICA.

    E esta a propriedade exata que a biblioteca garante, e e so ela que o
    algoritmo genetico precisa: nenhuma mutacao produz uma cadeia que o
    decodificador recuse.
    """
    rng = random.Random(0)
    falhas = []
    for _ in range(2000):
        cadeia = rep.aleatoria(rng, n_min=1, n_max=30)
        if not rep.decodifica_ok(cadeia):
            falhas.append(cadeia)
    assert not falhas, f"{len(falhas)} cadeias nao decodificam, ex.: {falhas[:3]}"


def test_cadeias_longas_dao_moleculas_com_atomos():
    """Decodificar sempre funciona, mas pode dar molecula VAZIA.

    Uma cadeia de uma ficha so pode ser degenerada (`[Ring1]` nao tem
    atomo). Com algumas fichas isso praticamente nao acontece -- e e por
    isso que a populacao inicial nunca usa n_min=1.
    """
    rng = random.Random(0)
    sem_atomo = [c for _ in range(500)
                 if not rep.valida(c := rep.aleatoria(rng, n_min=6, n_max=30))]
    assert len(sem_atomo) < 25, (
        f"{len(sem_atomo)}/500 cadeias sem atomo, ex.: {sem_atomo[:3]}")


def test_molecula_vazia_e_rejeitada():
    """A armadilha que causou um TypeError em producao.

    `Chem.MolFromSmiles("")` devolve uma molecula VAZIA, nao None. Se ela
    passar, o QED e o SA score recebem zero atomos e devolvem None, e
    `float(None)` levanta TypeError. `para_mol` tem de devolver None.
    """
    from rdkit import Chem

    vazia = Chem.MolFromSmiles("")
    assert vazia is not None, "o RDKit mudou: MolFromSmiles('') agora da None"
    assert vazia.GetNumAtoms() == 0

    # `[Ring1]` sozinha: uma ficha de fechamento de anel sem nenhum atomo.
    assert rep.decodifica_ok("[Ring1]")      # decodifica, como prometido
    assert not rep.valida("[Ring1]")         # mas nao tem atomo
    assert rep.para_mol("[Ring1]") is None
    assert rep.canonico("[Ring1]") is None   # nunca a cadeia vazia


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
