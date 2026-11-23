"""Tests des fonctions objectif et du budget.

L'erreur de signe dans la normalisation du SA est le bug le plus probable de
ce projet : les deux notes ont des directions opposees, et les intervertir
produit un programme qui tourne, converge et optimise exactement le
contraire de ce qui est voulu, sans prevenir. D'ou les tests de direction
ci-dessous.
"""

import pytest

from src import objectives as obj
from src import representation as rep

ASPIRINE = "CC(=O)Oc1ccccc1C(=O)O"
CAFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"


# --------------------------------------------------------------------------
# Normalisation du SA
# --------------------------------------------------------------------------

def test_sa_norm_aux_extremes():
    assert obj.normaliser_sa(1.0) == pytest.approx(1.0)
    assert obj.normaliser_sa(10.0) == pytest.approx(0.0)


def test_sa_norm_inverse_la_direction():
    # SA plus bas = plus facile a faire = sa_norm PLUS HAUT.
    assert obj.normaliser_sa(2.0) > obj.normaliser_sa(8.0)


def test_sa_norm_reste_dans_l_intervalle():
    assert obj.normaliser_sa(-5.0) == 1.0
    assert obj.normaliser_sa(99.0) == 0.0


# --------------------------------------------------------------------------
# Notes
# --------------------------------------------------------------------------

def test_notes_d_une_molecule_connue():
    n = obj.notes(rep.depuis_smiles(ASPIRINE))
    assert n.valide
    assert 0.0 <= n.qed <= 1.0
    assert 1.0 <= n.sa <= 10.0
    assert 0.0 <= n.sa_norm <= 1.0


def test_aspirine_est_facile_a_synthetiser():
    """Test de bon sens a reponse connue.

    L'aspirine est l'une des molecules les plus simples de la pharmacopee.
    Si son SA ressort eleve, c'est que sascorer a ete charge de travers.
    """
    n = obj.notes(rep.depuis_smiles(ASPIRINE))
    assert n.sa < 4.0, f"SA de l'aspirine = {n.sa}, attendu bas"
    assert n.sa_norm > 0.6


def test_objectifs_sont_le_couple_a_maximiser():
    n = obj.notes(rep.depuis_smiles(CAFEINE))
    assert n.objectifs == (n.qed, n.sa_norm)


def test_notes_invalides_ne_cassent_rien():
    assert obj.NOTES_INVALIDES.valide is False
    assert obj.NOTES_INVALIDES.qed == 0.0
    assert obj.NOTES_INVALIDES.sa_norm == 0.0


def test_chaine_degeneree_ne_leve_pas():
    """Le bug qui a interrompu une execution dans Colab.

    `[Ring1]` se decode en SMILES vide. `MolFromSmiles("")` renvoie une
    molecule sans atome, pas None ; le sascorer renvoyait None pour elle et
    `float(None)` levait TypeError au milieu de la troisieme generation.
    """
    n = obj.notes("[Ring1]")
    assert n.valide is False
    assert n.qed == 0.0
    assert n.sa_norm == 0.0


def test_aptitude_d_une_chaine_degeneree_est_nulle():
    """L'aptitude doit rester un nombre : le tournoi compare avec `max`."""
    a = obj.aptitude("[Ring1]")
    assert isinstance(a, float)
    assert a == 0.0


@pytest.mark.parametrize("chaine", ["[Ring1]", "[Ring2]", "[Branch1]", "[#Branch2]"])
def test_jetons_isoles_ne_cassent_jamais(chaine):
    """N'importe quel jeton isole doit pouvoir etre note sans exception."""
    n = obj.notes(chaine)
    assert isinstance(n.qed, float)
    assert isinstance(n.sa, float)
    assert isinstance(n.sa_norm, float)


# --------------------------------------------------------------------------
# Aptitude
# --------------------------------------------------------------------------

def test_aptitude_reste_entre_zero_et_un():
    chaine = rep.depuis_smiles(ASPIRINE)
    for w in (0.0, 0.25, 0.5, 0.75, 1.0):
        assert 0.0 <= obj.aptitude(chaine, poids_qed=w) <= 1.0


def test_poids_un_n_utilise_que_le_qed():
    chaine = rep.depuis_smiles(ASPIRINE)
    n = obj.notes(chaine)
    assert obj.aptitude(chaine, poids_qed=1.0) == pytest.approx(n.qed)


def test_poids_zero_n_utilise_que_la_synthese():
    chaine = rep.depuis_smiles(ASPIRINE)
    n = obj.notes(chaine)
    assert obj.aptitude(chaine, poids_qed=0.0) == pytest.approx(n.sa_norm)


# --------------------------------------------------------------------------
# Budget
# --------------------------------------------------------------------------

def test_budget_compte_chaque_evaluation():
    b = obj.Budget(3)
    chaine = rep.depuis_smiles("CCO")
    for _ in range(3):
        obj.notes(chaine, b)
    assert b.utilisees == 3
    assert b.epuise
    assert b.restantes == 0


def test_budget_commence_vide():
    b = obj.Budget(10)
    assert not b.epuise
    assert b.restantes == 10


def test_budget_ne_devient_pas_negatif():
    b = obj.Budget(2)
    b.depenser(5)
    assert b.restantes == 0


def test_le_cache_ne_change_pas_le_compte_du_budget():
    """Le cache economise du temps, pas du budget.

    S'il economisait du budget, deux methodes ayant la meme limite feraient
    des nombres differents d'evaluations reelles et la comparaison serait
    injuste.
    """
    obj.vider_cache()
    b = obj.Budget(100)
    chaine = rep.depuis_smiles(ASPIRINE)
    for _ in range(10):
        obj.notes(chaine, b)
    assert b.utilisees == 10
