"""Tests des operateurs genetiques.

Ce qu'il faut prouver : les operateurs ne produisent jamais de dechet, et la
taille de la molecule peut varier (sans quoi le GA ne fait que rearranger ce
qu'il avait deja). La validite n'est pas verifiee par les operateurs -- elle
est garantie par la representation -- c'est donc ici qu'elle devient un test.
"""

import random

from src import operators as ops
from src import representation as rep


def chaines_de_test(n=40, graine=0):
    rng = random.Random(graine)
    return [rep.aleatoire(rng) for _ in range(n)]


# --------------------------------------------------------------------------
# Mutation
# --------------------------------------------------------------------------

def test_mutation_rend_toujours_une_molecule_valide():
    rng = random.Random(0)
    for chaine in chaines_de_test():
        for _ in range(10):
            enfant = ops.muter(chaine, rng, p=0.3)
            assert rep.decodage_ok(enfant), f"{chaine} -> {enfant}"


def test_mutation_ne_rend_jamais_une_chaine_vide():
    rng = random.Random(0)
    # p=1.0 force une operation sur chaque jeton, suppressions en serie incluses.
    for chaine in chaines_de_test(n=60):
        assert ops.muter(chaine, rng, p=1.0) != ""


def test_mutation_avec_p_nul_ne_change_rien():
    rng = random.Random(0)
    for chaine in chaines_de_test(n=10):
        assert ops.muter(chaine, rng, p=0.0) == chaine


def test_mutation_change_quelque_chose_avec_p_eleve():
    rng = random.Random(0)
    originales = chaines_de_test(n=30)
    change = sum(1 for c in originales if ops.muter(c, rng, p=0.5) != c)
    assert change > 20


def test_mutation_modifie_la_taille_dans_les_deux_sens():
    """Inserer et supprimer doivent vraiment se produire.

    S'il n'y avait que la substitution, le nombre de jetons ne changerait
    jamais et le GA n'explorerait jamais de molecules plus grandes ou plus
    petites que les initiales.
    """
    rng = random.Random(7)
    base = rep.aleatoire(random.Random(3), n_min=12, n_max=12)
    tailles = {len(rep.jetons(ops.muter(base, rng, p=0.4))) for _ in range(200)}
    assert any(t > 12 for t in tailles), "aucune insertion en 200 essais"
    assert any(t < 12 for t in tailles), "aucune suppression en 200 essais"


def test_mutation_respecte_la_graine():
    base = rep.aleatoire(random.Random(1))
    a = ops.muter(base, random.Random(99), p=0.3)
    b = ops.muter(base, random.Random(99), p=0.3)
    assert a == b


# --------------------------------------------------------------------------
# Croisement
# --------------------------------------------------------------------------

def test_croisement_rend_deux_molecules_valides():
    rng = random.Random(0)
    chaines = chaines_de_test()
    for a, b in zip(chaines, chaines[1:]):
        e1, e2 = ops.croiser(a, b, rng)
        assert rep.decodage_ok(e1) and rep.decodage_ok(e2)


def test_croisement_ne_rend_jamais_une_chaine_vide():
    rng = random.Random(0)
    chaines = chaines_de_test()
    for a, b in zip(chaines, chaines[1:]):
        e1, e2 = ops.croiser(a, b, rng)
        assert e1 and e2


def test_croisement_fonctionne_avec_une_chaine_d_un_jeton():
    rng = random.Random(0)
    un = rep.alphabet()[0]
    e1, e2 = ops.croiser(un, rep.aleatoire(random.Random(2)), rng)
    assert rep.decodage_ok(e1) and rep.decodage_ok(e2)


def test_croisement_melange_les_deux_parents():
    rng = random.Random(5)
    a = rep.aleatoire(random.Random(10), n_min=15, n_max=15)
    b = rep.aleatoire(random.Random(20), n_min=15, n_max=15)
    differents = sum(1 for _ in range(50) if ops.croiser(a, b, rng)[0] not in (a, b))
    assert differents > 40


# --------------------------------------------------------------------------
# Tournoi
# --------------------------------------------------------------------------

def test_tournoi_avec_k_egal_a_la_population_choisit_le_meilleur():
    pop = ["a", "b", "c", "d"]
    apt = [0.1, 0.9, 0.3, 0.2]
    assert ops.tournoi(pop, apt, random.Random(0), k=4) == "b"


def test_tournoi_ne_choisit_jamais_hors_de_la_population():
    pop = ["a", "b", "c", "d", "e"]
    apt = [0.1, 0.9, 0.3, 0.2, 0.5]
    rng = random.Random(0)
    assert all(ops.tournoi(pop, apt, rng, k=3) in pop for _ in range(200))


def test_tournoi_ne_choisit_pas_toujours_le_meilleur():
    """Le jeu du tournoi est ce qui empeche l'effondrement des la premiere generation."""
    pop = ["a", "b", "c", "d", "e"]
    apt = [0.1, 0.9, 0.3, 0.2, 0.5]
    rng = random.Random(0)
    choisis = {ops.tournoi(pop, apt, rng, k=2) for _ in range(200)}
    assert len(choisis) > 1


def test_tournoi_avec_k_superieur_a_la_population():
    pop = ["a", "b"]
    apt = [0.1, 0.9]
    assert ops.tournoi(pop, apt, random.Random(0), k=10) == "b"


def test_tournoi_avec_aptitudes_egales():
    pop = ["a", "b", "c"]
    apt = [0.5, 0.5, 0.5]
    assert ops.tournoi(pop, apt, random.Random(0), k=3) in pop
