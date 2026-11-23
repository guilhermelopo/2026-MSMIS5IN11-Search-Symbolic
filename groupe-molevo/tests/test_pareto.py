"""Tests de dominance, front de Pareto et hypervolume.

Ce fichier n'importe ni RDKit ni selfies : pareto.py est du Python pur.
C'est donc le seul qui s'execute sur n'importe quelle machine, et aussi
celui dont les reponses sont les plus faciles a verifier a la main -- chaque
cas ci-dessous a ete calcule sur papier avant de devenir du code.
"""

from src.pareto import (
    couches,
    distance_agglomeration,
    domine,
    front,
    hypervolume_2d,
    indices_non_domines,
)


# --------------------------------------------------------------------------
# Dominance
# --------------------------------------------------------------------------

def test_domine_strictement_meilleur_partout():
    assert domine((0.8, 0.8), (0.5, 0.5))


def test_domine_meilleur_sur_un_et_egal_sur_l_autre():
    assert domine((0.8, 0.5), (0.5, 0.5))


def test_ne_domine_pas_quand_egal():
    assert not domine((0.5, 0.5), (0.5, 0.5))


def test_ne_domine_pas_quand_c_est_un_compromis():
    # Le cas qui donne son sens au projet : meilleur en QED, pire en synthese.
    assert not domine((0.8, 0.2), (0.2, 0.8))
    assert not domine((0.2, 0.8), (0.8, 0.2))


def test_domine_n_est_pas_symetrique():
    a, b = (0.9, 0.9), (0.1, 0.1)
    assert domine(a, b)
    assert not domine(b, a)


# --------------------------------------------------------------------------
# Front
# --------------------------------------------------------------------------

def test_front_d_un_seul_point():
    assert indices_non_domines([(0.5, 0.5)]) == [0]


def test_front_avec_un_domine_au_milieu():
    points = [
        (0.9, 0.1),   # 0  front
        (0.5, 0.5),   # 1  front
        (0.1, 0.9),   # 2  front
        (0.4, 0.4),   # 3  domine par 1
        (0.2, 0.2),   # 4  domine par tous
    ]
    assert indices_non_domines(points) == [0, 1, 2]


def test_front_quand_tous_sont_egaux():
    # Des points identiques ne se dominent pas : tous restent sur le front.
    points = [(0.5, 0.5)] * 4
    assert indices_non_domines(points) == [0, 1, 2, 3]


def test_front_quand_un_domine_tout_le_monde():
    points = [(0.1, 0.1), (0.2, 0.2), (1.0, 1.0)]
    assert indices_non_domines(points) == [2]


def test_front_renvoie_les_bons_points():
    # (0.4, 0.05) est domine par (0.9, 0.1) : pire sur les deux objectifs.
    # Attention : (0.4, 0.4) ne serait domine par aucun des extremes -- un
    # point median ne quitte le front que s'il perd sur TOUS les objectifs.
    points = [(0.9, 0.1), (0.4, 0.05), (0.1, 0.9)]
    assert front(points) == [(0.9, 0.1), (0.1, 0.9)]


# --------------------------------------------------------------------------
# Couches (NSGA-II)
# --------------------------------------------------------------------------

def test_couches_separe_en_deux():
    points = [
        (0.9, 0.1),   # 0  couche 0
        (0.1, 0.9),   # 1  couche 0
        (0.5, 0.05),  # 2  couche 1 (domine par 0)
        (0.05, 0.5),  # 3  couche 1 (domine par 1)
    ]
    assert couches(points) == [[0, 1], [2, 3]]


def test_couches_couvre_chaque_point_une_seule_fois():
    points = [(0.9, 0.1), (0.5, 0.5), (0.1, 0.9), (0.4, 0.4), (0.2, 0.2)]
    tous = [i for couche in couches(points) for i in couche]
    assert sorted(tous) == list(range(len(points)))


# --------------------------------------------------------------------------
# Distance d'agglomeration
# --------------------------------------------------------------------------

def test_agglomeration_les_extremes_sont_infinis():
    points = [(0.1, 0.9), (0.5, 0.5), (0.9, 0.1)]
    d = distance_agglomeration(points)
    assert d[0] == float("inf")
    assert d[2] == float("inf")
    assert 0 < d[1] < float("inf")


def test_agglomeration_prefere_le_point_le_plus_isole():
    # 1 occupe un grand vide ; 2 est colle a 3. 1 doit marquer davantage.
    points = [(0.0, 1.0), (0.5, 0.5), (0.9, 0.12), (0.92, 0.1)]
    d = distance_agglomeration(points)
    assert d[1] > d[2]


# --------------------------------------------------------------------------
# Hypervolume
# --------------------------------------------------------------------------

def test_hypervolume_d_un_point_est_le_rectangle():
    assert hypervolume_2d([(0.5, 0.4)]) == 0.5 * 0.4


def test_hypervolume_escalier_de_trois_points():
    # Calcul a la main :
    #   y de 0,0 a 0,3 -> x jusqu'a 0,8  ->  0,8 * 0,3 = 0,24
    #   y de 0,3 a 0,6 -> x jusqu'a 0,5  ->  0,5 * 0,3 = 0,15
    #   y de 0,6 a 0,9 -> x jusqu'a 0,2  ->  0,2 * 0,3 = 0,06
    #                                        total     = 0,45
    points = [(0.8, 0.3), (0.5, 0.6), (0.2, 0.9)]
    assert abs(hypervolume_2d(points) - 0.45) < 1e-9


def test_hypervolume_ignore_les_domines():
    front_pur = [(0.8, 0.3), (0.5, 0.6), (0.2, 0.9)]
    avec_bruit = front_pur + [(0.3, 0.2), (0.1, 0.1), (0.4, 0.5)]
    assert abs(hypervolume_2d(avec_bruit) - hypervolume_2d(front_pur)) < 1e-9


def test_hypervolume_croit_quand_le_front_s_ameliore():
    pire = [(0.5, 0.5)]
    meilleur = [(0.6, 0.5), (0.5, 0.6)]
    assert hypervolume_2d(meilleur) > hypervolume_2d(pire)


def test_hypervolume_est_nul_au_point_de_reference():
    assert hypervolume_2d([(0.0, 0.0)]) == 0.0


def test_hypervolume_d_un_ensemble_vide():
    assert hypervolume_2d([]) == 0.0
