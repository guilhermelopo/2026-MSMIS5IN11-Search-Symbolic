"""Dominance, front de Pareto et hypervolume.

Ce module est pur : il travaille sur des tuples de nombres et n'importe ni
RDKit ni selfies. C'est le plus facile a tester avec des reponses connues
calculees a la main -- et c'est pour cela que ses tests sont les plus
precieux du projet.

Convention dans tout le module : TOUS les objectifs sont a MAXIMISER.
Le SA score, qui est a minimiser, arrive deja converti en `sa_norm` par le
module objectives. Melanger les directions est l'erreur classique ici.
"""

from __future__ import annotations

from typing import Sequence

Point = Sequence[float]


def domine(a: Point, b: Point) -> bool:
    """A domine B : au moins aussi bon partout, strictement meilleur quelque part."""
    au_moins_egal = all(x >= y for x, y in zip(a, b))
    strictement_meilleur = any(x > y for x, y in zip(a, b))
    return au_moins_egal and strictement_meilleur


def indices_non_domines(points: Sequence[Point]) -> list[int]:
    """Indices du front de Pareto.

    Complexite quadratique. Pour les centaines ou milliers de points de ce
    projet, c'est sans importance ; au-dela de ~1e5, passer a un algorithme
    par balayage.
    """
    return [
        i
        for i, p in enumerate(points)
        if not any(domine(q, p) for j, q in enumerate(points) if j != i)
    ]


def front(points: Sequence[Point]) -> list[Point]:
    return [points[i] for i in indices_non_domines(points)]


def couches(points: Sequence[Point]) -> list[list[int]]:
    """Tri par couches de dominance, comme dans NSGA-II.

    La couche 0 est le front de Pareto ; la couche 1 est le front de ce qui
    reste ; et ainsi de suite. Couche plus basse, meilleur.
    """
    restants = list(range(len(points)))
    sortie: list[list[int]] = []
    while restants:
        sous = [points[i] for i in restants]
        locaux = set(indices_non_domines(sous))
        couche = [restants[i] for i in sorted(locaux)]
        sortie.append(couche)
        restants = [r for k, r in enumerate(restants) if k not in locaux]
    return sortie


def distance_agglomeration(points: Sequence[Point]) -> list[float]:
    """Distance d'agglomeration de NSGA-II, pour departager dans une couche.

    Les points situes dans les zones vides du front recoivent une distance
    plus grande et sont preferes, de sorte que le front s'etale au lieu de
    s'entasser. Les extremes recoivent l'infini et sont toujours preserves.
    """
    n = len(points)
    if n == 0:
        return []
    if n <= 2:
        return [float("inf")] * n

    n_obj = len(points[0])
    dist = [0.0] * n
    for m in range(n_obj):
        ordre = sorted(range(n), key=lambda i: points[i][m])
        dist[ordre[0]] = float("inf")
        dist[ordre[-1]] = float("inf")
        mini = points[ordre[0]][m]
        maxi = points[ordre[-1]][m]
        plage = maxi - mini
        if plage == 0:
            continue
        for k in range(1, n - 1):
            if dist[ordre[k]] == float("inf"):
                continue
            precedent = points[ordre[k - 1]][m]
            suivant = points[ordre[k + 1]][m]
            dist[ordre[k]] += (suivant - precedent) / plage
    return dist


def hypervolume_2d(points: Sequence[Point], reference: Point = (0.0, 0.0)) -> float:
    """Aire dominee par le front relativement a un point de reference.

    Un seul nombre pour comparer deux fronts entiers : plus grand est
    meilleur. Valable pour deux objectifs, ce qui est le cas de ce projet.

    Le point de reference doit etre le MEME dans toutes les comparaisons du
    rapport, sinon les nombres ne veulent rien dire. (0, 0) est le choix
    naturel ici, les deux objectifs vivant deja dans 0..1.
    """
    utiles = [
        p for p in front(points)
        if p[0] > reference[0] and p[1] > reference[1]
    ]
    if not utiles:
        return 0.0

    # Balayage de droite a gauche, en sommant des bandes rectangulaires.
    utiles.sort(key=lambda p: p[0], reverse=True)
    aire = 0.0
    y_precedent = reference[1]
    for x, y in utiles:
        if y > y_precedent:
            aire += (x - reference[0]) * (y - y_precedent)
            y_precedent = y
    return aire
