"""Operateurs genetiques sur chaines SELFIES.

Tous recoivent un `random.Random` explicite au lieu d'utiliser le module
`random` global. Sans cela, pas de reproductibilite : deux executions avec
la meme graine donneraient des resultats differents des qu'une autre
partie du programme tirerait un nombre au milieu.

La validite des enfants n'est verifiee nulle part -- elle est garantie par
la representation. C'est tout l'interet du projet.
"""

from __future__ import annotations

import random

from .representation import alphabet, jetons


def muter(chaine: str, rng: random.Random, p: float = 0.15) -> str:
    """Applique substitutions, insertions et suppressions, jeton par jeton.

    Chaque jeton subit, avec la probabilite p, une des trois operations
    (equiprobables entre elles) :

      substituer -> change un atome ou une liaison
      inserer    -> la molecule grandit
      supprimer  -> la molecule retrecit

    Inserer et supprimer sont indispensables. Sans elles, la taille de la
    molecule reste figee a la taille initiale, or la taille est un des
    descripteurs qui pesent le plus dans le QED.
    """
    alpha = alphabet()
    sortie: list[str] = []
    for jeton in jetons(chaine):
        r = rng.random()
        if r < p / 3:
            sortie.append(rng.choice(alpha))
        elif r < 2 * p / 3:
            sortie.append(rng.choice(alpha))
            sortie.append(jeton)
        elif r < p:
            continue
        else:
            sortie.append(jeton)
    # Tout supprimer est possible quand p est eleve ; renvoyer "" casserait
    # le decodeur, on retombe donc sur un jeton unique.
    return "".join(sortie) or rng.choice(alpha)


def croiser(a: str, b: str, rng: random.Random) -> tuple[str, str]:
    """Croisement en un point : coupe les deux listes et echange les queues."""
    ja, jb = jetons(a), jetons(b)
    i = rng.randrange(1, max(2, len(ja)))
    j = rng.randrange(1, max(2, len(jb)))
    e1 = "".join(ja[:i] + jb[j:])
    e2 = "".join(jb[:j] + ja[i:])
    return (e1 or a), (e2 or b)


def tournoi(
    pop: list[str],
    aptitudes: list[float],
    rng: random.Random,
    k: int = 3,
) -> str:
    """Tire k individus au sort et renvoie le meilleur des k.

    Prefere a la roulette pour trois raisons pratiques : il ne depend pas
    de l'echelle absolue des notes, la pression selective se regle par le
    seul k, et il n'y a pas de division par zero quand toutes les aptitudes
    sont egales.
    """
    k = min(k, len(pop))
    indices = rng.sample(range(len(pop)), k)
    meilleur = max(indices, key=lambda i: aptitudes[i])
    return pop[meilleur]
