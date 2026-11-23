"""Lignes de base : recherche aleatoire et hill climbing.

Sans elles, une courbe d'aptitude qui monte ne prouve rien -- n'importe
quelle methode qui garde son meilleur resultat produit une courbe qui monte.

Les deux renvoient le meme `Resultat` que le module ga et depensent le meme
`Budget`, de sorte que les trois methodes sont comparables sans ajustement.
L'axe horizontal des graphiques doit etre `historique_evaluations`, jamais
le nombre d'iterations : une iteration de methodes differentes ne coute pas
la meme chose.
"""

from __future__ import annotations

import random

from . import objectives as obj
from . import representation as rep
from .ga import Archive, Resultat, Trouvaille


def _enregistrer(
    archive: Archive,
    chaine: str,
    n: obj.Notes,
    apt: float,
    pas: int,
    methode: str,
    graine: int,
) -> None:
    if not n.valide:
        return
    smi = rep.canonique(chaine)
    if smi is None:
        return
    archive.enregistrer(
        Trouvaille(
            smiles=smi, selfies=chaine, qed=n.qed, sa=n.sa, sa_norm=n.sa_norm,
            aptitude=apt, generation=pas, methode=methode, graine=graine,
        )
    )


def recherche_aleatoire(
    graine: int,
    budget: obj.Budget,
    poids_qed: float = 0.5,
    pas_enregistrement: int = 100,
) -> Resultat:
    """Le plancher. Tire des chaines au hasard et garde la meilleure.

    Si l'algorithme genetique fait match nul avec ceci, il n'evolue rien.
    La cause habituelle est un bug : aptitude calculee mais non utilisee
    dans la selection, ou taux de mutation si eleve qu'il detruit tout
    heritage.
    """
    rng = random.Random(graine)
    archive = Archive()
    hist_meilleur: list[float] = []
    hist_moyenne: list[float] = []
    hist_eval: list[int] = []
    meilleur = 0.0
    somme, compte = 0.0, 0
    pas = 0

    while not budget.epuise:
        chaine = rep.aleatoire(rng)
        n = obj.notes(chaine, budget)
        apt = poids_qed * n.qed + (1 - poids_qed) * n.sa_norm
        _enregistrer(archive, chaine, n, apt, pas, "aleatoire", graine)
        meilleur = max(meilleur, apt)
        somme += apt
        compte += 1
        pas += 1
        if pas % pas_enregistrement == 0:
            hist_meilleur.append(meilleur)
            hist_moyenne.append(somme / compte)
            hist_eval.append(budget.utilisees)

    if not hist_meilleur or hist_eval[-1] != budget.utilisees:
        hist_meilleur.append(meilleur)
        hist_moyenne.append(somme / max(1, compte))
        hist_eval.append(budget.utilisees)

    return Resultat(
        methode="aleatoire", graine=graine, archive=archive,
        historique_meilleur=hist_meilleur, historique_moyenne=hist_moyenne,
        historique_evaluations=hist_eval, population_finale=[],
    )


def hill_climbing(
    graine: int,
    budget: obj.Budget,
    poids_qed: float = 0.5,
    p_mutation: float = 0.15,
    redemarrages: bool = False,
    patience: int = 500,
    pas_enregistrement: int = 100,
) -> Resultat:
    """Recherche locale gloutonne : mute et n'accepte que les ameliorations.

    Monte vite au debut puis se bloque. Arrive sur un optimum local, toute
    mutation degrade et il y reste -- c'est le contraste avec la population
    du GA, qui couvre plusieurs sommets a la fois.

    `redemarrages=True` repart d'une chaine neuve apres `patience` essais
    sans amelioration. C'est une ligne de base plus forte et plus honnete ;
    executer les deux et montrer la difference.
    """
    from . import operators as ops

    rng = random.Random(graine)
    archive = Archive()
    hist_meilleur: list[float] = []
    hist_moyenne: list[float] = []
    hist_eval: list[int] = []

    courant = rep.aleatoire(rng)
    n_courant = obj.notes(courant, budget)
    apt_courant = poids_qed * n_courant.qed + (1 - poids_qed) * n_courant.sa_norm
    _enregistrer(archive, courant, n_courant, apt_courant, 0, "hill", graine)

    meilleur = apt_courant
    somme, compte = apt_courant, 1
    pas, sans_amelioration = 0, 0

    while not budget.epuise:
        candidat = ops.muter(courant, rng, p_mutation)
        n = obj.notes(candidat, budget)
        apt = poids_qed * n.qed + (1 - poids_qed) * n.sa_norm
        _enregistrer(archive, candidat, n, apt, pas, "hill", graine)

        if apt > apt_courant:
            courant, apt_courant = candidat, apt
            sans_amelioration = 0
        else:
            sans_amelioration += 1
            if redemarrages and sans_amelioration >= patience:
                courant = rep.aleatoire(rng)
                n_courant = obj.notes(courant, budget)
                apt_courant = (poids_qed * n_courant.qed
                               + (1 - poids_qed) * n_courant.sa_norm)
                sans_amelioration = 0

        meilleur = max(meilleur, apt)
        somme += apt
        compte += 1
        pas += 1
        if pas % pas_enregistrement == 0:
            hist_meilleur.append(meilleur)
            hist_moyenne.append(somme / compte)
            hist_eval.append(budget.utilisees)

    if not hist_meilleur or hist_eval[-1] != budget.utilisees:
        hist_meilleur.append(meilleur)
        hist_moyenne.append(somme / max(1, compte))
        hist_eval.append(budget.utilisees)

    return Resultat(
        methode="hill" + ("_redemarrages" if redemarrages else ""),
        graine=graine, archive=archive,
        historique_meilleur=hist_meilleur, historique_moyenne=hist_moyenne,
        historique_evaluations=hist_eval, population_finale=[courant],
    )
