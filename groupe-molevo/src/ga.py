"""L'algorithme genetique.

La boucle principale est volontairement courte ; toute la complexite reelle
vit dans les operateurs et les fonctions objectif. Deux choix de conception
meritent d'etre signales :

1. L'*archive des trouvailles* conserve toute molecule distincte evaluee, et
   pas seulement la population finale. Cela regle l'essentiel du probleme
   d'effondrement de la diversite en trois lignes : meme si la population
   finit en clones, l'archive preserve tout ce qui a ete rencontre.

2. Le budget d'evaluations est partage avec les lignes de base, de sorte que
   la comparaison est juste par construction et non par discipline.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field

from . import objectives as obj
from . import operators as ops
from . import representation as rep


@dataclass
class Trouvaille:
    """Une molecule distincte rencontree a un moment de la recherche."""

    smiles: str
    selfies: str
    qed: float
    sa: float
    sa_norm: float
    aptitude: float
    generation: int
    methode: str
    graine: int

    def en_dict(self) -> dict:
        return self.__dict__.copy()


@dataclass
class Archive:
    """Collectionne les molecules distinctes, indexees par SMILES canonique.

    Indexer par SMILES canonique est essentiel : la meme molecule apparait
    sous plusieurs chaines SELFIES, et indexer par la chaine brute gonflerait
    le compte de « molecules distinctes » sans que rien de neuf ait ete
    decouvert.
    """

    items: dict[str, Trouvaille] = field(default_factory=dict)

    def enregistrer(self, trouvaille: Trouvaille) -> bool:
        """Renvoie True si la molecule etait inedite."""
        if trouvaille.smiles in self.items:
            return False
        self.items[trouvaille.smiles] = trouvaille
        return True

    def meilleures(self, n: int = 10) -> list[Trouvaille]:
        return sorted(self.items.values(), key=lambda t: t.aptitude, reverse=True)[:n]

    def __len__(self) -> int:
        return len(self.items)


@dataclass
class Resultat:
    """Ce qu'une execution renvoie."""

    methode: str
    graine: int
    archive: Archive
    historique_meilleur: list[float]      # meilleure aptitude jusqu'a chaque generation
    historique_moyenne: list[float]
    historique_evaluations: list[int]     # evaluations depensees a chaque generation
    population_finale: list[str]

    @property
    def meilleure_aptitude(self) -> float:
        return self.historique_meilleur[-1] if self.historique_meilleur else 0.0


def _evaluer(
    pop: list[str],
    poids_qed: float,
    budget: obj.Budget,
) -> list[obj.Notes]:
    return [obj.notes(c, budget) for c in pop]


def evoluer(
    graine: int,
    taille_pop: int = 100,
    generations: int = 100,
    elites: int = 2,
    k_tournoi: int = 3,
    p_croisement: float = 0.8,
    p_mutation: float = 0.15,
    poids_qed: float = 0.5,
    graines_smiles: list[str] | None = None,
    budget: obj.Budget | None = None,
) -> Resultat:
    """Execute l'algorithme genetique et renvoie le Resultat.

    Si `budget` n'est pas fourni, il en est cree un de taille_pop *
    generations, ce qui est exactement le cout d'une execution complete.
    """
    rng = random.Random(graine)
    if budget is None:
        budget = obj.Budget(taille_pop * generations)

    pop = rep.population_initiale(rng, taille_pop, graines_smiles)
    archive = Archive()
    hist_meilleur: list[float] = []
    hist_moyenne: list[float] = []
    hist_eval: list[int] = []
    meilleur_global = 0.0

    for g in range(generations):
        if budget.epuise:
            break

        notes_pop = _evaluer(pop, poids_qed, budget)
        apt = [poids_qed * n.qed + (1 - poids_qed) * n.sa_norm for n in notes_pop]

        for chaine, n, a in zip(pop, notes_pop, apt):
            if not n.valide:
                continue
            smi = rep.canonique(chaine)
            if smi is None:
                continue
            archive.enregistrer(
                Trouvaille(
                    smiles=smi, selfies=chaine, qed=n.qed, sa=n.sa,
                    sa_norm=n.sa_norm, aptitude=a, generation=g,
                    methode="ga", graine=graine,
                )
            )

        meilleur_global = max(meilleur_global, max(apt))
        hist_meilleur.append(meilleur_global)
        hist_moyenne.append(sum(apt) / len(apt))
        hist_eval.append(budget.utilisees)

        # Elitisme : les e meilleurs passent intacts. Sans cela, une mutation
        # malchanceuse peut detruire la meilleure trouvaille et la courbe peut
        # DESCENDRE.
        ordre = sorted(range(len(pop)), key=lambda i: apt[i], reverse=True)
        nouvelle = [pop[i] for i in ordre[:elites]]

        while len(nouvelle) < taille_pop:
            p1 = ops.tournoi(pop, apt, rng, k_tournoi)
            p2 = ops.tournoi(pop, apt, rng, k_tournoi)
            if rng.random() < p_croisement:
                e1, e2 = ops.croiser(p1, p2, rng)
            else:
                e1, e2 = p1, p2
            nouvelle.append(ops.muter(e1, rng, p_mutation))
            if len(nouvelle) < taille_pop:
                nouvelle.append(ops.muter(e2, rng, p_mutation))

        pop = nouvelle

    return Resultat(
        methode="ga",
        graine=graine,
        archive=archive,
        historique_meilleur=hist_meilleur,
        historique_moyenne=hist_moyenne,
        historique_evaluations=hist_eval,
        population_finale=pop,
    )
