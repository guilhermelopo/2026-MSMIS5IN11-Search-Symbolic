"""Fonctions objectif : QED, SA score et leur combinaison.

Deux notes aux directions opposees :

  QED  0..1   plus haut est meilleur   "a-t-elle un profil de medicament oral ?"
  SA   1..10  plus BAS est meilleur    "peut-on la fabriquer ?"

Pour les combiner, le SA est inverse et normalise vers 0..1 par
`normaliser_sa`. Toute l'arithmetique des objectifs vit ici ; aucun autre
module ne doit inverser une echelle de son cote.

Le budget d'evaluations est compte ici aussi, parce que c'est le seul
endroit par lequel toutes les evaluations passent. Comparer deux methodes
avec des budgets differents devient ainsi impossible par inadvertance.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from functools import lru_cache

from rdkit.Chem import QED, RDConfig

from . import representation as rep

# --------------------------------------------------------------------------
# SA score : il vit dans Contrib, pas dans le noyau de RDKit
# --------------------------------------------------------------------------

_sascorer = None


def _charger_sascorer():
    """Importe sascorer en rustinant sys.path une seule fois.

    Le SA score ne fait pas partie du RDKit installe par pip ; il arrive
    dans le dossier des contributions. Cette rustine est isolee ici a
    dessein -- disseminee dans le code, c'est la premiere chose qui casse
    sur la machine de quelqu'un d'autre.
    """
    global _sascorer
    if _sascorer is None:
        chemin = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if chemin not in sys.path:
            sys.path.append(chemin)
        import sascorer  # noqa: PLC0415 - import tardif volontaire

        _sascorer = sascorer
    return _sascorer


# --------------------------------------------------------------------------
# Budget
# --------------------------------------------------------------------------

class Budget:
    """Compte les evaluations d'aptitude et s'arrete quand la limite tombe.

    C'est la monnaie du projet : comparer le GA aux lignes de base n'est
    honnete qu'a nombre d'evaluations egal. Les trois methodes recoivent un
    Budget et s'arretent des que `epuise` devient vrai.
    """

    def __init__(self, total: int) -> None:
        self.total = total
        self.utilisees = 0

    def depenser(self, n: int = 1) -> None:
        self.utilisees += n

    @property
    def restantes(self) -> int:
        return max(0, self.total - self.utilisees)

    @property
    def epuise(self) -> bool:
        return self.utilisees >= self.total

    def __repr__(self) -> str:
        return f"Budget({self.utilisees}/{self.total})"


# --------------------------------------------------------------------------
# Notes
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Notes:
    """Les notes d'une molecule. `valide` faux signifie chaine illisible."""

    qed: float
    sa: float
    sa_norm: float
    valide: bool = True

    @property
    def objectifs(self) -> tuple[float, float]:
        """Couple a maximiser, utilise par le module pareto."""
        return (self.qed, self.sa_norm)


NOTES_INVALIDES = Notes(qed=0.0, sa=10.0, sa_norm=0.0, valide=False)


def normaliser_sa(sa: float) -> float:
    """SA 1..10 (a minimiser) -> 0..1 (a maximiser).

    SA 1  (trivial)   -> 1.0
    SA 10 (cauchemar) -> 0.0
    """
    return max(0.0, min(1.0, (10.0 - sa) / 9.0))


@lru_cache(maxsize=200_000)
def _notes_par_smiles(smiles_canonique: str) -> Notes:
    """Coeur du calcul, avec cache sur le SMILES canonique.

    Le cache porte sur le canonique et non sur la chaine SELFIES parce que
    plusieurs chaines decrivent la meme molecule -- le gain est important
    apres quelques generations, quand la population commence a repeter.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles_canonique)
    # `MolFromSmiles("")` renvoie une molecule VIDE, pas None. Sans le test
    # du nombre d'atomes, elle atteint le QED et le sascorer, qui renvoient
    # None -- et `float(None)` leve TypeError.
    if mol is None or mol.GetNumAtoms() == 0:
        return NOTES_INVALIDES

    try:
        sa_brut = _charger_sascorer().calculateScore(mol)
        qed_brut = QED.qed(mol)
    except Exception:  # noqa: BLE001 - molecule exotique refusee par les scores
        return NOTES_INVALIDES

    # Les deux peuvent renvoyer None sur une molecule degeneree, sans lever.
    if sa_brut is None or qed_brut is None:
        return NOTES_INVALIDES

    sa = float(sa_brut)
    return Notes(qed=float(qed_brut), sa=sa, sa_norm=normaliser_sa(sa))


def notes(chaine: str, budget: Budget | None = None) -> Notes:
    """Notes d'une chaine SELFIES. Depense une unite du budget."""
    if budget is not None:
        budget.depenser()
    can = rep.canonique(chaine)
    if can is None:
        return NOTES_INVALIDES
    return _notes_par_smiles(can)


# --------------------------------------------------------------------------
# Scalarisation
# --------------------------------------------------------------------------

def aptitude(chaine: str, poids_qed: float = 0.5, budget: Budget | None = None) -> float:
    """Somme ponderee des deux objectifs : le MVP.

    poids_qed = 1.0  -> QED seul (le plan B du catalogue)
    poids_qed = 0.0  -> facilite de synthese seule
    poids_qed = 0.5  -> compromis neutre

    La valeur du poids est arbitraire, et c'est precisement le defaut de la
    somme ponderee. La balayer de 0 a 1 dans le rapport et comparer au
    front de Pareto du module pareto.
    """
    n = notes(chaine, budget)
    return poids_qed * n.qed + (1.0 - poids_qed) * n.sa_norm


def vider_cache() -> None:
    """Vide les caches. Utile entre deux executions de graines differentes."""
    _notes_par_smiles.cache_clear()
    rep.vers_mol.cache_clear()
    rep.decodage_ok.cache_clear()
