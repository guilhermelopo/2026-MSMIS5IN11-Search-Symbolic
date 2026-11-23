"""Representation des molecules sous forme de chaines SELFIES.

Ce module isole tout ce qui depend du format textuel des molecules. Le
reste du projet ne manipule jamais de SMILES ni de SELFIES directement :
il passe par ici. Changer de representation (par exemple pour un
algorithme genetique sur graphes) ne touche donc qu'un seul fichier.

La propriete centrale de SELFIES est que *toute* sequence de jetons se
decode sans erreur. Le test
tests/test_representation.py::test_robustesse_alphabet le verifie.
"""

from __future__ import annotations

import random
from functools import lru_cache

import selfies as sf
from rdkit import Chem, RDLogger

# RDKit emet des avertissements en masse a la lecture de molecules
# exotiques, et l'algorithme genetique en produit beaucoup. Ce ne sont pas
# des erreurs ; sans cette ligne, une execution depasse les dizaines de
# milliers de lignes de sortie.
RDLogger.DisableLog("rdApp.*")


# --------------------------------------------------------------------------
# Alphabet
# --------------------------------------------------------------------------

_ALPHABET: tuple[str, ...] | None = None


def alphabet() -> tuple[str, ...]:
    """Jetons que la bibliotheque garantit surs, dans un ordre fixe.

    L'ordre fixe compte pour la reproductibilite : `rng.choice` sur un
    tuple trie donne toujours le meme resultat pour la meme graine.
    """
    global _ALPHABET
    if _ALPHABET is None:
        _ALPHABET = tuple(sorted(sf.get_semantic_robust_alphabet()))
    return _ALPHABET


def jetons(chaine: str) -> list[str]:
    """Decoupe une chaine SELFIES en sa liste de jetons."""
    return list(sf.split_selfies(chaine))


# --------------------------------------------------------------------------
# Conversions
# --------------------------------------------------------------------------

def vers_smiles(chaine: str) -> str:
    """SELFIES -> SMILES."""
    return sf.decoder(chaine)


def depuis_smiles(smiles: str) -> str:
    """SMILES -> SELFIES."""
    return sf.encoder(smiles)


@lru_cache(maxsize=200_000)
def decodage_ok(chaine: str) -> bool:
    """La chaine se decode-t-elle sans erreur ? (molecule vide incluse)

    C'est la propriete exacte que SELFIES garantit : le decodage n'echoue
    jamais. Utiliser `valide` quand il faut une molecule avec des atomes.
    """
    try:
        return Chem.MolFromSmiles(sf.decoder(chaine)) is not None
    except Exception:  # noqa: BLE001 - jeton hors de l'alphabet robuste
        return False


@lru_cache(maxsize=200_000)
def vers_mol(chaine: str):
    """SELFIES -> objet Mol de RDKit avec au moins un atome, sinon None.

    PIEGE : `Chem.MolFromSmiles("")` ne renvoie PAS None -- il renvoie une
    molecule VIDE, sans atome. Une chaine SELFIES degeneree (un jeton de
    cycle ou d'indice seul, par exemple `[Ring1]`) se decode en SMILES
    vide. Sans le test du nombre d'atomes ci-dessous, cette molecule vide
    traverse le `is None` et atteint le QED et le SA score, qui renvoient
    alors des valeurs indefinies ou None -- c'etait l'origine du
    `TypeError: float() argument must be ... not 'NoneType'`.
    """
    try:
        mol = Chem.MolFromSmiles(sf.decoder(chaine))
    except Exception:  # noqa: BLE001
        return None
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    return mol


def canonique(chaine: str) -> str | None:
    """SMILES canonique, la bonne cle pour comparer des molecules.

    Deux chaines SELFIES differentes peuvent decrire la meme molecule.
    Comparer les chaines brutes compte des doublons comme distincts.
    """
    mol = vers_mol(chaine)
    if mol is None:
        return None
    can = Chem.MolToSmiles(mol)
    return can or None          # ne renvoie jamais la chaine vide


def valide(chaine: str) -> bool:
    """La chaine decrit une molecule ayant au moins un atome."""
    return vers_mol(chaine) is not None


# --------------------------------------------------------------------------
# Generation
# --------------------------------------------------------------------------

def aleatoire(rng: random.Random, n_min: int = 6, n_max: int = 25) -> str:
    """Une chaine SELFIES aleatoire de n_min a n_max jetons."""
    alpha = alphabet()
    n = rng.randint(n_min, n_max)
    return "".join(rng.choice(alpha) for _ in range(n))


def population_initiale(
    rng: random.Random,
    n: int,
    graines: list[str] | None = None,
) -> list[str]:
    """Construit la population de depart.

    `graines` est une liste optionnelle de SMILES connus. Partir de
    molecules reelles converge souvent plus vite que partir de rien, mais
    biaise le resultat vers leur voisinage -- comparer les deux dans le
    rapport est une experience peu couteuse.
    """
    pop: list[str] = []
    if graines:
        for smi in graines:
            try:
                pop.append(depuis_smiles(smi))
            except Exception:  # noqa: BLE001 - SMILES d'entree invalide
                continue
    while len(pop) < n:
        pop.append(aleatoire(rng))
    return pop[:n]
