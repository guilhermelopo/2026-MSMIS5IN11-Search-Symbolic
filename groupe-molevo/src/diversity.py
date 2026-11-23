"""Mesures de diversite.

L'algorithme genetique converge, et ce n'est pas seulement une bonne chose :
apres quelques dizaines de generations, la population devient des copies du
meme gagnant. La courbe d'aptitude est excellente et le resultat cesse
d'etre une recherche.

Pour la decouverte de medicaments, c'est un echec deguise : l'objectif n'a
jamais ete *une* molecule, mais un ensemble de candidats differents, car la
plupart echoueront plus tard pour des raisons qu'aucune de ces formules ne
prevoit.

Reporter la courbe d'aptitude a cote de celle de diversite est l'une des
figures les plus fortes que ce projet puisse produire.
"""

from __future__ import annotations

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold


def empreinte(smiles: str, rayon: int = 2, bits: int = 2048):
    """Empreinte moleculaire de Morgan (ECFP4 pour rayon=2).

    Un vecteur de bits ou chaque bit indique la presence d'une sous-structure
    centree sur un atome, jusqu'a `rayon` liaisons de distance.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=rayon, nBits=bits)


def diversite_interne(liste_smiles: list[str]) -> float:
    """1 moins la similarite moyenne de Tanimoto sur toutes les paires.

    1.0 = tout different, 0.0 = des clones. Cout quadratique en nombre de
    molecules ; pour de grandes populations, echantillonner.
    """
    emp = [e for e in (empreinte(s) for s in liste_smiles) if e is not None]
    if len(emp) < 2:
        return 0.0
    sims: list[float] = []
    for i in range(len(emp)):
        sims.extend(DataStructs.BulkTanimotoSimilarity(emp[i], emp[i + 1:]))
    return 1.0 - sum(sims) / len(sims)


def squelette(smiles: str) -> str | None:
    """Squelette de Bemis-Murcko : la molecule sans ses chaines laterales."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        noyau = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(noyau)
    except Exception:  # noqa: BLE001 - les molecules sans cycle peuvent echouer
        return None


def compter_squelettes(liste_smiles: list[str]) -> int:
    """Combien de squelettes distincts existent dans l'ensemble.

    Mesure plus parlante pour un chimiste que le Tanimoto moyen : cent
    molecules reparties sur trois squelettes, c'est peu de diversite reelle,
    meme si la similarite moyenne ne parait pas si elevee.
    """
    return len({s for s in (squelette(x) for x in liste_smiles) if s})


def resume(liste_smiles: list[str]) -> dict:
    """Les trois mesures d'un coup, pretes a ecrire dans le registre."""
    distinctes = sorted(set(liste_smiles))
    return {
        "n_molecules": len(liste_smiles),
        "n_distinctes": len(distinctes),
        "diversite_interne": diversite_interne(distinctes),
        "n_squelettes": compter_squelettes(distinctes),
    }
