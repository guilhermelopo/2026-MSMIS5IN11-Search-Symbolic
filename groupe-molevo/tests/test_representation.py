"""Tests de la representation SELFIES.

Le test `test_robustesse_alphabet` est le plus important du projet : il
verifie la propriete sur laquelle tout le reste repose. S'il echoue,
l'algorithme genetique perd sa justification.
"""

import random

import pytest
from rdkit import Chem

from src import representation as rep

# Quelques medicaments connus, pour l'aller-retour.
CONNUS = [
    "CC(=O)Oc1ccccc1C(=O)O",            # aspirine
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",       # ibuprofene
    "CC(=O)Nc1ccc(O)cc1",               # paracetamol
    "Cn1cnc2c1c(=O)n(C)c(=O)n2C",       # cafeine
    "CCO",                              # ethanol
    "c1ccccc1",                         # benzene
]


def canoniser(smiles: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


# --------------------------------------------------------------------------
# Aller-retour
# --------------------------------------------------------------------------

@pytest.mark.parametrize("smiles", CONNUS)
def test_aller_retour_preserve_la_molecule(smiles):
    chaine = rep.depuis_smiles(smiles)
    retour = rep.vers_smiles(chaine)
    assert canoniser(retour) == canoniser(smiles)


def test_canonique_concorde_avec_rdkit():
    chaine = rep.depuis_smiles("CC(=O)Oc1ccccc1C(=O)O")
    assert rep.canonique(chaine) == canoniser("CC(=O)Oc1ccccc1C(=O)O")


def test_deux_chaines_peuvent_donner_la_meme_molecule():
    # C'est pourquoi le cache et l'archive indexent par SMILES canonique et
    # non par la chaine SELFIES.
    a = rep.depuis_smiles("CCO")
    b = rep.depuis_smiles("OCC")
    assert rep.canonique(a) == rep.canonique(b)


# --------------------------------------------------------------------------
# La propriete centrale
# --------------------------------------------------------------------------

def test_robustesse_alphabet():
    """La promesse de SELFIES : toute sequence de jetons SE DECODE.

    C'est exactement la propriete que la bibliotheque garantit, et c'est la
    seule dont l'algorithme genetique a besoin : aucune mutation ne produit
    une chaine que le decodeur refuse.
    """
    rng = random.Random(0)
    echecs = []
    for _ in range(2000):
        chaine = rep.aleatoire(rng, n_min=1, n_max=30)
        if not rep.decodage_ok(chaine):
            echecs.append(chaine)
    assert not echecs, f"{len(echecs)} chaines ne se decodent pas, ex. : {echecs[:3]}"


def test_chaines_longues_donnent_des_molecules_avec_atomes():
    """Le decodage reussit toujours, mais peut donner une molecule VIDE.

    Une chaine d'un seul jeton peut etre degeneree (`[Ring1]` n'a pas
    d'atome). Avec quelques jetons cela n'arrive quasiment plus -- et c'est
    pourquoi la population initiale n'utilise jamais n_min=1.
    """
    rng = random.Random(0)
    sans_atome = [c for _ in range(500)
                  if not rep.valide(c := rep.aleatoire(rng, n_min=6, n_max=30))]
    assert len(sans_atome) < 25, (
        f"{len(sans_atome)}/500 chaines sans atome, ex. : {sans_atome[:3]}")


def test_molecule_vide_est_rejetee():
    """Le piege qui a provoque un TypeError en execution reelle.

    `Chem.MolFromSmiles("")` renvoie une molecule VIDE, pas None. Si elle
    passe, le QED et le SA score recoivent zero atome et renvoient None, et
    `float(None)` leve TypeError. `vers_mol` doit renvoyer None.
    """
    vide = Chem.MolFromSmiles("")
    assert vide is not None, "RDKit a change : MolFromSmiles('') renvoie None"
    assert vide.GetNumAtoms() == 0

    # `[Ring1]` seul : un jeton de fermeture de cycle, sans aucun atome.
    assert rep.decodage_ok("[Ring1]")      # se decode, comme promis
    assert not rep.valide("[Ring1]")       # mais n'a pas d'atome
    assert rep.vers_mol("[Ring1]") is None
    assert rep.canonique("[Ring1]") is None   # jamais la chaine vide


def test_alphabet_non_vide_et_stable():
    a, b = rep.alphabet(), rep.alphabet()
    assert len(a) > 10
    assert a == b          # l'ordre fixe garantit la reproductibilite


def test_jetons_recomposent_la_chaine():
    rng = random.Random(1)
    chaine = rep.aleatoire(rng)
    assert "".join(rep.jetons(chaine)) == chaine


# --------------------------------------------------------------------------
# Generation
# --------------------------------------------------------------------------

def test_aleatoire_respecte_la_graine():
    assert rep.aleatoire(random.Random(42)) == rep.aleatoire(random.Random(42))


def test_population_initiale_a_la_taille_demandee():
    pop = rep.population_initiale(random.Random(0), 25)
    assert len(pop) == 25
    assert all(rep.valide(s) for s in pop)


def test_population_initiale_utilise_les_graines():
    pop = rep.population_initiale(random.Random(0), 10, graines=CONNUS[:3])
    canoniques = {rep.canonique(s) for s in pop}
    assert canoniser(CONNUS[0]) in canoniques
    assert len(pop) == 10


def test_population_initiale_ignore_une_graine_invalide():
    pop = rep.population_initiale(
        random.Random(0), 5, graines=["ceci n'est pas un smiles"])
    assert len(pop) == 5
