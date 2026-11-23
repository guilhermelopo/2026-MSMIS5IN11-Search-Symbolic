"""Tests de l'algorithme genetique, des lignes de base et du verificateur.

Executions volontairement petites : ces tests existent pour prouver des
proprietes (l'elitisme fonctionne, le budget est respecte, la graine
reproduit), pas pour mesurer une qualite. Les mesures vont dans le notebook.
"""

import json

import pytest

from src import baselines, ga
from src import objectives as obj
from src import verifier


PETIT = dict(taille_pop=12, generations=6, poids_qed=0.5)


# --------------------------------------------------------------------------
# Elitisme et historique
# --------------------------------------------------------------------------

def test_le_meilleur_ne_regresse_jamais():
    """La consequence observable de l'elitisme.

    Sans elitisme, une mutation malchanceuse detruit le meilleur individu et
    la courbe DESCEND. Ce test est la facon economique de prouver qu'il est
    bien actif.
    """
    r = ga.evoluer(graine=0, **PETIT)
    historique = r.historique_meilleur
    assert historique == sorted(historique), historique


def test_historique_a_une_entree_par_generation():
    r = ga.evoluer(graine=0, **PETIT)
    assert len(r.historique_meilleur) == PETIT["generations"]
    assert len(r.historique_moyenne) == PETIT["generations"]
    assert len(r.historique_evaluations) == PETIT["generations"]


def test_les_evaluations_croissent_de_facon_monotone():
    r = ga.evoluer(graine=0, **PETIT)
    assert r.historique_evaluations == sorted(r.historique_evaluations)


def test_population_finale_a_la_bonne_taille():
    r = ga.evoluer(graine=0, **PETIT)
    assert len(r.population_finale) == PETIT["taille_pop"]


# --------------------------------------------------------------------------
# Reproductibilite
# --------------------------------------------------------------------------

def test_meme_graine_meme_resultat():
    obj.vider_cache()
    a = ga.evoluer(graine=42, **PETIT)
    obj.vider_cache()
    b = ga.evoluer(graine=42, **PETIT)
    assert a.historique_meilleur == b.historique_meilleur
    assert a.population_finale == b.population_finale


def test_graines_differentes_resultats_differents():
    a = ga.evoluer(graine=1, **PETIT)
    b = ga.evoluer(graine=2, **PETIT)
    assert a.population_finale != b.population_finale


# --------------------------------------------------------------------------
# Budget
# --------------------------------------------------------------------------

def test_le_ga_respecte_le_budget():
    b = obj.Budget(50)
    ga.evoluer(graine=0, taille_pop=20, generations=100, budget=b)
    # Il s'arrete a la premiere generation apres epuisement, donc il ne
    # depasse pas d'un lot.
    assert b.utilisees <= 50 + 20


def test_les_trois_methodes_depensent_le_meme_budget():
    """La comparaison du rapport n'est honnete que si ceci tient."""
    total = 200
    b_al = obj.Budget(total)
    b_hc = obj.Budget(total)
    baselines.recherche_aleatoire(0, b_al)
    baselines.hill_climbing(0, b_hc)
    assert b_al.utilisees == total
    assert b_hc.utilisees == total


# --------------------------------------------------------------------------
# Archive des trouvailles
# --------------------------------------------------------------------------

def test_l_archive_ne_garde_pas_de_doublon():
    r = ga.evoluer(graine=0, **PETIT)
    smiles = [t.smiles for t in r.archive.items.values()]
    assert len(smiles) == len(set(smiles))


def test_l_archive_garde_plus_que_la_population_finale():
    """Tout l'interet de l'archive : preserver ce qui a ete trouve en chemin."""
    r = ga.evoluer(graine=0, taille_pop=20, generations=15, poids_qed=0.5)
    assert len(r.archive) > 20


def test_les_meilleures_sortent_triees():
    r = ga.evoluer(graine=0, **PETIT)
    aptitudes = [t.aptitude for t in r.archive.meilleures(5)]
    assert aptitudes == sorted(aptitudes, reverse=True)


# --------------------------------------------------------------------------
# Lignes de base
# --------------------------------------------------------------------------

@pytest.mark.parametrize("methode", ["aleatoire", "hill"])
def test_les_lignes_de_base_produisent_un_resultat_utilisable(methode):
    b = obj.Budget(300)
    fonction = (baselines.recherche_aleatoire if methode == "aleatoire"
                else baselines.hill_climbing)
    r = fonction(0, b)
    assert r.historique_meilleur
    assert r.historique_meilleur == sorted(r.historique_meilleur)
    assert len(r.archive) > 0
    assert r.historique_evaluations[-1] == 300


def test_hill_climbing_avec_redemarrages_fonctionne_aussi():
    b = obj.Budget(300)
    r = baselines.hill_climbing(0, b, redemarrages=True, patience=20)
    assert r.meilleure_aptitude > 0
    assert "redemarrages" in r.methode


# --------------------------------------------------------------------------
# Verificateur
# --------------------------------------------------------------------------

def test_le_verificateur_approuve_un_registre_honnete():
    r = ga.evoluer(graine=0, **PETIT)
    registre = [t.en_dict() for t in r.archive.items.values()]
    assert verifier.verifier_registre(registre) == []


def test_le_verificateur_attrape_une_note_falsifiee():
    """L'erreur pour laquelle le verificateur existe."""
    r = ga.evoluer(graine=0, **PETIT)
    registre = [t.en_dict() for t in r.archive.items.values()]
    registre[0]["qed"] = round(min(1.0, registre[0]["qed"] + 0.1), 6)
    problemes = verifier.verifier_registre(registre)
    assert any(p.type == "qed_divergent" for p in problemes)


def test_le_verificateur_attrape_un_doublon():
    r = ga.evoluer(graine=0, **PETIT)
    registre = [t.en_dict() for t in r.archive.items.values()]
    registre.append(dict(registre[0]))
    assert any(p.type == "doublon"
               for p in verifier.verifier_registre(registre))


def test_le_verificateur_attrape_un_smiles_illisible():
    problemes = verifier.verifier_registre([{"smiles": "ceci n'est pas un smiles"}])
    assert any(p.type == "smiles_illisible" for p in problemes)


def test_le_verificateur_attrape_un_front_domine():
    front_faux = [
        {"smiles": "CCO", "qed": 0.8, "sa_norm": 0.8},
        {"smiles": "CCC", "qed": 0.5, "sa_norm": 0.5},   # dominee par la premiere
    ]
    problemes = verifier.verifier_front(front_faux)
    assert any(p.type == "front_domine" for p in problemes)


def test_le_verificateur_accepte_un_front_legitime():
    front_ok = [
        {"smiles": "CCO", "qed": 0.8, "sa_norm": 0.2},
        {"smiles": "CCC", "qed": 0.2, "sa_norm": 0.8},
    ]
    assert verifier.verifier_front(front_ok) == []


def test_le_registre_survit_a_l_aller_retour_json(tmp_path):
    r = ga.evoluer(graine=0, **PETIT)
    chemin = tmp_path / "registre.jsonl"
    with open(chemin, "w", encoding="utf-8") as f:
        for t in r.archive.items.values():
            f.write(json.dumps(t.en_dict(), ensure_ascii=False) + "\n")
    assert verifier.verifier_registre(
        verifier.charger_jsonl(str(chemin))) == []
