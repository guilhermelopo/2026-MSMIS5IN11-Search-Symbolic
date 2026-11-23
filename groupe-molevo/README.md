# C3 — Conception évolutionnaire de molécules candidates

Algorithme génétique sur chaînes SELFIES pour explorer l'espace chimique à la
recherche de molécules au profil de médicament oral (QED) qui restent
plausibles à synthétiser (SA score).

## Groupe

| Membre | Rôle | Contribution principale |
|---|---|---|
| [@guilhermelopo](https://github.com/guilhermelopo) | Représentation et objectifs | `representation.py`, `objectives.py`, `diversity.py`, `verifier.py` — la couche « quoi » : comment une molécule est écrite et comment elle est notée |
| [@lucaslapenda-cpu](https://github.com/lucaslapenda-cpu) | Recherche et évaluation | `operators.py`, `ga.py`, `baselines.py`, `pareto.py` — la couche « comment » : la boucle évolutive, les lignes de base et le multiobjectif |

Chaque membre écrit les tests de ses propres modules et relit le code de
l'autre. Le `verifier.py` est écrit par le membre qui n'a **pas** écrit
l'algorithme génétique : c'est ce qui rend la vérification indépendante.

## Sujet

- **Identifiant du catalogue :** C3 — Conception évolutionnaire de molécules candidates
- **Problématique :** l'espace des molécules d'apparence médicamenteuse est
  estimé à 10⁶⁰ structures. Il n'est pas énumérable, et il n'y a pas de
  gradient à suivre, car une molécule est un objet discret. Reste la recherche
  heuristique.
- **MVP :** algorithme génétique sur SELFIES (mutations et croisements toujours
  valides), optimisation du QED et de la facilité de synthèse avec RDKit.

## Architecture

```
representation.py ──┐
                    ├──► ga.py ────────┐
objectives.py ──────┤                  ├──► registro.jsonl ──► diversity.py
                    └──► baselines.py ─┘                       pareto.py
operators.py ───────┘                                          verifier.py
                                                               notebooks/
```

La recherche ne fait qu'**écrire** dans le registre ; l'analyse ne fait que le
**lire**. Refaire une figure ne coûte pas une exécution complète, et le
notebook n'importe rien de la recherche.

| Module | Responsabilité |
|---|---|
| `src/representation.py` | SELFIES ↔ SMILES, alphabet robuste, génération de chaînes |
| `src/objectives.py` | QED, SA score, normalisation, somme pondérée, budget |
| `src/operators.py` | tournoi, croisement en un point, mutation |
| `src/ga.py` | boucle évolutive, élitisme, archive des trouvailles |
| `src/baselines.py` | recherche aléatoire et hill climbing, même budget |
| `src/diversity.py` | Tanimoto, diversité interne, squelettes de Murcko |
| `src/pareto.py` | dominance, fronts, couches, hypervolume (Python pur) |
| `src/verifier.py` | vérificateur indépendant des solutions |
| `src/run.py` | ligne de commande : exécute tout et écrit le registre |
| `medir.py` | produit toutes les mesures du rapport en une fois |

### Deux décisions à expliquer à la soutenance

**Pourquoi SELFIES et pas SMILES.** Une mutation d'un caractère dans un SMILES
produit presque toujours une chaîne qui ne correspond à aucune molécule :
parenthèse ouverte, cycle qui ne se ferme pas, carbone à cinq liaisons. En
SELFIES, toute séquence de jetons se décode en une molécule valide, parce que
la grammaire rabaisse toute demande impossible au lieu de la rejeter. Le
croisement et la mutation deviennent des opérations triviales sur des listes,
sans rejet.
Vérifié par `tests/test_representation.py::test_robustez_do_alfabeto`.

**Pourquoi le budget d'évaluations est un objet.** Comparer l'algorithme
génétique aux lignes de base n'est honnête qu'avec le même nombre d'appels à
la fonction d'aptitude. La classe `Orcamento` est partagée par les trois
méthodes : la comparaison est juste par construction et non par discipline.
L'axe horizontal des graphiques est en **évaluations**, jamais en
**générations**.

## Installation

Python 3.10 ou plus récent.

```bash
python -m venv .venv
.venv/Scripts/activate        # Windows ; sous Linux/macOS : source .venv/bin/activate
pip install -r requirements.txt
```

> **Windows.** RDKit charge des DLL natives. Si l'import échoue avec
> *« DLL load failed »*, deux causes fréquentes : le chemin du projet est trop
> long (gardez-le court, par exemple `C:\dev\c3`), ou une politique de
> contrôle d'application bloque les bibliothèques non signées. Dans ce second
> cas, utilisez WSL, Conda ou Google Colab.

## Utilisation

```bash
# exécution rapide, pour vérifier que tout fonctionne
python -m src.run --sementes 1 --pop 30 --geracoes 10

# exécution du rapport : 10 graines, 10 000 évaluations par méthode
python -m src.run --sementes 10 --pop 100 --geracoes 100

# plan B du catalogue : un seul objectif
python -m src.run --sementes 10 --peso-qed 1.0
```

Sortie dans `resultados/` : `registro.jsonl` (une molécule distincte par ligne,
avec ses notes, la méthode, la graine et la génération) et `resumo.json`
(courbes agrégées, front de Pareto, hypervolume, diversité).

### Vérificateur indépendant

```bash
python -m src.verifier resultados/registro.jsonl
```

Relit le registre et revérifie tout par un chemin de code propre, sans
réutiliser le cache ni les fonctions employées pendant la recherche : SMILES
lisibles et canoniques, aller-retour par SELFIES, QED et SA recalculés, absence
de doublons, front de Pareto réellement non dominé. Code de sortie 0 s'il n'y a
aucun problème.

## Tests

```bash
pytest -q
```

79 tests répartis en cinq fichiers.

| Fichier | Ce qu'il prouve |
|---|---|
| `test_pareto.py` | dominance, fronts, couches, hypervolume, avec des réponses calculées à la main. Ne dépend pas de RDKit. |
| `test_representation.py` | aller-retour SMILES↔SELFIES ; **toute** chaîne aléatoire est valide |
| `test_operators.py` | les opérateurs ne produisent jamais de chaîne vide ou invalide ; la taille varie dans les deux sens |
| `test_objectives.py` | sens de la normalisation du SA ; budget correctement compté |
| `test_ga.py` | élitisme (le meilleur ne régresse jamais), reproductibilité par graine, le vérificateur détecte une note falsifiée |

### État de la vérification

| | Tests | Résultat |
|---|---:|---|
| `test_pareto.py` | 20 | **20 passent** (exécuté) |
| les quatre autres | 59 | **pas encore exécutés** — RDKit ne se charge pas sur la machine où le squelette a été écrit |

**À COMPLÉTER :** exécuter `pytest -q` sur une machine où RDKit fonctionne et
coller la sortie ici. Certains de ces 59 tests peuvent échouer : le code a été
écrit sans être exécuté.

## Résultats et évaluation

**À COMPLÉTER.** Aucune molécule n'a encore été générée ; aucun QED ni SA n'a
été calculé. `python medir.py` produit en une fois toutes les mesures
attendues :

| Mesure | Alimente |
|---|---|
| Validité sous mutation : SMILES contre SELFIES, 2 000 mutations chacun | la justification du projet |
| QED et SA de six médicaments connus | l'échelle réelle des deux notes |
| GA contre recherche aléatoire et hill climbing, 5 graines, même budget | le résultat principal |
| Diversité interne et nombre de squelettes | la section diversité |
| Front de Pareto et hypervolume | le multiobjectif |
| Les cinq molécules au QED le plus élevé | la section limites |

Protocole attendu : budget en évaluations, graines fixées et consignées,
moyenne et écart-type sur au moins cinq graines, comparaison à une ligne de
base simple.

## Limites et plan B

- **QED et SA sont des approximations statistiques, pas des prédictions.** Le
  QED mesure la ressemblance du profil physique avec celui des médicaments
  oraux connus ; le SA mesure la fréquence de fragments dans des bases de
  composés, sans simuler aucune réaction. Rien ici ne dit qu'une molécule
  fonctionne ou qu'elle est sûre.
- **Loi de Goodhart.** Optimiser le QED jusqu'à la limite produit des molécules
  à note élevée et à chimie absurde. Les médicaments approuvés dépassent
  rarement 0,8 ; un résultat à 0,95 est un signal d'alerte, pas un succès.
- **SELFIES garantit la validité, pas la plausibilité.** Les molécules sont
  bien formées (valences respectées) ; pas nécessairement stables ni
  existantes.
- **Le SA pénalise la nouveauté.** Une structure nouvelle et facile à faire
  reçoit une mauvaise note simplement parce qu'elle est rare dans les bases.
- **Effondrement de la diversité.** La population converge vers des clones.
  Atténué par l'archive des trouvailles, qui conserve tout ce qui a été
  rencontré en chemin.
- **Plan B :** objectif unique (QED), avec `--peso-qed 1.0`.

## Répartition du travail

Découpage en tranches verticales : chacun possède une moitié du pipeline de
bout en bout, mène une des deux expériences principales et présente sa moitié
de la soutenance.

| | [@guilhermelopo](https://github.com/guilhermelopo) | [@lucaslapenda-cpu](https://github.com/lucaslapenda-cpu) |
|---|---|---|
| **Modules** | `representation.py`, `objectives.py`, `diversity.py`, `verifier.py` | `operators.py`, `ga.py`, `baselines.py`, `pareto.py` |
| **Expérience** | validité sous mutation : SMILES contre SELFIES | GA contre recherche aléatoire et hill climbing |
| **Notebook** | validité, notes des médicaments connus, diversité | convergence, front de Pareto, grille des meilleures molécules |
| **Soutenance** | problème, représentation, objectifs, limites | algorithme, protocole, résultats, diversité, Pareto |
| **Relecture** | relit la couche recherche | relit la couche représentation |

En commun : `run.py`, ce README, et les slides d'ouverture et de conclusion.

Chacun doit pouvoir expliquer l'architecture, les choix, les mesures et les
limites de **tout** le projet, pas seulement de sa moitié.

## Contribution d'agents de code

Le squelette du code et ce README ont été écrits avec l'aide d'un agent de code
(Claude Code), conformément à la charte `AGENTS.md` du dépôt du cours. Les
mesures, le notebook et les slides sont produits et relus par le groupe. La
déclaration est obligatoire ; une contribution déclarée ne retire aucun point.

## Licence et références

Code sous licence MIT, comme le dépôt du cours. RDKit (BSD-3-Clause),
selfies (Apache-2.0).

1. Krenn, M. et al. (2020). Self-referencing embedded strings (SELFIES): a 100%
   robust molecular string representation. *Machine Learning: Science and
   Technology*, 1(4), 045024.
2. Bickerton, G. R. et al. (2012). Quantifying the chemical beauty of drugs.
   *Nature Chemistry*, 4, 90–98.
3. Ertl, P., Schuffenhauer, A. (2009). Estimation of synthetic accessibility
   score of drug-like molecules based on molecular complexity and fragment
   contributions. *Journal of Cheminformatics*, 1, 8.
4. Brown, N. et al. (2019). GuacaMol: benchmarking models for de novo molecular
   design. *J. Chem. Inf. Model.*, 59(3), 1096–1108.
5. Deb, K. et al. (2002). A fast and elitist multiobjective genetic algorithm:
   NSGA-II. *IEEE Trans. Evol. Comput.*, 6(2), 182–197.
6. Rogers, D., Hahn, M. (2010). Extended-connectivity fingerprints.
   *J. Chem. Inf. Model.*, 50(5), 742–754.
7. Nigam, A. et al. (2020). Augmenting genetic algorithms with deep neural
   networks for exploring the chemical space. *ICLR 2020*.
