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
objectives.py ──────┤                  ├──► registre.jsonl ──► diversity.py
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
la fonction d'aptitude. La classe `Budget` est partagée par les trois
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
python -m src.run --graines 1 --pop 30 --generations 10

# exécution du rapport : 10 graines, 10 000 évaluations par méthode
python -m src.run --graines 10 --pop 100 --generations 100

# plan B du catalogue : un seul objectif
python -m src.run --graines 10 --poids-qed 1.0
```

Sortie dans `resultats/` : `registre.jsonl` (une molécule distincte par ligne,
avec ses notes, la méthode, la graine et la génération) et `resume.json`
(courbes agrégées, front de Pareto, hypervolume, diversité).

### Vérificateur indépendant

```bash
python -m src.verifier resultats/registre.jsonl
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

**La suite complète passe.** Exécutée le 23 novembre 2026 sur Google Colab,
RDKit 2026.09.2, selfies 2.1.1 : `pytest -q` sort avec le **code 0**, aucun
échec ni erreur.

Relancée depuis un clone vierge **après la traduction complète en français**,
qui a renommé l'intégralité des identifiants : code 0 de nouveau. Le
renommage n'a donc rien casse.

Reproduction depuis un environnement vierge :

```bash
git clone -b c3-evolucao-moleculas <url-du-fork> c3
cd c3/groupe-molevo
pip install -r requirements.txt
pytest -q
```

Les huit tests couvrant la molécule vide (`test_molecula_vazia_e_rejeitada`,
`test_cadeia_degenerada_nao_levanta` et les suivants) ont été écrits **après**
le `TypeError` rencontré en exécution : ils verrouillent un défaut réel, pas un
cas imaginé.

## Résultats et évaluation

Exécution du **23 novembre 2026** : 5 graines, population 60, 60 générations,
3 600 évaluations par méthode et par graine, poids QED 0,5.
RDKit 2026.09.2, selfies 2.1.1. Tous les chiffres ci-dessous viennent de
`notebooks/RESULTATS.md`, régénéré par le notebook.

### 1. La représentation tient sa promesse

| Représentation | Mutations donnant encore une molécule |
|---|---:|
| SMILES | **13,7 %** |
| SELFIES | **100,0 %** |

Facteur **×7,3**, sur 2 000 substitutions d'une unité appliquées à six
médicaments connus. Les 100 % de SELFIES ne sont pas une performance de
l'algorithme : c'est une propriété de la grammaire, obtenue par construction.
Le chiffre de SMILES dépend de l'alphabet de substitution, donné dans la
cellule correspondante du notebook.

### 2. La recherche fonctionne, modestement

| Méthode | Moyenne | Écart-type | Max | Molécules distinctes |
|---|---:|---:|---:|---:|
| **Algorithme génétique** | **0,7163** | 0,0305 | 0,7417 | 5 572 |
| Hill climbing | 0,6640 | 0,0900 | 0,7514 | 4 977 |
| Recherche aléatoire | 0,6262 | 0,0263 | 0,6656 | 12 908 |

Écart GA − aléatoire : **+0,0901**, soit environ 3,4 écarts-types de la ligne
de base. Le gain est réel et n'est pas du bruit.

Deux nuances que la moyenne cache, et qu'il faut dire :

- le **hill climbing a atteint un maximum supérieur** (0,7514 contre 0,7417)
  sur une graine, avec trois fois l'écart-type du GA. Il est instable, pas
  meilleur : le GA gagne en fiabilité, pas en pic ;
- la **recherche aléatoire a trouvé plus du double de molécules distinctes**
  (12 908 contre 5 572) à budget identique. Elle explore davantage ; elle ne
  sait simplement pas où aller.

### 3. Le test le plus sévère : les médicaments réels

| Médicament | Molécules du front qu'il domine | Dominé par |
|---|---:|---:|
| Paracétamol | **13 / 25** | **0** |
| Ibuprofène | 7 | **0** |
| Pénicilline G | 2 | **0** |
| Aspirine | 1 | **0** |
| Caféine | 0 | 9 |
| Atorvastatine | 0 | 23 |

**Quatre des six médicaments ne sont dominés par aucune molécule du front.**
Après 18 000 évaluations, la recherche n'a produit aucune structure qui soit
à la fois plus proche d'un médicament oral *et* plus facile à synthétiser que
le paracétamol — lequel domine à lui seul plus de la moitié du front.

Ce n'est pas un défaut du code : c'est la mesure honnête de ce que vaut le
générateur, et elle situe la difficulté réelle du problème.

### 4. Multiobjectif

- Molécules distinctes évaluées par le GA : **5 572**
- Front de Pareto : **25 molécules**
- Hypervolume, référence (0, 0) : **0,674**

Chiffres issus de l'exécution après correction des doublons entre graines
(voir « Défauts révélés par l'exécution » ci-dessous). L'hypervolume, lui,
n'était pas affecté : une molécule répétée n'ajoute aucune aire dominée.

Les figures, le détail complet et les données brutes sont versionnés dans
`notebooks/` — `RESULTATS.md`, `mesures.json` et les sept PNG. Le notebook
est commité **avec ses sorties** : les figures et les tableaux sont visibles
directement sur GitHub, sans rien exécuter.

### 5. Où l'algorithme est réellement allé

| Critère | Molécule championne | QED | SA |
|---|---|---:|---:|
| Aptitude combinée | `CCCCCCCCO` — octanol | 0,547 | 1,57 |
| QED seul | `NC=[P-]CC(=[P-])COBr` | 0,739 | 6,53 |

La loi de Goodhart s'est manifestée, mais pas comme attendu. Le QED maximal
généré (0,739) reste **en dessous** de l'ibuprofène (0,822) : l'optimiseur n'a
pas fait exploser la métrique. Il a fait deux choses plus intéressantes :

- avec les deux objectifs, il converge vers des **alcools gras simples**
  (octanol, heptanol, octylamine). Ce sont des solvants, pas des médicaments.
  Un poids de 0,5 sur le SA récompense tant la facilité de synthèse que
  l'optimum devient « faites quelque chose de trivial » ;
- avec le QED seul, il produit des **structures chimiquement absurdes**
  (phosphore anionique, aromatiques impossibles) avec un SA supérieur à 6 :
  valides selon la grammaire, irréalisables en laboratoire.

### 6. Diversité : pas d'effondrement

| Mesure | Début | Fin |
|---|---:|---:|
| Diversité interne | 0,984 | 0,860 |
| Molécules inédites par génération (moyenne sur 10) | 21,2 | 20,2 |

Contrairement à ce qu'on attendait, l'exploration reste quasi constante de la
première à la dernière génération. Deux hypothèses testables : soit le taux de
mutation de 0,15 reconstitue la variation assez vite, soit 60 générations sont
trop peu pour que la convergence s'installe. Relancer avec 200 générations
tranche en quelques minutes.

Le nombre de squelettes de la population finale est **0** : ce n'est pas un
bug, les championnes sont des chaînes acycliques et une molécule sans cycle
n'a pas de squelette de Bemis–Murcko.

### Défauts révélés par l'exécution, et corrigés

| Défaut | Symptôme mesuré | Correction |
|---|---|---|
| Molécule vide acceptée comme valide | `TypeError` à la troisième génération : `MolFromSmiles("")` renvoie une molécule sans atome, pas `None` | test du nombre d'atomes dans `para_mol` et dans le calcul des notes |
| Doublons entre graines | front annoncé à 38 entrées pour 25 molécules ; vérificateur signalant 785 « doublons » sur 6 357, soit exactement 6 357 − 5 572 | déduplication par SMILES canonique avant l'analyse — après correction : front de 25, et **0 problème** signalé sur 5 572 molécules |

Le second est instructif : **le vérificateur avait raison**. C'est le notebook
qui lui transmettait la même molécule trouvée par des graines différentes.
C'est exactement le rôle d'un vérificateur indépendant que de faire remonter
ce genre d'écart.

### Protocole

Budget compté en évaluations de la fonction d'aptitude et partagé par les
trois méthodes via la classe `Budget`. Graines 0 à 4, fixées et
consignées. Moyenne et écart-type sur les cinq graines. Comparaison
systématique à deux lignes de base. Les figures et `mesures.json` sont
régénérés ensemble, donc ils ne peuvent pas diverger.

### Ce qui reste à mesurer

- les benchmarks GuacaMol : non abordés ;
- l'effet du nombre de générations sur la diversité : l'effondrement attendu
  n'apparaît pas à 60 générations, une exécution à 200 trancherait ;
- l'effet du poids `peso_qed` : à 0,5 l'optimum est un alcool gras, et une
  balayage de 0 à 1 montrerait à partir de quel poids le résultat cesse
  d'être trivial.

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
- **Plan B :** objectif unique (QED), avec `--poids-qed 1.0`.

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
