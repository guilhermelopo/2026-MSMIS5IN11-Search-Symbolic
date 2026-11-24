# Slides de la soutenance

Soutenance : mercredi 25 novembre 2026, 13:15-17:00, salle A202 C. Evaluation
par les pairs et par l'enseignant, sur quatre criteres notes de 0 a 10 :
presentation, theorie, technique et organisation.

Le PDF exporte depuis la presentation est depose a cote de ce fichier.

## Registre visuel

Cahier de laboratoire : papier vieilli, filets d'un pixel au lieu de cadres,
releves en monospace, chiffres en serif. Deux familles seulement, chacune avec
un role : **Newsreader** pour les titres, les chiffres et la prose,
**IBM Plex Mono** pour tout ce qui est mesure, code ou etiquette.

Deux couleurs portent le sens d'un bout a l'autre : **petrole** `#1D5260` pour
ce qui tient (SELFIES, l'algorithme genetique, le verificateur), **vermillon**
`#BE4228` pour ce qui casse (SMILES, les mutations, les limites).

Chaque diapositive porte deux rails : en haut le numero et la section, en bas
le releve de mesure correspondant.

Le logo de l'ecole est en haut a gauche des treize diapositives
(`source/assets/epf.png`, charge dans la presentation comme ressource). Sur
les diapositives a fond encre, le violet du logo tomberait a un contraste de
1,5:1 : il y est donc pose sur une plaque creme. Le logo garde ses couleurs
d'origine dans les deux cas, il n'est jamais recolorie.

## Source

`source/` contient la presentation telle qu'elle est projetee : une
diapositive par fichier HTML, toutes les regles communes dans un seul CSS.

```
source/deck.json        index : titre, ordre des diapositives, polices
source/deck.css         palette, echelle typographique, rails, mouvement
source/slides/*.html    une diapositive par fichier, dans l'ordre ci-dessous
source/assets/epf.png   le logo de l'ecole
```

Derriere le contenu de chaque feuille passe une couche organique : de larges
taches d'encre en courbes fermees, a 5 % d'opacite sur le papier et 16 a 22 %
sur l'encre, plus une onde derriere le titre de la diapositive 5 et une courbe
de tendance lissee par-dessus l'histogramme de la diapositive 8.

Chaque fichier est un seul `<section>` sur une toile fixe de 1920x1080, avec
les notes de l'orateur en dernier dans un `<aside>`.

### Ce qui bouge

- **Trois animations vivantes** (`<x-embed>`, canvas, une page isolee par
  diapositive) : un champ de graphes moleculaires qui derivent et mutent sur
  la couverture ; une machine a mutations SELFIES qui fait defiler cinq
  mutations **reelles** de l'aspirine sur la diapositive 3 ; la population de
  60 individus qui monte generation apres generation sur la diapositive 5,
  elites en or et mutants en vermillon.
- **Des etapes au clic** (`fragment`) sur les diapositives 7, 9 et 10 : les
  courbes se revelent une par une, puis le front de Pareto, puis les
  medicaments reels par-dessus. L'orateur mene le rythme.
- **Des traces SVG** animes au SMIL : les barres de diversite, la chaine de
  l'octanol, les moyennes en pointilles.

Les graphiques sont ecrits a la main a partir de `../notebooks/mesures.json` :
les coordonnees sont calculees depuis le JSON produit par le notebook, aucun
chiffre n'est ressaisi.

Les cinq mutations SELFIES de la diapositive 3 ont ete generees avec
`selfies 2.1.1` a partir de l'aspirine, et chacune decode bien vers la
molecule affichee. La cassure SMILES montree a cote vient de la mesure.

## Les 15 diapositives

| # | Titre | Ce qu'elle montre |
|---|---|---|
| 1 | Conception evolutionnaire de molecules candidates | couverture |
| 2 | Pourquoi une recherche evolutionnaire | ni enumeration, ni gradient, deux objectifs opposes |
| 3 | Une grammaire qui ne casse pas | la meme mutation sur les deux formats : **13,7 % contre 100 %** |
| 4 | Deux notes qui se contredisent | QED et SA sur deux regles graduees, avec les six medicaments |
| 5 | La boucle, en cinq etapes | les cinq etapes et la population qui monte |
| 6 | L'operateur de mutation, en entier | `muter()` en entier : la preuve en code de ce qu'annonce la 3 |
| 7 | Ce qu'il faut fixer pour pouvoir conclure | budget, graines, lignes de base, verificateur |
| 8 | A budget egal, sur cinq graines | **0,7163 contre 0,6262**, ecart de 3,4 ecarts-types |
| 9 | L'effondrement n'a pas eu lieu | la diversite tient, et deux hypotheses testables |
| 10 | Le front de Pareto | 25 molecules, hypervolume 0,674, medicaments sur le meme plan |
| 11 | **Face aux medicaments reels** | quatre des six ne sont domines par rien |
| 12 | **Ce que l'optimiseur a vraiment trouve** | l'octanol gagne ; le QED seul produit du phosphore anionique |
| 13 | Deux bugs | la molecule vide de RDKit, et le verificateur qui avait raison |
| 14 | Comment nous avons travaille | repartition, historique Git, documentation, risques |
| 15 | Reproduire | cinq commandes, plan B versionne, repartition |

La diapositive 6 montre `src/operators.py` tel quel, seule la docstring etant
repliee dans la marge. Elle se lit pour ce qui **n'y est pas** : aucun
`try/except`, aucune verification de validite, aucune boucle de reessai. C'est
la preuve en code de ce que la diapositive 3 annonce en mesure.

La diapositive 14 repond au critere Organisation point par point : repartition
du travail, historique Git, documentation et gestion des risques. Elle porte
aussi la declaration d'usage d'un agent de code, exigee par `AGENTS.md`.

Les travaux cites figurent dans les rails : Krenn et al. 2020 pour SELFIES
(diapositive 3), Bickerton et al. 2012 pour le QED et Ertl et Schuffenhauer
2009 pour le SA score (diapositive 4). La liste complete est dans le README
du dossier.

Les diapositives 10 et 11 valent souvent plus que la 7 : montrer ou la methode
echoue demontre plus de maitrise que n'importe quel bon chiffre.

## Repartition a l'oral

- **Guilherme** : probleme, representation, objectifs, limites, bugs.
- **Lucas** : algorithme, operateurs, protocole, resultats, diversite, Pareto.

Chaque membre doit pouvoir expliquer l'architecture, les choix, les mesures et
les limites de **tout** le projet. « C'est l'agent qui l'a ecrit » n'est pas
une reponse recevable.
