# Slides de la soutenance

Soutenance : mercredi 25 novembre 2026, 13:15-17:00, salle A202 C. Evaluation
par les pairs et par l'enseignant, sur quatre criteres notes de 0 a 10 :
presentation, theorie, technique et organisation.

**Le PDF est `MolEvo-soutenance-C3.pdf`**, a cote de ce fichier.

> **A reexporter.** Ce PDF a seize pages et date de la version precedente du
> jeu. La presentation en compte dix-sept depuis l'ajout de la diapositive 3
> (« Une molecule, ecrite en une ligne »). Les seize pages restent justes une
> a une, mais la troisieme manque et la numerotation imprimee est decalee a
> partir de la quatrieme.

Le PDF sert de livrable et de plan B si le reseau manque en salle. Il ne
remplace pas la presentation projetee : l'export fige les trois animations et
les revelations au clic, donc les courbes de la diapositive 10, le front de la
12 et les medicaments de la 13 y apparaissent d'un coup au lieu d'etre
devoiles un par un.

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

Le logo de l'ecole est en haut a gauche de chaque diapositive
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
sur l'encre, plus une onde derriere le titre de la diapositive 6 et une courbe
de tendance lissee par-dessus l'histogramme de la diapositive 9.

Chaque fichier est un seul `<section>` sur une toile fixe de 1920x1080, avec
les notes de l'orateur en dernier dans un `<aside>`.

### Ce qui bouge

- **Trois animations vivantes** (`<x-embed>`, canvas, une page isolee par
  diapositive) : un champ de graphes moleculaires qui derivent et mutent sur
  la couverture ; une machine a mutations SELFIES qui fait defiler cinq
  mutations **reelles** de l'aspirine sur la diapositive 5 ; la population de
  60 individus qui monte generation apres generation sur la diapositive 7,
  elites en or et mutants en vermillon.
- **Des etapes au clic** (`fragment`) sur les diapositives 10, 12 et 13 : les
  courbes se revelent une par une, puis le front de Pareto, puis les
  medicaments reels par-dessus. L'orateur mene le rythme.
- **Des traces SVG** animes au SMIL : la chaine carbonee du plan, les barres
  de diversite, la chaine de l'octanol, les moyennes en pointilles.

Les graphiques sont ecrits a la main a partir de `../notebooks/mesures.json` :
les coordonnees sont calculees depuis le JSON produit par le notebook, aucun
chiffre n'est ressaisi.

Les cinq mutations SELFIES de la diapositive 5 ont ete generees avec
`selfies 2.1.1` a partir de l'aspirine, et chacune decode bien vers la
molecule affichee. La cassure SMILES montree a cote vient de la mesure.

## Les 17 diapositives

| # | Titre | Ce qu'elle montre |
|---|---|---|
| 1 | Conception evolutionnaire de molecules candidates | couverture |
| 2 | Ce que nous allons montrer | le plan dessine comme une chaine carbonee a cinq atomes |
| 3 | **Une molecule, ecrite en une ligne** | les regles d'ecriture des deux formats, cote a cote |
| 4 | Pourquoi une recherche evolutionnaire | ni enumeration, ni gradient, deux objectifs opposes |
| 5 | Une grammaire qui ne casse pas | la meme mutation sur les deux formats : **13,7 % contre 100 %** |
| 6 | Deux notes qui se contredisent | QED et SA sur deux regles graduees, avec les six medicaments |
| 7 | La boucle, en cinq etapes | les cinq etapes et la population qui monte |
| 8 | L'operateur de mutation, en entier | `muter()` en entier : la preuve en code de ce qu'annonce la 5 |
| 9 | Ce qu'il faut fixer pour pouvoir conclure | budget, graines, lignes de base, verificateur |
| 10 | A budget egal, sur cinq graines | **0,7163 contre 0,6262**, ecart de 3,4 ecarts-types |
| 11 | L'effondrement n'a pas eu lieu | la diversite tient, et deux hypotheses testables |
| 12 | Le front de Pareto | 25 molecules, hypervolume 0,674, medicaments sur le meme plan |
| 13 | **Face aux medicaments reels** | quatre des six ne sont domines par rien |
| 14 | **Ce que l'optimiseur a vraiment trouve** | l'octanol gagne ; le QED seul produit du phosphore anionique |
| 15 | Deux bugs | la molecule vide de RDKit, et le verificateur qui avait raison |
| 16 | Comment nous avons travaille | repartition, historique Git, documentation, risques |
| 17 | Reproduire | cinq commandes, plan B versionne, repartition |

La diapositive 2 donne le plan avant tout le reste : cinq etapes posees comme
les atomes d'une chaine carbonee, les liaisons se dessinant l'une apres
l'autre et un pulse parcourant la chaine. La cinquieme, la critique, est en
vermillon parce qu'elle compte autant que les quatre autres.

La diapositive 3 pose le vocabulaire en une minute, parce que tout le reste en
depend. En haut SMILES : atomes en lettres, liaisons en signes, branches entre
parentheses, cycles ouverts et refermes par le meme chiffre, et l'aspirine
decoupee en groupes fonctionnels. En bas SELFIES : les quatre familles de
symboles, le point autoreferent — `[Branch1]` et `[Ring1]` ne portent pas leur
taille, ils la lisent dans les symboles suivants — et les trois regles qui
rendent le decodage total : liaison abaissee, cycle raccourci, branche
tronquee. C'est ce qui rend vraie la mesure de la diapositive 5.

La diapositive 8 montre `src/operators.py` tel quel, seule la docstring etant
repliee dans la marge. Elle se lit pour ce qui **n'y est pas** : aucun
`try/except`, aucune verification de validite, aucune boucle de reessai. C'est
la preuve en code de ce que les diapositives 3 et 5 annoncent, l'une en
grammaire et l'autre en mesure.

La diapositive 16 repond au critere Organisation point par point : repartition
du travail, historique Git, documentation et gestion des risques. Elle porte
aussi la declaration d'usage d'un agent de code, exigee par `AGENTS.md`.

Les travaux cites figurent dans les rails : Weininger 1988 pour SMILES et
Krenn et al. 2020 pour SELFIES (diapositives 3 et 5), Bickerton et al. 2012
pour le QED et Ertl et Schuffenhauer 2009 pour le SA score (diapositive 6).
La liste complete est dans le README du dossier.

Les diapositives 13 et 14 valent souvent plus que la 8 : montrer ou la methode
echoue demontre plus de maitrise que n'importe quel bon chiffre.

## Repartition a l'oral

- **Guilherme** : probleme, ecritures, representation, objectifs, limites, bugs.
- **Lucas** : algorithme, operateurs, protocole, resultats, diversite, Pareto.

Chaque membre doit pouvoir expliquer l'architecture, les choix, les mesures et
les limites de **tout** le projet. « C'est l'agent qui l'a ecrit » n'est pas
une reponse recevable.
