# Donnees

## Provenance

Ce projet **ne depend d'aucun jeu de donnees telecharge**. Les molecules sont
generees par l'algorithme lui-meme, a partir de l'alphabet robuste de la
bibliotheque `selfies`.

Les deux fonctions objectif apportent leurs propres donnees embarquees :

| Source | Contenu | Licence |
|---|---|---|
| `rdkit.Chem.QED` | parametres des fonctions de desirabilite, ajustes par Bickerton et al. (2012) sur des medicaments oraux approuves | BSD-3-Clause, avec RDKit |
| `rdkit Contrib/SA_Score` | table de frequence de fragments derivee de PubChem, d'Ertl et Schuffenhauer (2009) | BSD-3-Clause, avec RDKit |

Aucune donnee personnelle et aucune donnee de sante reelle n'est utilisee a
aucune etape.

## Si vous utilisez des molecules de depart

La fonction `representation.population_initiale` accepte une liste de SMILES
pour demarrer la recherche depuis des molecules reelles plutot que depuis des
chaines aleatoires.

Si vous employez un jeu public pour cela (ChEMBL, ZINC, le jeu de reference de
GuacaMol) :

1. ecrivez ici la source exacte, la version et la licence ;
2. ajoutez un script de telechargement dans `src/` ;
3. **ne commitez pas les fichiers telecharges** -- le `.gitignore` couvre deja
   `data/` ;
4. verifiez que la licence autorise l'usage ; plusieurs interdisent la
   redistribution.

## Comparer les deux points de depart est une bonne experience

Partir de molecules reelles converge souvent plus vite, mais biaise le
resultat vers leur voisinage. Executer les deux configurations et comparer la
diversite finale coute peu et fournit une section entiere du rapport.
