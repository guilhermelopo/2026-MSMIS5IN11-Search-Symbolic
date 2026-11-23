# Notebook d'analyse

`analise.ipynb` -- 35 cellules, autonome.

Il installe ses dependances, execute les experiences et produit les six
figures du rapport. **Aucun fichier de resultats prealable n'est necessaire** :
ouvrir, « Executer tout », et les mesures comme les figures sortent ensemble.

## Comment l'executer

**Google Colab** (recommande si RDKit ne se charge pas en local) : televerser
le dossier `groupe-molevo`, ouvrir `notebooks/analise.ipynb`, « Executer tout ».

**En local**, depuis `groupe-molevo/` :

```bash
pip install -r requirements.txt
jupyter lab notebooks/analise.ipynb
```

Duree : environ 3 a 6 minutes avec les parametres par defaut. Tous les
parametres sont dans une seule cellule (« Parametres ») : un seul endroit a
modifier pour refaire toutes les mesures.

## Ce qu'il produit

| Figure | Contenu | Fichier |
|---|---|---|
| 1 | validite sous mutation : SMILES contre SELFIES | `fig1_validite.png` |
| 2 | QED et SA de six medicaments reels | `fig2_medicaments.png` |
| 3 | convergence GA / aleatoire / hill climbing, bande d'ecart-type | `fig3_convergence.png` |
| 4 | front de Pareto et hypervolume | `fig4_pareto.png` |
| 5 | effondrement de la diversite | `fig5_diversite.png` |
| 6 | grilles des meilleures molecules (aptitude, puis QED seul) | affichee |

Plus `mesures.json`, pour que le README et les slides citent exactement les
memes nombres que le notebook.

## Deux choix a expliquer

**Il est autonome et non dependant de `src/run.py`.** Un notebook qui exige
qu'on lance d'abord un script separe se casse des qu'on le partage. Celui-ci
tourne dans Colab sans rien installer localement -- ce qui est aussi le plan B
pour la demonstration si une machine fait defaut le jour de la soutenance.

**Il verifie son propre verificateur.** La derniere cellule falsifie
volontairement une note et confirme que `verifier.py` la detecte. Un
verificateur qui ne trouve jamais rien n'est pas une preuve de qualite : il
peut simplement etre casse.

## A completer apres execution

La section « Conclusions » a la fin du notebook contient cinq phrases a
trous. Les remplir avec les nombres obtenus, puis reporter les memes chiffres
dans la section « Resultats et evaluation » du README du groupe.
