"""Verificateur independant des solutions produites.

Livrable obligatoire du cours. Le mot qui compte est *independant* : ce
module lit le registre ecrit par la recherche et revérifie tout par un
chemin de code propre, sans reutiliser le cache ni les fonctions d'aptitude
employees pendant l'execution.

L'erreur qu'il existe pour attraper est banale et devastatrice : rapporter
une molecule avec une note qui n'est pas la sienne. Cela arrive par un
indice echange lors d'un tri, par un cache perime, ou par un enregistrement
ecrit avant la derniere mutation. Le graphique reste beau et le resultat
est faux.

Utilisation :
    python -m src.verifier registre.jsonl
"""

from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass

from rdkit import Chem, RDLogger
from rdkit.Chem import QED, RDConfig

RDLogger.DisableLog("rdApp.*")

TOLERANCE = 1e-6


@dataclass(frozen=True)
class Probleme:
    type: str
    detail: str
    ligne: int

    def __str__(self) -> str:
        return f"[ligne {self.ligne}] {self.type} : {self.detail}"


def _sascorer_propre():
    """Charge sascorer sans passer par le module objectives.

    A dessein : si objectives avait un bug dans la normalisation ou dans le
    cache, un verificateur qui l'importerait heriterait du meme bug et ne
    verrait rien.
    """
    chemin = os.path.join(RDConfig.RDContribDir, "SA_Score")
    if chemin not in sys.path:
        sys.path.append(chemin)
    import sascorer  # noqa: PLC0415

    return sascorer


def verifier_registre(registre: list[dict]) -> list[Probleme]:
    """Verifie un registre de molecules. Liste vide signifie tout bon."""
    sascorer = _sascorer_propre()
    problemes: list[Probleme] = []
    vus: dict[str, int] = {}

    for n, ligne in enumerate(registre, start=1):
        smiles = ligne.get("smiles")
        if not smiles:
            problemes.append(Probleme("sans_smiles", repr(ligne)[:80], n))
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            problemes.append(Probleme("smiles_illisible", smiles, n))
            continue

        canonique = Chem.MolToSmiles(mol)
        if canonique != smiles:
            problemes.append(
                Probleme("smiles_non_canonique", f"{smiles} -> {canonique}", n)
            )
        if canonique in vus:
            problemes.append(
                Probleme("doublon", f"{canonique} (deja ligne {vus[canonique]})", n)
            )
        else:
            vus[canonique] = n

        # Aller-retour par la representation, si la chaine SELFIES est la.
        chaine = ligne.get("selfies")
        if chaine:
            try:
                import selfies as sf

                retour = Chem.MolFromSmiles(sf.decoder(chaine))
                if retour is None or Chem.MolToSmiles(retour) != canonique:
                    problemes.append(Probleme("selfies_divergent", canonique, n))
            except Exception as erreur:  # noqa: BLE001
                problemes.append(
                    Probleme("selfies_illisible", f"{chaine} : {erreur}", n))

        # Recalcul des notes par un chemin propre.
        if "qed" in ligne:
            qed = float(QED.qed(mol))
            if abs(qed - float(ligne["qed"])) > TOLERANCE:
                problemes.append(
                    Probleme("qed_divergent",
                             f"{canonique} : {ligne['qed']} != {qed:.6f}", n))
        if "sa" in ligne:
            sa = float(sascorer.calculateScore(mol))
            if abs(sa - float(ligne["sa"])) > TOLERANCE:
                problemes.append(
                    Probleme("sa_divergent",
                             f"{canonique} : {ligne['sa']} != {sa:.6f}", n))
        if "sa" in ligne and "sa_norm" in ligne:
            attendu = max(0.0, min(1.0, (10.0 - float(ligne["sa"])) / 9.0))
            if abs(attendu - float(ligne["sa_norm"])) > TOLERANCE:
                problemes.append(
                    Probleme("sa_norm_divergent",
                             f"{canonique} : {ligne['sa_norm']} != {attendu:.6f}", n))

    return problemes


def verifier_front(front_declare: list[dict]) -> list[Probleme]:
    """Verifie que le front de Pareto declare est bien non domine."""
    from .pareto import domine

    points = [(float(d["qed"]), float(d["sa_norm"])) for d in front_declare]
    problemes: list[Probleme] = []
    for i, p in enumerate(points):
        for j, q in enumerate(points):
            if i != j and domine(q, p):
                problemes.append(
                    Probleme(
                        "front_domine",
                        f"{front_declare[i]['smiles']} est dominee par "
                        f"{front_declare[j]['smiles']}",
                        i + 1,
                    )
                )
                break
    return problemes


def charger_jsonl(chemin: str) -> list[dict]:
    with open(chemin, encoding="utf-8") as f:
        return [json.loads(ligne) for ligne in f if ligne.strip()]


def main(argv: list[str] | None = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    if not argv:
        print("usage : python -m src.verifier <registre.jsonl>")
        return 2

    registre = charger_jsonl(argv[0])
    problemes = verifier_registre(registre)

    print(f"Verificateur independant -- {len(registre)} molecules verifiees")
    if not problemes:
        print("OK : 0 probleme.")
        return 0

    par_type: dict[str, int] = {}
    for p in problemes:
        par_type[p.type] = par_type.get(p.type, 0) + 1
    print(f"ECHEC : {len(problemes)} problemes.")
    for type_, combien in sorted(par_type.items()):
        print(f"  {type_} : {combien}")
    for p in problemes[:20]:
        print(f"  {p}")
    if len(problemes) > 20:
        print(f"  ... et {len(problemes) - 20} de plus")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
