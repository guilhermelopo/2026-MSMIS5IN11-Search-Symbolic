"""Ligne de commande : execute les trois methodes et ecrit le registre.

Exemples :

    # execution rapide, pour verifier que tout fonctionne
    python -m src.run --graines 1 --pop 30 --generations 10

    # execution du rapport : 10 graines, budget de 10 000 evaluations
    python -m src.run --graines 10 --pop 100 --generations 100

    # plan B du catalogue : un seul objectif
    python -m src.run --poids-qed 1.0

La sortie va dans resultats/ : un JSONL avec toutes les molecules distinctes
et un JSON avec les courbes. Le notebook lit ces fichiers ; il n'appelle
jamais la recherche. Refaire une figure ne coute donc pas une execution.
"""

from __future__ import annotations

import argparse
import json
import statistics
import time
from pathlib import Path

from . import baselines, diversity, ga
from . import objectives as obj
from . import pareto


def _courbes(resultats: list[ga.Resultat]) -> dict:
    """Agrege les courbes de plusieurs graines en moyenne et ecart-type.

    Une execution unique est du bruit, pas un resultat. Le rapport donne
    moyenne et ecart-type sur au moins cinq graines, dix de preference.
    """
    if not resultats:
        return {}
    n = min(len(r.historique_meilleur) for r in resultats)
    meilleur_par_pas = [
        [r.historique_meilleur[i] for r in resultats] for i in range(n)
    ]
    return {
        "evaluations": resultats[0].historique_evaluations[:n],
        "meilleur_moyenne": [statistics.fmean(v) for v in meilleur_par_pas],
        "meilleur_ecart": [
            statistics.pstdev(v) if len(v) > 1 else 0.0 for v in meilleur_par_pas
        ],
        "n_graines": len(resultats),
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Projet C3 -- conception evolutionnaire de molecules")
    p.add_argument("--graines", type=int, default=5,
                   help="nombre de graines aleatoires")
    p.add_argument("--pop", type=int, default=100, help="taille de la population")
    p.add_argument("--generations", type=int, default=100)
    p.add_argument("--poids-qed", type=float, default=0.5, dest="poids_qed",
                   help="poids du QED dans la somme ponderee (1.0 = plan B)")
    p.add_argument("--p-mutation", type=float, default=0.15, dest="p_mutation")
    p.add_argument("--p-croisement", type=float, default=0.8, dest="p_croisement")
    p.add_argument("--k-tournoi", type=int, default=3, dest="k_tournoi")
    p.add_argument("--elites", type=int, default=2)
    p.add_argument("--sortie", default="resultats")
    args = p.parse_args(argv)

    budget_total = args.pop * args.generations
    destination = Path(args.sortie)
    destination.mkdir(parents=True, exist_ok=True)

    print(f"Budget par graine et par methode : {budget_total} evaluations")
    print(f"Graines : {args.graines} | poids_qed = {args.poids_qed}\n")

    executions: dict[str, list[ga.Resultat]] = {
        "ga": [], "aleatoire": [], "hill": []}
    toutes: list[ga.Trouvaille] = []

    for graine in range(args.graines):
        t0 = time.time()
        obj.vider_cache()

        r_ga = ga.evoluer(
            graine=graine, taille_pop=args.pop, generations=args.generations,
            elites=args.elites, k_tournoi=args.k_tournoi,
            p_croisement=args.p_croisement, p_mutation=args.p_mutation,
            poids_qed=args.poids_qed,
        )
        r_al = baselines.recherche_aleatoire(
            graine, obj.Budget(budget_total), args.poids_qed
        )
        r_hc = baselines.hill_climbing(
            graine, obj.Budget(budget_total), args.poids_qed, args.p_mutation
        )

        for cle, r in (("ga", r_ga), ("aleatoire", r_al), ("hill", r_hc)):
            executions[cle].append(r)
            toutes.extend(r.archive.items.values())

        print(
            f"  graine {graine} : "
            f"GA {r_ga.meilleure_aptitude:.4f} | "
            f"aleatoire {r_al.meilleure_aptitude:.4f} | "
            f"hill {r_hc.meilleure_aptitude:.4f} "
            f"({time.time() - t0:.1f} s)"
        )

    # Registre : une molecule distincte par ligne, pour le verificateur et
    # le notebook.
    chemin_registre = destination / "registre.jsonl"
    vus: set[str] = set()
    with open(chemin_registre, "w", encoding="utf-8") as f:
        for t in toutes:
            cle = f"{t.methode}|{t.graine}|{t.smiles}"
            if cle in vus:
                continue
            vus.add(cle)
            f.write(json.dumps(t.en_dict(), ensure_ascii=False) + "\n")

    # Front de Pareto sur tout ce que le GA a evalue, calcule a posteriori.
    # Cela n'exige pas de modifier l'algorithme et donne deja la figure.
    # Deduplication entre graines : l'archive ne deduplique que dans une
    # execution, et la meme molecule trouvee par trois graines compterait
    # trois fois.
    par_smiles: dict[str, ga.Trouvaille] = {}
    for t in toutes:
        if t.methode != "ga":
            continue
        precedent = par_smiles.get(t.smiles)
        if precedent is None or t.aptitude > precedent.aptitude:
            par_smiles[t.smiles] = t
    du_ga = list(par_smiles.values())

    points = [(t.qed, t.sa_norm) for t in du_ga]
    indices = pareto.indices_non_domines(points) if points else []
    front = [du_ga[i].en_dict() for i in indices]

    resume = {
        "parametres": vars(args),
        "budget_par_execution": budget_total,
        "courbes": {k: _courbes(v) for k, v in executions.items()},
        "meilleur_par_methode": {
            k: max((r.meilleure_aptitude for r in v), default=0.0)
            for k, v in executions.items()
        },
        "molecules_distinctes_par_methode": {
            k: len({t.smiles for t in toutes if t.methode.startswith(k)})
            for k in executions
        },
        "front_pareto": front,
        "hypervolume": pareto.hypervolume_2d(
            [(d["qed"], d["sa_norm"]) for d in front]),
    }

    smiles_ga = [t.smiles for t in du_ga][:200]
    if smiles_ga:
        resume["diversite_ga"] = diversity.resume(smiles_ga)

    with open(destination / "resume.json", "w", encoding="utf-8") as f:
        json.dump(resume, f, ensure_ascii=False, indent=2)

    print(f"\nRegistre : {chemin_registre} ({len(vus)} lignes)")
    print(f"Resume   : {destination / 'resume.json'}")
    print(f"Front de Pareto : {len(front)} molecules | "
          f"hypervolume {resume['hypervolume']:.4f}")
    print("\nVerifier avec :  python -m src.verifier " + str(chemin_registre))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
