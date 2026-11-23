"""Linha de comando: roda os tres metodos e grava o registro.

Exemplos:

    # rodada rapida, para ver se tudo funciona
    python -m src.run --sementes 1 --pop 30 --geracoes 10

    # rodada do relatorio: 10 sementes, orcamento de 10 000 avaliacoes
    python -m src.run --sementes 10 --pop 100 --geracoes 100

    # plano B do catalogo: um objetivo so
    python -m src.run --peso-qed 1.0

A saida vai para resultados/: um JSONL com todas as moleculas distintas e
um JSON com as curvas. O notebook le esses arquivos; nunca chama a busca.
Assim, refazer uma figura nao exige rodar o GA de novo.
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


def _curvas(resultados: list[ga.Resultado]) -> dict:
    """Agrega as curvas de varias sementes em media e desvio.

    Uma execucao unica e ruido, nao resultado. O relatorio reporta media e
    desvio sobre pelo menos cinco sementes, de preferencia dez.
    """
    if not resultados:
        return {}
    n = min(len(r.historico_melhor) for r in resultados)
    melhor_por_passo = [
        [r.historico_melhor[i] for r in resultados] for i in range(n)
    ]
    return {
        "avaliacoes": resultados[0].historico_avaliacoes[:n],
        "melhor_media": [statistics.fmean(v) for v in melhor_por_passo],
        "melhor_desvio": [
            statistics.pstdev(v) if len(v) > 1 else 0.0 for v in melhor_por_passo
        ],
        "n_sementes": len(resultados),
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Projeto C3 -- evolucao de moleculas")
    p.add_argument("--sementes", type=int, default=5, help="quantas sementes aleatorias")
    p.add_argument("--pop", type=int, default=100, help="tamanho da populacao")
    p.add_argument("--geracoes", type=int, default=100)
    p.add_argument("--peso-qed", type=float, default=0.5,
                   help="1.0 = so QED (plano B); 0.5 = compromisso")
    p.add_argument("--p-mutacao", type=float, default=0.15)
    p.add_argument("--p-cruzamento", type=float, default=0.8)
    p.add_argument("--k-torneio", type=int, default=3)
    p.add_argument("--elites", type=int, default=2)
    p.add_argument("--saida", default="resultados")
    args = p.parse_args(argv)

    orcamento_total = args.pop * args.geracoes
    destino = Path(args.saida)
    destino.mkdir(parents=True, exist_ok=True)

    print(f"Orcamento por semente e por metodo: {orcamento_total} avaliacoes")
    print(f"Sementes: {args.sementes} | peso_qed = {args.peso_qed}\n")

    execucoes: dict[str, list[ga.Resultado]] = {"ga": [], "aleatoria": [], "hill": []}
    todos: list[ga.Achado] = []

    for semente in range(args.sementes):
        t0 = time.time()
        obj.limpar_cache()

        r_ga = ga.evoluir(
            semente=semente, tamanho_pop=args.pop, geracoes=args.geracoes,
            elites=args.elites, k_torneio=args.k_torneio,
            p_cruzamento=args.p_cruzamento, p_mutacao=args.p_mutacao,
            peso_qed=args.peso_qed,
        )
        r_al = baselines.busca_aleatoria(
            semente, obj.Orcamento(orcamento_total), args.peso_qed
        )
        r_hc = baselines.hill_climbing(
            semente, obj.Orcamento(orcamento_total), args.peso_qed, args.p_mutacao
        )

        for chave, r in (("ga", r_ga), ("aleatoria", r_al), ("hill", r_hc)):
            execucoes[chave].append(r)
            todos.extend(r.arquivo.itens.values())

        print(
            f"  semente {semente}: "
            f"GA {r_ga.melhor_aptidao:.4f} | "
            f"aleatoria {r_al.melhor_aptidao:.4f} | "
            f"hill {r_hc.melhor_aptidao:.4f} "
            f"({time.time() - t0:.1f}s)"
        )

    # Registro: uma molecula distinta por linha, para o verificador e o notebook.
    caminho_registro = destino / "registro.jsonl"
    vistos: set[str] = set()
    with open(caminho_registro, "w", encoding="utf-8") as f:
        for a in todos:
            chave = f"{a.metodo}|{a.semente}|{a.smiles}"
            if chave in vistos:
                continue
            vistos.add(chave)
            f.write(json.dumps(a.como_dict(), ensure_ascii=False) + "\n")

    # Frente de Pareto sobre tudo o que o GA avaliou, calculada a posteriori.
    # Nao exige mudar o algoritmo e ja rende a figura do capitulo 10.
    do_ga = [a for a in todos if a.metodo == "ga"]
    pontos = [(a.qed, a.sa_norm) for a in do_ga]
    indices = pareto.indices_nao_dominados(pontos) if pontos else []
    frente = [do_ga[i].como_dict() for i in indices]

    resumo = {
        "parametros": vars(args),
        "orcamento_por_execucao": orcamento_total,
        "curvas": {k: _curvas(v) for k, v in execucoes.items()},
        "melhor_por_metodo": {
            k: max((r.melhor_aptidao for r in v), default=0.0)
            for k, v in execucoes.items()
        },
        "moleculas_distintas_por_metodo": {
            k: len({a.smiles for a in todos if a.metodo.startswith(k)})
            for k in execucoes
        },
        "frente_pareto": frente,
        "hipervolume": pareto.hipervolume_2d([(d["qed"], d["sa_norm"]) for d in frente]),
    }

    populacao_final = [
        a.smiles for a in execucoes["ga"][-1].arquivo.itens.values()
    ][:200]
    if populacao_final:
        resumo["diversidade_ga"] = diversity.resumo(populacao_final)

    with open(destino / "resumo.json", "w", encoding="utf-8") as f:
        json.dump(resumo, f, ensure_ascii=False, indent=2)

    print(f"\nRegistro: {caminho_registro} ({len(vistos)} linhas)")
    print(f"Resumo:   {destino / 'resumo.json'}")
    print(f"Frente de Pareto: {len(frente)} moleculas | "
          f"hipervolume {resumo['hipervolume']:.4f}")
    print("\nConfira com:  python -m src.verifier " + str(caminho_registro))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
