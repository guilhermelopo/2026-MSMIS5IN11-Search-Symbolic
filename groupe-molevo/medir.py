"""Produz TODOS os numeros que o material explicativo precisa.

Rode numa maquina onde o RDKit funcione. No Google Colab basta:

    !pip -q install rdkit selfies
    !python medir.py

A saida e um bloco compacto, pronto para colar. Cada numero vem com o
parametro que o produziu, para que qualquer pessoa reproduza.

Tempo: cerca de 2 a 5 minutos com os valores padrao.
"""

from __future__ import annotations

import json
import random
import statistics
import sys
import time

# --------------------------------------------------------------------------
# Parametros das medicoes (mude aqui, nao no meio do codigo)
# --------------------------------------------------------------------------
SEMENTES = 5
POP = 60
GERACOES = 60
PESO_QED = 0.5
N_VALIDADE = 2000          # mutacoes por representacao, no teste de validade
SEMENTE_VALIDADE = 12345

FARMACOS = {
    "aspirina": "CC(=O)Oc1ccccc1C(=O)O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "ibuprofeno": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "cafeina": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "penicilina G": "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O",
    "atorvastatina": "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
}


def cabecalho(titulo: str) -> None:
    print("\n" + "=" * 68)
    print(titulo)
    print("=" * 68)


# --------------------------------------------------------------------------
# 1. Validade sob mutacao: SMILES contra SELFIES
# --------------------------------------------------------------------------

def medir_validade() -> dict:
    """A justificativa do projeto inteiro, em dois numeros.

    Aplica a MESMA perturbacao conceitual nas duas representacoes:
    uma alteracao de uma unidade (um carater no SMILES, uma ficha no
    SELFIES) e conta quantas vezes o resultado ainda e uma molecula.
    """
    from rdkit import Chem, RDLogger
    import selfies as sf

    RDLogger.DisableLog("rdApp.*")
    rng = random.Random(SEMENTE_VALIDADE)
    alfabeto_selfies = sorted(sf.get_semantic_robust_alphabet())
    # Alfabeto SMILES: os caracteres que de fato aparecem em SMILES comuns.
    alfabeto_smiles = list("CNOSFPcnos()[]=#123456+-")

    partidas = list(FARMACOS.values())

    val_smiles = 0
    val_selfies = 0
    exemplos_quebrados: list[tuple[str, str]] = []

    for i in range(N_VALIDADE):
        base_smi = partidas[i % len(partidas)]

        # --- SMILES: troca um carater por outro do alfabeto ---
        pos = rng.randrange(len(base_smi))
        mutado = base_smi[:pos] + rng.choice(alfabeto_smiles) + base_smi[pos + 1:]
        if Chem.MolFromSmiles(mutado) is not None:
            val_smiles += 1
        elif len(exemplos_quebrados) < 5:
            exemplos_quebrados.append((base_smi, mutado))

        # --- SELFIES: troca uma ficha por outra do alfabeto robusto ---
        base_sel = list(sf.split_selfies(sf.encoder(base_smi)))
        pos = rng.randrange(len(base_sel))
        base_sel[pos] = rng.choice(alfabeto_selfies)
        if Chem.MolFromSmiles(sf.decoder("".join(base_sel))) is not None:
            val_selfies += 1

    return {
        "n_mutacoes": N_VALIDADE,
        "validade_smiles_pct": 100.0 * val_smiles / N_VALIDADE,
        "validade_selfies_pct": 100.0 * val_selfies / N_VALIDADE,
        "exemplos_quebrados": exemplos_quebrados,
    }


# --------------------------------------------------------------------------
# 2. QED e SA de farmacos conhecidos
# --------------------------------------------------------------------------

def medir_farmacos() -> dict:
    from src import objectives as obj
    from src import representation as rep

    saida = {}
    for nome, smi in FARMACOS.items():
        n = obj.notas(rep.de_smiles(smi))
        saida[nome] = {
            "qed": round(n.qed, 4),
            "sa": round(n.sa, 4),
            "sa_norm": round(n.sa_norm, 4),
        }
    return saida


# --------------------------------------------------------------------------
# 3. GA contra as linhas de base
# --------------------------------------------------------------------------

def medir_busca() -> dict:
    from src import baselines, diversity, ga
    from src import objectives as obj
    from src import pareto

    orcamento_total = POP * GERACOES
    resultados: dict[str, list] = {"ga": [], "aleatoria": [], "hill": []}
    achados_ga = []

    for semente in range(SEMENTES):
        obj.limpar_cache()
        t0 = time.time()

        r_ga = ga.evoluir(
            semente=semente, tamanho_pop=POP, geracoes=GERACOES, peso_qed=PESO_QED
        )
        r_al = baselines.busca_aleatoria(
            semente, obj.Orcamento(orcamento_total), PESO_QED
        )
        r_hc = baselines.hill_climbing(
            semente, obj.Orcamento(orcamento_total), PESO_QED
        )

        resultados["ga"].append(r_ga)
        resultados["aleatoria"].append(r_al)
        resultados["hill"].append(r_hc)
        achados_ga.extend(r_ga.arquivo.itens.values())

        print(
            f"  semente {semente}: GA {r_ga.melhor_aptidao:.4f} | "
            f"aleatoria {r_al.melhor_aptidao:.4f} | "
            f"hill {r_hc.melhor_aptidao:.4f}  ({time.time() - t0:.1f}s)",
            flush=True,
        )

    resumo = {}
    for nome, lista in resultados.items():
        melhores = [r.melhor_aptidao for r in lista]
        resumo[nome] = {
            "melhor_media": round(statistics.fmean(melhores), 4),
            "melhor_desvio": round(statistics.pstdev(melhores), 4),
            "melhor_max": round(max(melhores), 4),
            "distintas": len({a.smiles for r in lista for a in r.arquivo.itens.values()}),
        }

    # Curva do GA em alguns pontos, para a figura de convergencia.
    n = min(len(r.historico_melhor) for r in resultados["ga"])
    curva_ga = [
        round(statistics.fmean([r.historico_melhor[i] for r in resultados["ga"]]), 4)
        for i in range(n)
    ]

    # Diversidade na populacao final da ultima semente.
    pop_final = resultados["ga"][-1].populacao_final
    from src import representation as rep

    smiles_final = [s for s in (rep.canonico(c) for c in pop_final) if s]
    div = diversity.resumo(smiles_final) if smiles_final else {}

    # Frente de Pareto sobre tudo o que o GA avaliou.
    pontos = [(a.qed, a.sa_norm) for a in achados_ga]
    indices = pareto.indices_nao_dominados(pontos) if pontos else []
    frente = [achados_ga[i] for i in indices]
    hv = pareto.hipervolume_2d([(a.qed, a.sa_norm) for a in frente])

    melhores_qed = sorted(achados_ga, key=lambda a: a.qed, reverse=True)[:5]

    return {
        "orcamento_por_metodo": orcamento_total,
        "sementes": SEMENTES,
        "resumo": resumo,
        "curva_ga_media": curva_ga,
        "diversidade_populacao_final": div,
        "n_frente_pareto": len(frente),
        "hipervolume": round(hv, 4),
        "frente_amostra": [
            {"smiles": a.smiles, "qed": round(a.qed, 3), "sa": round(a.sa, 2)}
            for a in sorted(frente, key=lambda a: a.qed, reverse=True)[:8]
        ],
        "top_qed": [
            {"smiles": a.smiles, "qed": round(a.qed, 3), "sa": round(a.sa, 2)}
            for a in melhores_qed
        ],
        "total_avaliadas_distintas": len({a.smiles for a in achados_ga}),
    }


# --------------------------------------------------------------------------
def main() -> int:
    print("Medicoes do projeto C3")
    print(f"Parametros: {SEMENTES} sementes, pop {POP}, {GERACOES} geracoes, "
          f"peso_qed {PESO_QED}")

    try:
        import rdkit
        import selfies
        print(f"rdkit {rdkit.__version__} | selfies {selfies.__version__}")
    except ImportError as erro:
        print(f"\nFALTA INSTALAR: {erro}")
        print("pip install rdkit selfies")
        return 1

    tudo = {}

    cabecalho("1. VALIDADE SOB MUTACAO  (a justificativa do projeto)")
    v = medir_validade()
    tudo["validade"] = v
    print(f"  {v['n_mutacoes']} mutacoes de uma unidade em cada representacao")
    print(f"  SMILES  : {v['validade_smiles_pct']:.1f}% continuam moleculas validas")
    print(f"  SELFIES : {v['validade_selfies_pct']:.1f}% continuam moleculas validas")
    if v["exemplos_quebrados"]:
        print("  exemplos de SMILES quebrados:")
        for antes, depois in v["exemplos_quebrados"][:3]:
            print(f"    {antes}  ->  {depois}")

    cabecalho("2. QED E SA DE FARMACOS CONHECIDOS")
    f = medir_farmacos()
    tudo["farmacos"] = f
    print(f"  {'molecula':<16}{'QED':>8}{'SA':>8}{'sa_norm':>10}")
    for nome, d in f.items():
        print(f"  {nome:<16}{d['qed']:>8.3f}{d['sa']:>8.2f}{d['sa_norm']:>10.3f}")

    cabecalho("3. GA CONTRA AS LINHAS DE BASE")
    b = medir_busca()
    tudo["busca"] = b
    print(f"\n  orcamento: {b['orcamento_por_metodo']} avaliacoes por metodo "
          f"e por semente")
    print(f"  {'metodo':<12}{'media':>9}{'desvio':>9}{'max':>9}{'distintas':>11}")
    for nome, d in b["resumo"].items():
        print(f"  {nome:<12}{d['melhor_media']:>9.4f}{d['melhor_desvio']:>9.4f}"
              f"{d['melhor_max']:>9.4f}{d['distintas']:>11}")
    print(f"\n  moleculas distintas avaliadas pelo GA: {b['total_avaliadas_distintas']}")
    print(f"  frente de Pareto: {b['n_frente_pareto']} moleculas, "
          f"hipervolume {b['hipervolume']}")
    if b["diversidade_populacao_final"]:
        d = b["diversidade_populacao_final"]
        print(f"  populacao final: {d.get('n_distintas')} distintas de "
              f"{d.get('n_moleculas')}, diversidade interna "
              f"{d.get('diversidade_interna', 0):.3f}, "
              f"{d.get('n_esqueletos')} esqueletos")
    print("\n  melhores por QED (atencao: QED alto costuma ser sinal de alerta):")
    for m in b["top_qed"]:
        print(f"    QED {m['qed']:.3f}  SA {m['sa']:.2f}  {m['smiles'][:52]}")

    with open("medicoes.json", "w", encoding="utf-8") as fh:
        json.dump(tudo, fh, ensure_ascii=False, indent=2)
    cabecalho("PRONTO")
    print("  Gravado em medicoes.json -- mande este arquivo para colar no material.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
