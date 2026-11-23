"""Verificador independente das solucoes produzidas.

Entregavel obrigatorio do curso. A palavra que importa e *independente*:
este modulo le o registro gravado pela busca e confere tudo por um
caminho de codigo proprio, sem reaproveitar o cache nem as funcoes de
aptidao usadas durante a execucao.

O erro que ele existe para pegar e banal e devastador: reportar uma
molecula com uma nota que nao e a dela. Acontece por um indice trocado ao
reordenar a populacao, por um cache desatualizado, ou por um registro
gravado antes da ultima mutacao. O grafico fica bonito e o resultado e
falso.

Uso:
    python -m src.verifier registro.jsonl
"""

from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass

from rdkit import Chem, RDLogger
from rdkit.Chem import QED, RDConfig

RDLogger.DisableLog("rdApp.*")

TOLERANCIA = 1e-6


@dataclass(frozen=True)
class Problema:
    tipo: str
    detalhe: str
    linha: int

    def __str__(self) -> str:
        return f"[linha {self.linha}] {self.tipo}: {self.detalhe}"


def _sascorer_proprio():
    """Carrega o sascorer sem passar pelo modulo objectives.

    De proposito: se objectives tiver um bug na normalizacao ou no cache,
    um verificador que o importasse herdaria o mesmo bug e nao veria nada.
    """
    caminho = os.path.join(RDConfig.RDContribDir, "SA_Score")
    if caminho not in sys.path:
        sys.path.append(caminho)
    import sascorer  # noqa: PLC0415

    return sascorer


def verificar(registro: list[dict]) -> list[Problema]:
    """Confere um registro de moleculas. Lista vazia significa tudo certo."""
    sascorer = _sascorer_proprio()
    problemas: list[Problema] = []
    vistos: dict[str, int] = {}

    for n, linha in enumerate(registro, start=1):
        smiles = linha.get("smiles")
        if not smiles:
            problemas.append(Problema("sem_smiles", repr(linha)[:80], n))
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            problemas.append(Problema("smiles_ilegivel", smiles, n))
            continue

        canonico = Chem.MolToSmiles(mol)
        if canonico != smiles:
            problemas.append(
                Problema("smiles_nao_canonico", f"{smiles} -> {canonico}", n)
            )
        if canonico in vistos:
            problemas.append(
                Problema("duplicado", f"{canonico} (ja na linha {vistos[canonico]})", n)
            )
        else:
            vistos[canonico] = n

        # Ida e volta pela representacao, se a cadeia SELFIES foi gravada.
        cadeia = linha.get("selfies")
        if cadeia:
            try:
                import selfies as sf

                volta = Chem.MolFromSmiles(sf.decoder(cadeia))
                if volta is None or Chem.MolToSmiles(volta) != canonico:
                    problemas.append(Problema("selfies_divergente", canonico, n))
            except Exception as erro:  # noqa: BLE001
                problemas.append(Problema("selfies_ilegivel", f"{cadeia}: {erro}", n))

        # Recalculo das notas por caminho proprio.
        if "qed" in linha:
            qed = float(QED.qed(mol))
            if abs(qed - float(linha["qed"])) > TOLERANCIA:
                problemas.append(
                    Problema("qed_divergente", f"{canonico}: {linha['qed']} != {qed:.6f}", n)
                )
        if "sa" in linha:
            sa = float(sascorer.calculateScore(mol))
            if abs(sa - float(linha["sa"])) > TOLERANCIA:
                problemas.append(
                    Problema("sa_divergente", f"{canonico}: {linha['sa']} != {sa:.6f}", n)
                )
        if "sa" in linha and "sa_norm" in linha:
            esperado = max(0.0, min(1.0, (10.0 - float(linha["sa"])) / 9.0))
            if abs(esperado - float(linha["sa_norm"])) > TOLERANCIA:
                problemas.append(
                    Problema("sa_norm_divergente", f"{canonico}: {linha['sa_norm']} != {esperado:.6f}", n)
                )

    return problemas


def verificar_frente(frente_declarada: list[dict]) -> list[Problema]:
    """Confere que a frente de Pareto declarada e mesmo nao dominada."""
    from .pareto import domina

    pontos = [(float(d["qed"]), float(d["sa_norm"])) for d in frente_declarada]
    problemas: list[Problema] = []
    for i, p in enumerate(pontos):
        for j, q in enumerate(pontos):
            if i != j and domina(q, p):
                problemas.append(
                    Problema(
                        "frente_dominada",
                        f"{frente_declarada[i]['smiles']} e dominada por "
                        f"{frente_declarada[j]['smiles']}",
                        i + 1,
                    )
                )
                break
    return problemas


def carregar_jsonl(caminho: str) -> list[dict]:
    with open(caminho, encoding="utf-8") as f:
        return [json.loads(linha) for linha in f if linha.strip()]


def main(argv: list[str] | None = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    if not argv:
        print("uso: python -m src.verifier <registro.jsonl>")
        return 2

    registro = carregar_jsonl(argv[0])
    problemas = verificar(registro)

    print(f"Verificador independente -- {len(registro)} moleculas conferidas")
    if not problemas:
        print("OK: 0 problemas.")
        return 0

    por_tipo: dict[str, int] = {}
    for p in problemas:
        por_tipo[p.tipo] = por_tipo.get(p.tipo, 0) + 1
    print(f"FALHOU: {len(problemas)} problemas.")
    for tipo, quantos in sorted(por_tipo.items()):
        print(f"  {tipo}: {quantos}")
    for p in problemas[:20]:
        print(f"  {p}")
    if len(problemas) > 20:
        print(f"  ... e mais {len(problemas) - 20}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
