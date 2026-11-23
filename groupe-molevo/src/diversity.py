"""Medidas de diversidade.

O algoritmo genetico converge, e isso nao e so bom: depois de algumas
dezenas de geracoes a populacao vira copias do mesmo vencedor. A curva de
aptidao fica otima e o resultado deixa de ser uma busca.

Para descoberta de farmacos isso e um fracasso disfarcado: o objetivo
nunca foi uma molecula, e sim um conjunto de candidatos diferentes entre
si, porque a maioria vai falhar depois por motivos que nenhuma destas
formulas preve.

Reportar a curva de aptidao junto com a de diversidade e uma das figuras
mais fortes que este projeto pode produzir.
"""

from __future__ import annotations

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold


def impressao(smiles: str, raio: int = 2, bits: int = 2048):
    """Impressao digital de Morgan (ECFP4 para raio=2).

    Um vetor de bits em que cada bit indica a presenca de uma subestrutura
    centrada num atomo, ate `raio` ligacoes de distancia.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=raio, nBits=bits)


def diversidade_interna(lista_smiles: list[str]) -> float:
    """1 menos a semelhanca media de Tanimoto entre todos os pares.

    1.0 = tudo diferente, 0.0 = clones. Custo quadratico no numero de
    moleculas; para populacoes grandes, amostre.
    """
    fps = [f for f in (impressao(s) for s in lista_smiles) if f is not None]
    if len(fps) < 2:
        return 0.0
    sims: list[float] = []
    for i in range(len(fps)):
        sims.extend(DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:]))
    return 1.0 - sum(sims) / len(sims)


def esqueleto(smiles: str) -> str | None:
    """Esqueleto de Bemis-Murcko: a molecula sem as cadeias laterais."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        nucleo = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(nucleo)
    except Exception:  # noqa: BLE001 - moleculas sem anel podem falhar
        return None


def contar_esqueletos(lista_smiles: list[str]) -> int:
    """Quantos esqueletos distintos existem no conjunto.

    Medida mais interpretavel para um quimico do que o Tanimoto medio:
    cem moleculas com tres esqueletos e pouca diversidade real, por mais
    que a semelhanca media nao pareca tao alta.
    """
    return len({e for e in (esqueleto(s) for s in lista_smiles) if e})


def resumo(lista_smiles: list[str]) -> dict:
    """As tres medidas de uma vez, pronto para gravar no registro."""
    distintas = sorted(set(lista_smiles))
    return {
        "n_moleculas": len(lista_smiles),
        "n_distintas": len(distintas),
        "diversidade_interna": diversidade_interna(distintas),
        "n_esqueletos": contar_esqueletos(distintas),
    }
