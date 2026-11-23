"""Funcoes objetivo: QED, SA score e a sua combinacao.

Duas notas com direcoes opostas:

  QED  0..1   maior e melhor   "tem perfil de farmaco oral?"
  SA   1..10  MENOR e melhor   "da para fabricar?"

Para combinar, o SA e invertido e normalizado para 0..1 com
`normalizar_sa`. Toda a aritmetica de objetivos mora aqui; nenhum outro
modulo deve inverter escalas por conta propria.

O orcamento de avaliacoes tambem e contado aqui, porque este e o unico
lugar por onde todas as avaliacoes passam. Isso torna impossivel comparar
dois metodos com orcamentos diferentes por descuido.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from functools import lru_cache

from rdkit.Chem import QED, RDConfig

from . import representation as rep

# --------------------------------------------------------------------------
# SA score: mora em Contrib, nao no nucleo do RDKit
# --------------------------------------------------------------------------

_sascorer = None


def _carregar_sascorer():
    """Importa o sascorer, remendando o sys.path uma unica vez.

    O SA score nao faz parte do RDKit instalado por pip; vem na pasta de
    contribuicoes. Este remendo fica isolado aqui de proposito -- espalhado
    pelo codigo, e a primeira coisa que quebra na maquina de outra pessoa.
    """
    global _sascorer
    if _sascorer is None:
        caminho = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if caminho not in sys.path:
            sys.path.append(caminho)
        import sascorer  # noqa: PLC0415 - import tardio e proposital

        _sascorer = sascorer
    return _sascorer


# --------------------------------------------------------------------------
# Orcamento
# --------------------------------------------------------------------------

class Orcamento:
    """Conta avaliacoes de aptidao e para quando o limite acaba.

    E a moeda do projeto: comparar o GA com as linhas de base so e honesto
    com o mesmo numero de avaliacoes. Todos os metodos recebem um Orcamento
    e param quando `esgotado` fica verdadeiro.
    """

    def __init__(self, total: int) -> None:
        self.total = total
        self.usadas = 0

    def gastar(self, n: int = 1) -> None:
        self.usadas += n

    @property
    def restantes(self) -> int:
        return max(0, self.total - self.usadas)

    @property
    def esgotado(self) -> bool:
        return self.usadas >= self.total

    def __repr__(self) -> str:
        return f"Orcamento({self.usadas}/{self.total})"


# --------------------------------------------------------------------------
# Notas
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Notas:
    """As notas de uma molecula. `valida` falso significa cadeia ilegivel."""

    qed: float
    sa: float
    sa_norm: float
    valida: bool = True

    @property
    def objetivos(self) -> tuple[float, float]:
        """Par para maximizar, usado pelo modulo pareto."""
        return (self.qed, self.sa_norm)


NOTAS_INVALIDAS = Notas(qed=0.0, sa=10.0, sa_norm=0.0, valida=False)


def normalizar_sa(sa: float) -> float:
    """SA 1..10 (minimizar) -> 0..1 (maximizar).

    SA 1  (trivial)  -> 1.0
    SA 10 (pesadelo) -> 0.0
    """
    return max(0.0, min(1.0, (10.0 - sa) / 9.0))


@lru_cache(maxsize=200_000)
def _notas_por_smiles(smiles_canonico: str) -> Notas:
    """Nucleo do calculo, com cache pelo SMILES canonico.

    O cache e pelo canonico e nao pela cadeia SELFIES porque varias cadeias
    diferentes descrevem a mesma molecula -- o ganho e grande depois de
    algumas geracoes, quando a populacao comeca a repetir.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles_canonico)
    if mol is None:
        return NOTAS_INVALIDAS
    sa = float(_carregar_sascorer().calculateScore(mol))
    return Notas(qed=float(QED.qed(mol)), sa=sa, sa_norm=normalizar_sa(sa))


def notas(cadeia: str, orcamento: Orcamento | None = None) -> Notas:
    """Notas de uma cadeia SELFIES. Gasta uma unidade do orcamento."""
    if orcamento is not None:
        orcamento.gastar()
    can = rep.canonico(cadeia)
    if can is None:
        return NOTAS_INVALIDAS
    return _notas_por_smiles(can)


# --------------------------------------------------------------------------
# Escalarizacao
# --------------------------------------------------------------------------

def aptidao(cadeia: str, peso_qed: float = 0.5, orcamento: Orcamento | None = None) -> float:
    """Soma ponderada dos dois objetivos: o MVP.

    peso_qed = 1.0  -> so QED (o plano B do catalogo)
    peso_qed = 0.0  -> so facilidade de sintese
    peso_qed = 0.5  -> compromisso neutro

    O valor do peso e arbitrario, e esse e justamente o defeito da soma
    ponderada. Varra-o de 0 a 1 no relatorio e compare com a frente de
    Pareto do modulo pareto.
    """
    n = notas(cadeia, orcamento)
    return peso_qed * n.qed + (1.0 - peso_qed) * n.sa_norm


def limpar_cache() -> None:
    """Zera os caches. Util entre execucoes com sementes diferentes."""
    _notas_por_smiles.cache_clear()
    rep.para_mol.cache_clear()
