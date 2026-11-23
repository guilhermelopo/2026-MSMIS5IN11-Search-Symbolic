"""Representacao de moleculas como cadeias SELFIES.

Este modulo isola tudo o que depende do formato textual das moleculas.
O resto do projeto nunca manipula SMILES ou SELFIES diretamente: pede
aqui. Assim, trocar a representacao (por exemplo, para um algoritmo
genetico sobre grafos) mexe num arquivo so.

A propriedade central do SELFIES e que *toda* sequencia de fichas
decodifica para uma molecula valida. O teste
tests/test_representation.py::test_robustez verifica isso.
"""

from __future__ import annotations

import random
from functools import lru_cache

import selfies as sf
from rdkit import Chem, RDLogger

# O RDKit emite avisos em massa ao ler moleculas estranhas, e o algoritmo
# genetico produz muitas. Os avisos nao sao erros; sem isto, a saida de uma
# execucao passa de dezenas de milhares de linhas.
RDLogger.DisableLog("rdApp.*")


# --------------------------------------------------------------------------
# Alfabeto
# --------------------------------------------------------------------------

_ALFABETO: tuple[str, ...] | None = None


def alfabeto() -> tuple[str, ...]:
    """Fichas que a biblioteca garante seguras, em ordem fixa.

    A ordem fixa importa para a reprodutibilidade: `rng.choice` sobre uma
    tupla ordenada da sempre o mesmo resultado para a mesma semente.
    """
    global _ALFABETO
    if _ALFABETO is None:
        _ALFABETO = tuple(sorted(sf.get_semantic_robust_alphabet()))
    return _ALFABETO


def fichas(cadeia: str) -> list[str]:
    """Quebra uma cadeia SELFIES na sua lista de fichas."""
    return list(sf.split_selfies(cadeia))


# --------------------------------------------------------------------------
# Conversoes
# --------------------------------------------------------------------------

def para_smiles(cadeia: str) -> str:
    """SELFIES -> SMILES."""
    return sf.decoder(cadeia)


def de_smiles(smiles: str) -> str:
    """SMILES -> SELFIES."""
    return sf.encoder(smiles)


@lru_cache(maxsize=200_000)
def decodifica_ok(cadeia: str) -> bool:
    """A cadeia decodifica sem erro? (molecula vazia conta como sucesso.)

    E esta a propriedade que o SELFIES garante: a decodificacao nunca
    falha. Use `valida` quando precisar de uma molecula com atomos.
    """
    try:
        return Chem.MolFromSmiles(sf.decoder(cadeia)) is not None
    except Exception:  # noqa: BLE001 - ficha fora do alfabeto robusto
        return False


@lru_cache(maxsize=200_000)
def para_mol(cadeia: str):
    """SELFIES -> objeto Mol do RDKit com ao menos um atomo, ou None.

    ARMADILHA: `Chem.MolFromSmiles("")` NAO devolve None -- devolve uma
    molecula VAZIA, sem atomos. Uma cadeia SELFIES degenerada (uma ficha
    de anel ou de indice sozinha, por exemplo `[Ring1]`) decodifica para
    um SMILES vazio. Sem o teste de numero de atomos abaixo, essa
    molecula vazia atravessa o `is None` e chega ao QED e ao SA score,
    que devolvem valores indefinidos ou None -- era a origem do
    `TypeError: float() argument must be ... not 'NoneType'`.
    """
    try:
        mol = Chem.MolFromSmiles(sf.decoder(cadeia))
    except Exception:  # noqa: BLE001
        return None
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    return mol


def canonico(cadeia: str) -> str | None:
    """SMILES canonico, a chave certa para comparar moleculas.

    Duas cadeias SELFIES diferentes podem descrever a mesma molecula.
    Comparar as cadeias brutas conta duplicatas como se fossem distintas.
    """
    mol = para_mol(cadeia)
    if mol is None:
        return None
    can = Chem.MolToSmiles(mol)
    return can or None          # nunca devolve a cadeia vazia


def valida(cadeia: str) -> bool:
    """A cadeia descreve uma molecula com ao menos um atomo."""
    return para_mol(cadeia) is not None


# --------------------------------------------------------------------------
# Geracao
# --------------------------------------------------------------------------

def aleatoria(rng: random.Random, n_min: int = 6, n_max: int = 25) -> str:
    """Uma cadeia SELFIES aleatoria de n_min a n_max fichas."""
    alfa = alfabeto()
    n = rng.randint(n_min, n_max)
    return "".join(rng.choice(alfa) for _ in range(n))


def populacao_inicial(
    rng: random.Random,
    n: int,
    sementes: list[str] | None = None,
) -> list[str]:
    """Monta a populacao de partida.

    `sementes` e uma lista opcional de SMILES conhecidos. Partir de
    moleculas reais costuma convergir mais rapido do que partir do nada,
    mas enviesa o resultado para a vizinhanca delas -- compare os dois
    no relatorio, e um experimento barato.
    """
    pop: list[str] = []
    if sementes:
        for smi in sementes:
            try:
                pop.append(de_smiles(smi))
            except Exception:  # noqa: BLE001 - SMILES de entrada invalido
                continue
    while len(pop) < n:
        pop.append(aleatoria(rng))
    return pop[:n]
