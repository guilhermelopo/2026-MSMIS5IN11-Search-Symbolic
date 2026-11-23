"""O algoritmo genetico.

O laco principal e curto de proposito; toda a complexidade real esta nos
operadores e nas funcoes objetivo. Duas decisoes de projeto merecem
destaque:

1. O *arquivo de achados* guarda toda molecula distinta avaliada, nao so a
   populacao final. Isso resolve a maior parte do problema de colapso de
   diversidade por tres linhas de codigo: mesmo que a populacao vire
   clones, o arquivo preserva tudo o que foi encontrado no caminho.

2. O orcamento de avaliacoes e compartilhado com as linhas de base, de
   modo que a comparacao e justa por construcao e nao por disciplina.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field

from . import objectives as obj
from . import operators as ops
from . import representation as rep


@dataclass
class Achado:
    """Uma molecula distinta encontrada em algum momento da busca."""

    smiles: str
    selfies: str
    qed: float
    sa: float
    sa_norm: float
    aptidao: float
    geracao: int
    metodo: str
    semente: int

    def como_dict(self) -> dict:
        return self.__dict__.copy()


@dataclass
class Arquivo:
    """Coleciona moleculas distintas, indexadas pelo SMILES canonico.

    Guardar por SMILES canonico e essencial: a mesma molecula aparece sob
    varias cadeias SELFIES diferentes, e indexar pela cadeia bruta inflaria
    a contagem de "moleculas distintas" sem que nada de novo tivesse sido
    descoberto.
    """

    itens: dict[str, Achado] = field(default_factory=dict)

    def registrar(self, achado: Achado) -> bool:
        """Devolve True se a molecula era inedita."""
        if achado.smiles in self.itens:
            return False
        self.itens[achado.smiles] = achado
        return True

    def melhores(self, n: int = 10) -> list[Achado]:
        return sorted(self.itens.values(), key=lambda a: a.aptidao, reverse=True)[:n]

    def __len__(self) -> int:
        return len(self.itens)


@dataclass
class Resultado:
    """O que uma execucao devolve."""

    metodo: str
    semente: int
    arquivo: Arquivo
    historico_melhor: list[float]      # melhor aptidao ate cada geracao
    historico_media: list[float]
    historico_avaliacoes: list[int]    # avaliacoes gastas ate cada geracao
    populacao_final: list[str]

    @property
    def melhor_aptidao(self) -> float:
        return self.historico_melhor[-1] if self.historico_melhor else 0.0


def _avaliar(
    pop: list[str],
    peso_qed: float,
    orcamento: obj.Orcamento,
) -> list[obj.Notas]:
    return [obj.notas(s, orcamento) for s in pop]


def evoluir(
    semente: int,
    tamanho_pop: int = 100,
    geracoes: int = 100,
    elites: int = 2,
    k_torneio: int = 3,
    p_cruzamento: float = 0.8,
    p_mutacao: float = 0.15,
    peso_qed: float = 0.5,
    sementes_smiles: list[str] | None = None,
    orcamento: obj.Orcamento | None = None,
) -> Resultado:
    """Roda o algoritmo genetico e devolve o Resultado.

    Se `orcamento` nao for passado, e criado um de tamanho_pop * geracoes,
    que e exatamente o custo de uma execucao completa.
    """
    rng = random.Random(semente)
    if orcamento is None:
        orcamento = obj.Orcamento(tamanho_pop * geracoes)

    pop = rep.populacao_inicial(rng, tamanho_pop, sementes_smiles)
    arquivo = Arquivo()
    hist_melhor: list[float] = []
    hist_media: list[float] = []
    hist_aval: list[int] = []
    melhor_global = 0.0

    for g in range(geracoes):
        if orcamento.esgotado:
            break

        notas_pop = _avaliar(pop, peso_qed, orcamento)
        apt = [peso_qed * n.qed + (1 - peso_qed) * n.sa_norm for n in notas_pop]

        for cadeia, n, a in zip(pop, notas_pop, apt):
            if not n.valida:
                continue
            smi = rep.canonico(cadeia)
            if smi is None:
                continue
            arquivo.registrar(
                Achado(
                    smiles=smi, selfies=cadeia, qed=n.qed, sa=n.sa,
                    sa_norm=n.sa_norm, aptidao=a, geracao=g,
                    metodo="ga", semente=semente,
                )
            )

        melhor_global = max(melhor_global, max(apt))
        hist_melhor.append(melhor_global)
        hist_media.append(sum(apt) / len(apt))
        hist_aval.append(orcamento.usadas)

        # Elitismo: os e melhores passam intactos. Sem isto, uma mutacao
        # azarada pode destruir o melhor achado e a curva pode CAIR.
        ordem = sorted(range(len(pop)), key=lambda i: apt[i], reverse=True)
        nova = [pop[i] for i in ordem[:elites]]

        while len(nova) < tamanho_pop:
            p1 = ops.torneio(pop, apt, rng, k_torneio)
            p2 = ops.torneio(pop, apt, rng, k_torneio)
            if rng.random() < p_cruzamento:
                f1, f2 = ops.cruzar(p1, p2, rng)
            else:
                f1, f2 = p1, p2
            nova.append(ops.mutar(f1, rng, p_mutacao))
            if len(nova) < tamanho_pop:
                nova.append(ops.mutar(f2, rng, p_mutacao))

        pop = nova

    return Resultado(
        metodo="ga",
        semente=semente,
        arquivo=arquivo,
        historico_melhor=hist_melhor,
        historico_media=hist_media,
        historico_avaliacoes=hist_aval,
        populacao_final=pop,
    )
