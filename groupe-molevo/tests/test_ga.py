"""Testes do algoritmo genetico, das linhas de base e do verificador.

Rodadas pequenas de proposito: estes testes existem para provar
propriedades (o elitismo funciona, o orcamento e respeitado, a semente
reproduz), nao para medir qualidade. As medidas vao no notebook.
"""

import json

import pytest

from src import baselines, ga
from src import objectives as obj
from src import verifier


PEQUENO = dict(tamanho_pop=12, geracoes=6, peso_qed=0.5)


# --------------------------------------------------------------------------
# Elitismo e historico
# --------------------------------------------------------------------------

def test_melhor_nunca_piora():
    """A consequencia observavel do elitismo.

    Sem elitismo, uma mutacao azarada destroi o melhor individuo e a curva
    CAI. Este teste e a forma barata de provar que ele esta ligado.
    """
    r = ga.evoluir(semente=0, **PEQUENO)
    historico = r.historico_melhor
    assert historico == sorted(historico), historico


def test_historico_tem_uma_entrada_por_geracao():
    r = ga.evoluir(semente=0, **PEQUENO)
    assert len(r.historico_melhor) == PEQUENO["geracoes"]
    assert len(r.historico_media) == PEQUENO["geracoes"]
    assert len(r.historico_avaliacoes) == PEQUENO["geracoes"]


def test_avaliacoes_crescem_de_forma_monotonica():
    r = ga.evoluir(semente=0, **PEQUENO)
    assert r.historico_avaliacoes == sorted(r.historico_avaliacoes)


def test_populacao_final_tem_o_tamanho_certo():
    r = ga.evoluir(semente=0, **PEQUENO)
    assert len(r.populacao_final) == PEQUENO["tamanho_pop"]


# --------------------------------------------------------------------------
# Reprodutibilidade
# --------------------------------------------------------------------------

def test_mesma_semente_mesmo_resultado():
    obj.limpar_cache()
    a = ga.evoluir(semente=42, **PEQUENO)
    obj.limpar_cache()
    b = ga.evoluir(semente=42, **PEQUENO)
    assert a.historico_melhor == b.historico_melhor
    assert a.populacao_final == b.populacao_final


def test_sementes_diferentes_resultados_diferentes():
    a = ga.evoluir(semente=1, **PEQUENO)
    b = ga.evoluir(semente=2, **PEQUENO)
    assert a.populacao_final != b.populacao_final


# --------------------------------------------------------------------------
# Orcamento
# --------------------------------------------------------------------------

def test_ga_respeita_o_orcamento():
    o = obj.Orcamento(50)
    ga.evoluir(semente=0, tamanho_pop=20, geracoes=100, orcamento=o)
    # Para na primeira geracao apos esgotar, logo nao passa de um lote.
    assert o.usadas <= 50 + 20


def test_os_tres_metodos_gastam_o_mesmo_orcamento():
    """A comparacao do relatorio so e honesta se isto valer."""
    total = 200
    o_al = obj.Orcamento(total)
    o_hc = obj.Orcamento(total)
    baselines.busca_aleatoria(0, o_al)
    baselines.hill_climbing(0, o_hc)
    assert o_al.usadas == total
    assert o_hc.usadas == total


# --------------------------------------------------------------------------
# Arquivo de achados
# --------------------------------------------------------------------------

def test_arquivo_nao_guarda_duplicatas():
    r = ga.evoluir(semente=0, **PEQUENO)
    smiles = [a.smiles for a in r.arquivo.itens.values()]
    assert len(smiles) == len(set(smiles))


def test_arquivo_guarda_mais_do_que_a_populacao_final():
    """O ponto do arquivo: preservar o que foi achado no caminho."""
    r = ga.evoluir(semente=0, tamanho_pop=20, geracoes=15, peso_qed=0.5)
    assert len(r.arquivo) > 20


def test_melhores_vem_ordenados():
    r = ga.evoluir(semente=0, **PEQUENO)
    aptidoes = [a.aptidao for a in r.arquivo.melhores(5)]
    assert aptidoes == sorted(aptidoes, reverse=True)


# --------------------------------------------------------------------------
# Linhas de base
# --------------------------------------------------------------------------

@pytest.mark.parametrize("metodo", ["aleatoria", "hill"])
def test_linhas_de_base_produzem_resultado_utilizavel(metodo):
    o = obj.Orcamento(300)
    r = (baselines.busca_aleatoria if metodo == "aleatoria" else baselines.hill_climbing)(0, o)
    assert r.historico_melhor
    assert r.historico_melhor == sorted(r.historico_melhor)
    assert len(r.arquivo) > 0
    assert r.historico_avaliacoes[-1] == 300


def test_hill_climbing_com_reinicios_tambem_funciona():
    o = obj.Orcamento(300)
    r = baselines.hill_climbing(0, o, reinicios=True, paciencia=20)
    assert r.melhor_aptidao > 0
    assert "reinicios" in r.metodo


# --------------------------------------------------------------------------
# Verificador
# --------------------------------------------------------------------------

def test_verificador_aprova_um_registro_honesto():
    r = ga.evoluir(semente=0, **PEQUENO)
    registro = [a.como_dict() for a in r.arquivo.itens.values()]
    assert verifier.verificar(registro) == []


def test_verificador_pega_nota_adulterada():
    """O erro que o verificador existe para encontrar."""
    r = ga.evoluir(semente=0, **PEQUENO)
    registro = [a.como_dict() for a in r.arquivo.itens.values()]
    registro[0]["qed"] = 0.99
    problemas = verifier.verificar(registro)
    assert any(p.tipo == "qed_divergente" for p in problemas)


def test_verificador_pega_duplicata():
    r = ga.evoluir(semente=0, **PEQUENO)
    registro = [a.como_dict() for a in r.arquivo.itens.values()]
    registro.append(dict(registro[0]))
    assert any(p.tipo == "duplicado" for p in verifier.verificar(registro))


def test_verificador_pega_smiles_ilegivel():
    problemas = verifier.verificar([{"smiles": "isto nao e um smiles"}])
    assert any(p.tipo == "smiles_ilegivel" for p in problemas)


def test_verificador_pega_frente_dominada():
    frente_errada = [
        {"smiles": "CCO", "qed": 0.8, "sa_norm": 0.8},
        {"smiles": "CCC", "qed": 0.5, "sa_norm": 0.5},   # dominada pela primeira
    ]
    problemas = verifier.verificar_frente(frente_errada)
    assert any(p.tipo == "frente_dominada" for p in problemas)


def test_verificador_aceita_frente_legitima():
    frente_ok = [
        {"smiles": "CCO", "qed": 0.8, "sa_norm": 0.2},
        {"smiles": "CCC", "qed": 0.2, "sa_norm": 0.8},
    ]
    assert verifier.verificar_frente(frente_ok) == []


def test_registro_sobrevive_a_ida_e_volta_pelo_json(tmp_path):
    r = ga.evoluir(semente=0, **PEQUENO)
    caminho = tmp_path / "registro.jsonl"
    with open(caminho, "w", encoding="utf-8") as f:
        for a in r.arquivo.itens.values():
            f.write(json.dumps(a.como_dict(), ensure_ascii=False) + "\n")
    assert verifier.verificar(verifier.carregar_jsonl(str(caminho))) == []
