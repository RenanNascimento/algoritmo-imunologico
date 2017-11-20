#encoding: utf-8

from selecao_clonal import SelecaoClonal

tam_populacao, inicio, fim, constante_clones, tax_clonagem, ro, max_it = (200, -512, 512, 50, 0.1, 2, 100)

sc = SelecaoClonal(tam_populacao, inicio, fim,
                   constante_clones, tax_clonagem, ro, max_it)
sc.executar()
