#encoding: utf-8

import numpy as np
import random
import copy
import math
import matplotlib.pyplot as plt

from anticorpo import Anticorpo


class SelecaoClonal:

    def __init__(self, tam_populacao, inicio, fim, constante_clones, tax_clonagem, ro, max_it):
        self.tam_populacao = tam_populacao
        self.inicio = inicio
        self.fim = fim
        self.constante_clones = constante_clones
        self.tax_clonagem = tax_clonagem
        self.ro = ro
        self.max_it = max_it
        self.QTD_ALELOS = 2
        self.populacao = []
        self.clones = []
        self.OFFSET = 1500
        self.delta_init = -2
        self.delta_end = 2

    '''
    ' Inicializa a lista de população
    '''

    def popular(self):
        for c in range(0, self.tam_populacao):
            ant = Anticorpo(self.QTD_ALELOS)
            for a in range(0, len(ant.alelos)):
                ant.alelos[a] = random.randrange(self.inicio, self.fim + 1)
            # Calcula a aptidao do anticorpo
            ant.aptidao = self.calcAptidao(ant.alelos[0], ant.alelos[1])
            self.populacao.append(ant)

    '''
    ' Maior aptidao da populacao
    '''

    def maiorAptidao(self, populacao):
        populacao_ordenada = sorted(
            populacao, reverse=True, key=lambda Anticorpo: Anticorpo.aptidao)
        return populacao_ordenada[0]

    '''
    ' Calcula a função de aptidão de toda a população
    '''

    def calcAptidao(self, x1, x2):
        return self.OFFSET + (x2 + 47) * np.sin(np.sqrt(abs(x2 + x1 / 2 + 47))) + x1 * np.sin(np.sqrt(abs(x1 - (x2 + 47))))

    '''
    ' Calcula a quantidade de clones a ser gerada para cada anticorpo
    '''

    def qtdeClones(self):
        return int(self.tax_clonagem * self.constante_clones)

    '''
    ' Realiza a mutacao de um clone de um anticorpo
    '''

    def maturacao(self, ant, clones):
        afinidade_normalizada = ant.aptidao / \
            self.maiorAptidao(self.populacao).aptidao
        tax_mutacao = math.exp(-self.ro * afinidade_normalizada)
        clones_mutados = []

        for clone in clones:
            for i in range(0, len(clone.alelos)):
                # Sorteia valor para decidir se o alelo sera mutado
                num = random.uniform(0, 1)
                # Muta o alelo
                if num <= tax_mutacao:
                    delta = random.uniform(self.delta_init, self.delta_end)
                    clone.alelos[i] += delta
                    # Poda valores maiores que o dominio estipulado da funcao
                    if clone.alelos[i] < self.inicio:
                        clone.alelos[i] = self.inicio
                    elif clone.alelos[i] > self.fim:
                        clone.alelos[i] = self.fim
            # Calcula aptidao do clone e o insere na lista de clones mutados
            clone.aptidao = self.calcAptidao(clone.alelos[0], clone.alelos[1])
            clones_mutados.append(clone)

        return clones_mutados

    '''
    ' Seleciona o melhor clone para se tornar o novo anticorpo da populacao
    '''

    def selecao(self):
        # Retira todos os individuos da populacao atual
        # ja que a populacao sera formada
        # pelos melhores clones dos clusters
        self.populacao[:] = []

        for clones in self.clones:
            self.populacao.append(self.maiorAptidao(clones))

    '''
    ' Todos os anticorpos sao selecionados e gerados qtdeClones para cada
    '''

    def clonagem(self):
        qtde_clones = self.qtdeClones()
        self.clones[:] = []

        for ant in self.populacao:
            ant_clones = []
            for i in range(0, qtde_clones):
                clone = copy.deepcopy(ant)
                clone.aptidao = 0.0
                ant_clones.append(clone)
            # Realiza a maturacao dos clones de um anticorpo
            clones_mutados = self.maturacao(ant, ant_clones)
            self.clones.append(clones_mutados)

    '''
    ' Retorna mínimo, média e máximo
    '''

    def retornaMedidas(self):
        total_aptidao = 0
        minimo = self.populacao[0].aptidao
        maximo = self.populacao[0].aptidao
        for ant in self.populacao:
            total_aptidao += ant.aptidao
            if ant.aptidao < minimo:
                minimo = ant.aptidao
            elif ant.aptidao > maximo:
                maximo = ant.aptidao
        return minimo, total_aptidao / self.tam_populacao, maximo

    '''
    ' Gera o gráfico Erro x Iteração
    '''

    def gerarGrafico(self, tempo, minimo, media, maximo):
        plt.plot(tempo, minimo, label='Mínimo')
        plt.plot(tempo, media, label='Média')
        plt.plot(tempo, maximo, label='Máximo')
        plt.ylabel('Aptidão')
        plt.xlabel('Iteração')
        plt.title('Gráfico Aptidão x Iteração')
        plt.legend()
        plt.show()

    '''
    ' Melhores anticorpos daquela populacao
    '''

    def melhoresAnticorpos(self):
        n_melhores = 3
        melhores = sorted(self.populacao, reverse=True,
                          key=lambda Anticorpo: Anticorpo.aptidao)
        print('Melhores anticorpos:')
        for ant in melhores[:n_melhores]:
            print(str(ant.alelos) + '\t=\t' + str(ant.aptidao - self.OFFSET))

    '''
    ' Executa algoritmo imunologico
    '''

    def executar(self):
        self.popular()
        # Inicaliza os vetores de tempo e erro com o estado inicial
        mini, med, maxi = self.retornaMedidas()
        tempo = [0]
        minimo = [mini]
        media = [med]
        maximo = [maxi]
        t = 1
        # Executa o algoritmo até atinger o máximo de iteraçẽos
        while t <= self.max_it:
            self.clonagem()
            self.selecao()
            mini, med, maxi = self.retornaMedidas()
            tempo.append(t)
            minimo.append(mini)
            media.append(med)
            maximo.append(maxi)
            t += 1
        # Gera o gráfico de saída
        self.gerarGrafico(tempo, minimo, media, maximo)
        self.melhoresAnticorpos()
