#encoding: utf-8

import numpy as np

class Anticorpo:
    
    def __init__(self, qtd_alelos):
        self.alelos = np.zeros(qtd_alelos)
        self.aptidao = 0.0
