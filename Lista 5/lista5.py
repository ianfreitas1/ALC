import sys
sys.path.append('../Lista 2/')
from lista2 import jacobiIterativo
import math

def integracaoPolinomial(funcao, a, b, numPontos):
    deltaX = float(abs(b-a)) / (numPontos - 1)
    pontos = []
    vetorB = []
    for i in range(1,numPontos+1):
        pontos.append(a+(i-1)*deltaX)
        vetorB.append(float((b**i-a**i))/i)
    matrizDeVandermonde = [[0.0 for i in range(numPontos)] for j in range(numPontos)]
    for i in range(numPontos):
        for j in range(numPontos):
            matrizDeVandermonde[i][j] = pontos[j]**i
    pesos = jacobiIterativo(matrizDeVandermonde, vetorB)
    soma = 0
    for i in range(numPontos):
        soma += pesos[i]*funcao(pontos[i])
    print(soma)

def funcao(x):
    return math.exp(-x*x)

integracaoPolinomial(funcao, 0, 1, 5)