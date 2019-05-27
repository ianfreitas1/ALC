import math
import numpy as np

tolerancia = 10**(-4)
numMaxDeIteracoes = 100

#Funcao que encontra uma unica raiz de uma equaçao nao linear pelo metodo da bissecao
#Os argumentos intervalMin e intervalMax sao os extremos do intervalo conhecido [a, b]
def metodoDaBissecao(intervalMin, intervalMax):
    numIteracoes = 0
    diferencaDoIntervalo = abs(intervalMax - intervalMin)
    while diferencaDoIntervalo > tolerancia:
        pontoMedio = float((intervalMax + intervalMin)/2)
        y = funcao(pontoMedio)
        if (y > 0.0):
            intervalMax = pontoMedio
        else:
            intervalMin = pontoMedio
        diferencaDoIntervalo = abs(intervalMax - intervalMin)
        numIteracoes += 1
    return "Ponto em que f(x) = 0: " + str(pontoMedio) + "\nNumero de iteracoes: " + str(numIteracoes) + "\n"

#Funcao que encontra uma unica raiz de uma equaçao nao linear pelo metodo de Newton
#Chute inicial eh um numero inteiro
def metodoDeNewton(chuteInicial):
    numIteracoes = 0
    xAnterior = chuteInicial
    for i in range(numMaxDeIteracoes):
        xAtual = xAnterior - funcao(xAnterior)/derivada(funcao, xAnterior)
        diferencaDosValores = abs(xAtual - xAnterior)
        if (diferencaDosValores > tolerancia):
            xAnterior = xAtual
            numIteracoes += 1
        else:
            break
    if (numIteracoes == numMaxDeIteracoes-1):
        print("Convergencia nao alcancada")
    return "Ponto em que f(x) = 0: " + str(xAtual) + "\nNumero de iteracoes: " + str(numIteracoes)

def derivada(funcao, valorDerivada):
    h = 10**(-10)
    topo = funcao(valorDerivada + h) - funcao(valorDerivada)
    inferior = h
    resultado = topo/inferior

    return resultado

#Funcao que encontra uma unica raiz de uma equaçao nao linear pelo metodo Secante
#Chute inicial eh um numero inteiro
def metodoSecante(chuteInicial):
    numIteracoes = 0
    deltaX = 0.001
    xAnterior = chuteInicial
    xAtual = xAnterior + deltaX
    funcaoAnterior = funcao(xAnterior)
    for i in range(numMaxDeIteracoes):
        funcaoAtual = funcao(xAtual)
        xProx = xAtual - secantHelper(xAnterior, xAtual, funcaoAnterior, funcaoAtual)
        diferencaDosValores = abs(xProx - xAtual)
        if (diferencaDosValores > tolerancia):
            funcaoAnterior = funcaoAtual
            xAnterior = xAtual
            xAtual = xProx
            numIteracoes += 1
        else:
            break

    if (numIteracoes == numMaxDeIteracoes-1):
        print("Convergencia nao alcancada")

    return "Ponto em que f(x) = 0: " + str(xProx) + "\nNumero de iteracoes: " + str(numIteracoes)

#Funcao auxiliar do metodo Secante utilizado em cada iteracao para calculo do proximo valor de X
def secantHelper(xAnterior, xAtual, funcaoAnterior, funcaoAtual):
    return float((funcaoAtual*(xAtual - xAnterior))/(funcaoAtual - funcaoAnterior))

#Funcao que encontra uma unica raiz de uma equaçao nao linear pelo metodo da interpolacao inversa
#x1, x2 e x3 sao numeros inteiros em que x1 < x2 < x3
def metodoDaInterpolacaoInversa(x1, x2, x3):
    numIteracoes = 0
    xAnterior = 10**(36)
    x = [x1, x2, x3]
    y1 = funcao(x1)
    y2 = funcao(x2)
    y3 = funcao(x3)
    y = [y1, y2, y3]
    for i in range(numMaxDeIteracoes):
        xAtual = inverseInterpolationHelper(x[0], x[1], x[2], y[0], y[1], y[2])
        diferencaDosValores = abs(xAtual - xAnterior)
        print(diferencaDosValores)
        if (diferencaDosValores > tolerancia):
            yMax = max(y)
            i = y.index(yMax)
            x[i] = xAtual
            x.sort()
            y[i] = funcao(xAtual)
            y.sort()
            xAnterior = xAtual
            numIteracoes += 1
        else:
            break

    if (numIteracoes == numMaxDeIteracoes-1):
        print("Convergencia nao alcancada")

    return "Ponto em que f(x) = 0: " + str(xAtual) + "\nNumero de iteracoes: " + str(numIteracoes)

#Funcao auxiliar do metodo da interpolacao inversa utilizada em cada iteracao para calculo do X
def inverseInterpolationHelper(x1, x2, x3, y1, y2, y3):
    a = (y2*y3*x1) / ((y1 - y2)*(y1 - y3))
    b = (y1*y3*x2) / ((y2 - y1)*(y2 - y3))
    c = (y1*y2*x3) / ((y3 - y1)*(y3 - y2))
    return a + b + c

#Solucao de um sistema de equacoes nao lineares pelo metodo de Newton
def metodoDeNewtonParaEquacoesNaoLineares(listaDeFuncoes, vetorSolucaoInicial):
    xAtual = vetorSolucaoInicial
    numIteracoes = 0
    for i in range(numMaxDeIteracoes):
        matrizJacobiano = criarMatrizJacobiano(listaDeFuncoes, xAtual)
        vetorF = criarVetorF(listaDeFuncoes, xAtual)
        inversaJacobiano = np.linalg.inv(matrizJacobiano)
        deltaX = multiplicacaoMatrizVetor(inversaJacobiano, vetorF)
        deltaX = [x*(-1) for x in deltaX]
        xAtual = somaDeVetores(xAtual, deltaX)
        tolk = norma(deltaX)/norma(xAtual)
        numIteracoes += 1
        if (tolk < tolerancia):
            break

    if (numIteracoes == numMaxDeIteracoes-1):
        print("Convergencia nao alcancada")

    return xAtual

#Solucao de um sistema de equacoes nao lineares pelo metodo de Broyden
def metodoDeBroyden(listaDeFuncoes, vetorSolucaoInicial):
    numIteracoes = 0
    dimensao = len(vetorSolucaoInicial)
    xAtual = vetorSolucaoInicial
    jacobianoB = criarMatrizJacobiano(listaDeFuncoes, vetorSolucaoInicial)
    jacobiano = [[0 for i in range(dimensao)] for j in range(dimensao)]
    for i in range(numMaxDeIteracoes):
        for k in range(dimensao):
            for j in range(dimensao):
                jacobiano[k][j] = jacobianoB[k][j]
        vetorF = criarVetorF(listaDeFuncoes, xAtual)
        inversaJacobiano = np.linalg.inv(jacobiano)
        deltaX = multiplicacaoMatrizVetor(inversaJacobiano, vetorF)
        deltaX = [x*(-1) for x in deltaX]
        xAtual = somaDeVetores(xAtual, deltaX)
        vetorF2 = criarVetorF(listaDeFuncoes, xAtual)   
        vetorFnegativo = [x*(-1) for x in vetorF]
        yK = somaDeVetores(vetorF2, vetorFnegativo)
        tolk = norma(deltaX)/norma(xAtual)
        if (tolk < tolerancia):
            break
        else:
            resultMultiplicacaoJacobianoDeltaX = multiplicacaoMatrizVetor(jacobianoB, deltaX)
            resultMultiplicacaoJacobianoDeltaX = [x*(-1) for x in resultMultiplicacaoJacobianoDeltaX]
            bottom = multiplicacaoVetores(deltaX, deltaX)
            top = somaDeVetores(yK, resultMultiplicacaoJacobianoDeltaX)
            top3 = broydenHelper(top, deltaX)
            for k in range(len(jacobianoB)):
                for j in range(len(jacobianoB)):
                    jacobianoB[k][j] = top3[k][j]/bottom
            jacobianoB = somaDeMatrizes(jacobiano, jacobianoB)
    
    if (numIteracoes == numMaxDeIteracoes-1):
        print("Convergencia nao alcancada")
    return xAtual

#Cria a matriz jacobiana, necessaria para os metodos de Newton e Broyden
def criarMatrizJacobiano(listaDeFuncoes, vetorSolucaoInicial):
    dimensao1 = len(listaDeFuncoes)
    dimensao2 = len(vetorSolucaoInicial)
    jacobiano = [[0 for i in range(dimensao2)] for j in range(dimensao1)]
    for i in range(dimensao1):
        for j in range(dimensao2):
            jacobiano[i][j] = derivadaParcial(listaDeFuncoes[i], vetorSolucaoInicial, j)
    return jacobiano

def derivadaParcial(funcaoDeVariasVariaveis, vetorSolucaoInicial, indiceArgumentoDerivada):
    h = 10**(-10)
    aux = funcaoDeVariasVariaveis(vetorSolucaoInicial)
    novoVetor = vetorSolucaoInicial[:]
    novoVetor[indiceArgumentoDerivada] += h
    top = funcaoDeVariasVariaveis(novoVetor) - aux
    bottom = h
    return top/bottom

def criarVetorF(listaDeFuncoes, vetorSolucaoInicial):
    dimensao = len(listaDeFuncoes)
    vetorF = [0 for i in range(dimensao)]
    for i in range(dimensao):
        vetorF[i] = listaDeFuncoes[i](vetorSolucaoInicial)
    return vetorF

def norma(vetor):
    dimensao = len(vetor)
    norma = 0
    for i in range(dimensao):
        norma += vetor[i]*vetor[i]

    return norma**(0.5)

#Funcao auxiliar do metodo de Broyden que multiplica dois vetores
def broydenHelper(vetorA, vetorB):
    dimensao = len(vetorA)
    resultado = [[0 for i in range(dimensao)] for j in range(dimensao)]
    for i in range(dimensao):
        for j in range(dimensao):
            resultado[i][j] = vetorA[i]*vetorB[j]
    return resultado

def funcaoDeVariasVariaveis1(listaDeVariaveis):
    return listaDeVariaveis[0] + 2*listaDeVariaveis[1] - 2

def funcaoDeVariasVariaveis2(listaDeVariaveis):
    return listaDeVariaveis[0]*listaDeVariaveis[0] + 4*listaDeVariaveis[1]*listaDeVariaveis[1] - 4

def funcao(x):
    return x*x - 4*math.cos(x)

#Ajuste de uma funcao pelo metodo dos minimos quadrados
#listaDePontosXY é uma lista de listas que contem os pontos X e os respectivos valores de Y
def minimosQuadrados(listaDePontosXY, vetorSolucaoInicialB):
    numIteracoes = 0
    xAtual = vetorSolucaoInicialB
    listaDeFuncoes = [funcaoRegressao1, funcaoRegressao2, funcaoRegressao3]
    for i in range (numMaxDeIteracoes):
        jacobiano = criarMatrizJacobiano(listaDeFuncoes, xAtual)
        jacobianoTransposto = map(list, zip(*jacobiano))
        #jacobianoTransposto = transposta(jacobiano)
        vetorF = criarVetorF(listaDeFuncoes,xAtual)
        a = np.linalg.inv(multiplyMatrices(jacobianoTransposto,jacobiano))
        b = multiplicacaoMatrizEscalar(a,-1)
        c = multiplicacaoMatrizVetor(jacobianoTransposto,vetorF)
        deltaB = multiplicacaoMatrizVetor(b,c)
        xAtual = somaDeVetores(xAtual, deltaB)
        tolk = norma(deltaB)/norma(xAtual)
        numIteracoes += 1
        if (tolk < tolerancia):
            break
    if (numIteracoes == numMaxDeIteracoes-1):
        print("Convergencia nao alcancada")
    return xAtual

def funcaoRegressao1(listaDeVariaveis):
    ponto = [1, 1.995]
    return math.exp((ponto[0]**listaDeVariaveis[0])/listaDeVariaveis[1]) - ponto[1]

def funcaoRegressao2(listaDeVariaveis):
    ponto = [2, 1.410]
    return math.exp((ponto[0]**listaDeVariaveis[0])/listaDeVariaveis[1]) - ponto[1]

def funcaoRegressao3(listaDeVariaveis):
    ponto = [3, 1.260]
    return math.exp((ponto[0]**listaDeVariaveis[0])/listaDeVariaveis[1]) - ponto[1]

# Retorna a transposta de uma matriz
def transposta(matriz):
    n1 = len(matriz) #3
    n2 = len(matriz[0]) #2
    resp = [[0.0]*n2 for i in range(n1)]
    for i in range(n2):
        for j in range(n1):
            resp[j][i] = matriz[i][j]
    return resp

# Multiplica duas matrizes
def multiplyMatrices(X, Y):
    result = [[0.0 for i in range(len(X))] for j in range(len(Y[0]))]
    for i in range(len(X)):
    # iterate through columns of Y
        for j in range(len(Y[0])):
        # iterate through rows of Y
            for k in range(len(Y)):
                result[i][j] += X[i][k] * Y[k][j]
    return result

# Multiplica uma matriz por um vetor
def multiplicacaoMatrizVetor(m, v):
    rows = len(m)
    w = [0]*rows
    irange = range(len(v))
    sum = 0
    for j in range(rows):
        r = m[j]
        for i in irange:
            sum += r[i]*v[i]
        w[j],sum = sum,0
    return w

# Soma de duas matrizes
def somaDeMatrizes(matrizA, matrizB):
    dimensao = len(matrizA)
    resultado = [[0 for i in range(dimensao)] for j in range(dimensao)]
    for i in range(dimensao):
        for j in range(dimensao):
            resultado[i][j] = matrizA[i][j] + matrizB[i][j]

    return resultado

# Multiplica dois vetores
def multiplicacaoVetores(A, B):
    s = 0
    for i in range(len(A)):
        for j in range(len(B)):
            if (i==j):
                s += A[i]*B[i]
    return s

# Soma dois vetores
def somaDeVetores(vetor1, vetor2):
    dimensao = len(vetor1)
    resultado = [0 for i in range(dimensao)]
    for i in range(dimensao):
        resultado[i] = vetor1[i] + vetor2[i]
    return resultado

# Multiplica matriz por um escalar
def multiplicacaoMatrizEscalar(X, escalar):
    n = len(X)
    resultado = [[0.0 for j in range(n)] for i in range(n)]
    for i in range(len(X)):
        for j in range (len(X)):
            resultado[i][j] = X[i][j] * escalar
    return resultado

#print(metodoDaBissecao(0,10))
#print(metodoDeNewton(10))
#print(metodoSecante(10))
#print(metodoDaInterpolacaoInversa(3,5,10))
#print(derivadaParcial(funcaoDeVariasVariaveis2, [2,3], 1))
#print(funcaoDeVariasVariaveis2([2,3]))
#criarMatrizJacobiano([funcaoDeVariasVariaveis1, funcaoDeVariasVariaveis2], [2,3])
#criarVetorF([funcaoDeVariasVariaveis1, funcaoDeVariasVariaveis2], [2,3])
#print(metodoDeNewtonParaEquacoesNaoLineares([funcaoDeVariasVariaveis1, funcaoDeVariasVariaveis2], [2,3]))
#print(minimosQuadrados([[1,1.995],[2,1.410],[3,1.260]], [0,1]))
#print(metodoDeBroyden([funcaoDeVariasVariaveis1, funcaoDeVariasVariaveis2], [2,3]))
#print(derivadaParcial(funcaoRegressao1, [0,1], 1))

