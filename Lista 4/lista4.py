import math
import numpy as np

tolerancia = 10**(-4)
numMaxDeIteracoes = 100

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
    if (numIteracoes == 99):
        print("Convergencia nao alcancada")
    return "Ponto em que f(x) = 0: " + str(xAtual) + "\nNumero de iteracoes: " + str(numIteracoes)

def derivada(funcao, valorDerivada):
    h = 10**(-10)
    topo = funcao(valorDerivada + h) - funcao(valorDerivada)
    inferior = h
    resultado = topo/inferior

    return resultado

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

    if (numIteracoes == 99):
        print("Convergencia nao alcancada")

    return "Ponto em que f(x) = 0: " + str(xProx) + "\nNumero de iteracoes: " + str(numIteracoes)

def secantHelper(xAnterior, xAtual, funcaoAnterior, funcaoAtual):
    return float((funcaoAtual*(xAtual - xAnterior))/(funcaoAtual - funcaoAnterior))

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

    return "Ponto em que f(x) = 0: " + str(xAtual) + "\nNumero de iteracoes: " + str(numIteracoes)

def inverseInterpolationHelper(x1, x2, x3, y1, y2, y3):
    a = (y2*y3*x1) / ((y1 - y2)*(y1 - y3))
    b = (y1*y3*x2) / ((y2 - y1)*(y2 - y3))
    c = (y1*y2*x3) / ((y3 - y1)*(y3 - y2))
    return a + b + c

def metodoDeNewtonParaEquacoesNaoLineares(listaDeFuncoes, vetorSolucaoInicial):
    xAtual = vetorSolucaoInicial
    for i in range(numMaxDeIteracoes):
        matrizJacobiano = criarMatrizJacobiano(listaDeFuncoes, xAtual)
        vetorF = criarVetorF(listaDeFuncoes, xAtual)
        inversaJacobiano = np.linalg.inv(matrizJacobiano)
        deltaX = multiplicacaoMatrizVetor(inversaJacobiano, vetorF)
        deltaX = [x*(-1) for x in deltaX]
        xAtual = somaDeVetores(xAtual, deltaX)
        tolk = norma(deltaX)/norma(xAtual)
        if (tolk < tolerancia):
            break
    return xAtual

def metodoDeBroyden(listaDeFuncoes, vetorSolucaoInicial):
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
    
    return xAtual


def criarMatrizJacobiano(listaDeFuncoes, vetorSolucaoInicial):
    dimensao = len(vetorSolucaoInicial)
    jacobiano = [[0 for i in range(dimensao)] for j in range(dimensao)]
    for i in range(dimensao):
        for j in range(dimensao):
            jacobiano[i][j] = derivadaParcial(listaDeFuncoes[i], vetorSolucaoInicial, j)
    return jacobiano

def derivadaParcial(funcaoDeVariasVariaveis, vetorSolucaoInicial, indiceArgumentoDerivada):
    h = 10**(-10)
    aux = funcaoDeVariasVariaveis(vetorSolucaoInicial)
    novoVetor = vetorSolucaoInicial
    novoVetor[indiceArgumentoDerivada] += h
    top = funcaoDeVariasVariaveis(novoVetor) - aux
    bottom = h
    return top/bottom

def criarVetorF(listaDeFuncoes, vetorSolucaoInicial):
    dimensao = len(vetorSolucaoInicial)
    vetorF = [0 for i in range(dimensao)]
    for i in range(dimensao):
        vetorF[i] = listaDeFuncoes[i](vetorSolucaoInicial)
    return vetorF

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

def somaDeVetores(vetor1, vetor2):
    dimensao = len(vetor1)
    resultado = [0 for i in range(dimensao)]
    for i in range(dimensao):
        resultado[i] = vetor1[i] + vetor2[i]
    return resultado

def norma(vetor):
    dimensao = len(vetor)
    norma = 0
    for i in range(dimensao):
        norma += vetor[i]*vetor[i]

    return norma**(0.5)

def broydenHelper(vetorA, vetorB):
    dimensao = len(vetorA)
    resultado = [[0 for i in range(dimensao)] for j in range(dimensao)]
    for i in range(dimensao):
        for j in range(dimensao):
            resultado[i][j] = vetorA[i]*vetorB[j]
    return resultado

def somaDeMatrizes(matrizA, matrizB):
    dimensao = len(matrizA)
    resultado = [[0 for i in range(dimensao)] for j in range(dimensao)]
    for i in range(dimensao):
        for j in range(dimensao):
            resultado[i][j] = matrizA[i][j] + matrizB[i][j]

    return resultado

def multiplicacaoVetores(A, B):
    s = 0
    for i in range(len(A)):
        for j in range(len(B)):
            if (i==j):
                s += A[i]*B[i]
    return s

def funcaoDeVariasVariaveis1(listaDeVariaveis):
    return listaDeVariaveis[0] + 2*listaDeVariaveis[1] - 2

def funcaoDeVariasVariaveis2(listaDeVariaveis):
    return listaDeVariaveis[0]*listaDeVariaveis[0] + 4*listaDeVariaveis[1]*listaDeVariaveis[1] - 4

def funcao(x):
    return x*x - 4*math.cos(x)

# Retorna a transposta de uma matriz
def transposta(matriz):
    n = len(matriz)
    resp = [[0.0]*n for i in range(n)]
    for i in range(n):
        for j in range(n):
            resp[j][i] = matriz[i][j]
    return resp

# Multiplica duas matrizes
def multiplyMatrices(X, Y):
    n = len(X)
    result = [[0.0 for j in range(n)] for i in range(n)]
    for i in range(len(X)):
        for j in range(len(Y[0])):
            for k in range(len(Y)):
                result[i][j] += X[i][k] * Y[k][j]
    return result

def multiplicacaoMatrizEscalar(X, escalar):
    n = len(X)
    resultado = [[0.0 for j in range(n)] for i in range(n)]
    for i in range(len(X)):
        for j in range (len(X))
            resultado[i][j] = X[i][j] * escalar
    print resultado



def minimosQuadrados(listaDeFuncoes, vetorSolucaoInicial):
    dimensao = len(vetorSolucaoInicial)
    xAtual = vetorSolucaoInicial
    for i in range (numMaxDeIteracoes):
        jacobiano = criarMatrizJacobiano(listaDeFuncoes, vetorSolucaoInicial)
        jacobianoTransposto = transposta(jacobiano)
        vetorF = criarVetorF(listaDeFuncoes,vetorSolucaoInicial)
        a = np.linalg.inv(multiplyMatrices(jacobianoTransposto,jacobiano))
        b = multiplicacaoMatrizVetor(a,-1)
        c = multiplicacaoMatrizVetor(jacobianoTransposto,vetorF)
        deltaB = multiplyMatrices(b,c)
        xAtual = somaDeVetores(xAtual, deltaX)
        tolk = norma(deltaX)/norma(xAtual)
        if (tolk < tolerancia):
            break
    return xAtual


#print(metodoDaBissecao(0,10))
#print(metodoDeNewton(10))
#print(metodoSecante(10))
#print(metodoDaInterpolacaoInversa(3,5,10))
#print(derivadaParcial(funcaoDeVariasVariaveis2, [2,3], 1))
#print(funcaoDeVariasVariaveis2([2,3]))
#criarMatrizJacobiano([funcaoDeVariasVariaveis1, funcaoDeVariasVariaveis2], [2,3])
#criarVetorF([funcaoDeVariasVariaveis1, funcaoDeVariasVariaveis2], [2,3])
#print(metodoDeNewtonParaEquacoesNaoLineares([funcaoDeVariasVariaveis1, funcaoDeVariasVariaveis2], [2,3]))
<<<<<<< HEAD
#print(minimosQuadrados(,)[0,1])
=======
print(metodoDeBroyden([funcaoDeVariasVariaveis1, funcaoDeVariasVariaveis2], [2,3]))
>>>>>>> 389dcfc7cc6ff0b73f9545793fb03d3b5301feeb
