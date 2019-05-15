import math

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

def metodoDeNewtonParaEquacoesNaoLineares(vetorSolucaoInicial):
    pass
    
def funcao(x):
    return x*x - 4*math.cos(x)

#print(metodoDaBissecao(0,10))
#print(metodoDeNewton(10))
#print(metodoSecante(10))
#print(metodoDaInterpolacaoInversa(3,5,10))