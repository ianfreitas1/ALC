def metodoEuler(funcaoDiferencial, tInicial, tFinal, deltaT, resultadoCondicaoInicial):
    incognitasX = [resultadoCondicaoInicial]
    incognitasT = [tInicial]
    numIteracoes = int((tFinal - tInicial)/(deltaT))
    for i in range(numIteracoes):
        incognitasT.append((i+1)*deltaT)
        incognitasX.append(incognitasX[i]+deltaT*funcaoDiferencial(incognitasT[i], incognitasX[i]))
    return incognitasX

def rungeKutta2(funcaoDiferencial, tInicial, tFinal, deltaT, resultadoCondicaoInicial):
    incognitasX = [resultadoCondicaoInicial]
    incognitasT = [tInicial]
    numIteracoes = int((tFinal - tInicial)/deltaT)
    for i in range(numIteracoes):
        incognitasT.append((i+1)*deltaT)
        K1 = funcaoDiferencial(incognitasT[i], incognitasX[i])
        K2 = funcaoDiferencial(incognitasT[i] + deltaT, incognitasX[i] + deltaT*K1)
        incognitasX.append(incognitasX[i] + deltaT/2*(K1 + K2))
    return incognitasX

def rungeKutta4(funcaoDiferencial, tInicial, tFinal, deltaT, resultadoCondicaoInicial):
    incognitasX = [resultadoCondicaoInicial]
    incognitasT = [tInicial]
    numIteracoes = int((tFinal - tInicial)/deltaT)
    for i in range(numIteracoes):
        incognitasT.append((i+1)*deltaT)
        K1 = funcaoDiferencial(incognitasT[i], incognitasX[i])
        K2 = funcaoDiferencial(incognitasT[i] + deltaT/2, incognitasX[i] + deltaT/2*K1)
        K3 = funcaoDiferencial(incognitasT[i] + deltaT/2, incognitasX[i] + deltaT/2*K2)
        K4 = funcaoDiferencial(incognitasT[i] + deltaT, incognitasX[i] + deltaT*K3)
        incognitasX.append(incognitasX[i] + deltaT/6*(K1 + 2*K2 + 2*K3 + K4))
    return incognitasX

def taylor2Ordem(funcaoDiferencial2Ordem, tInicial, tFinal, deltaT, xInicial, xLinhaInicial):
    incognitasX = [xInicial]
    xLinhaAtual = xLinhaInicial
    incognitasT = [tInicial]
    numIteracoes = int((tFinal - tInicial)/deltaT)
    for i in range(numIteracoes):
        incognitasT.append((i+1)*deltaT)
        x2Linhas = funcaoDiferencial2Ordem(incognitasT[i], incognitasX[i], xLinhaAtual)
        incognitasX.append(incognitasX[i] + xLinhaAtual*deltaT + x2Linhas*deltaT*deltaT/2)
        xLinhaAtual = xLinhaAtual + x2Linhas*deltaT

    return incognitasX

def rungeKuttaNystron(funcaoDiferencial2Ordem, tInicial, tFinal, deltaT, xInicial, xLinhaInicial):
    incognitasX = [xInicial]
    xLinhaAtual = xLinhaInicial
    incognitasT = [tInicial]
    numIteracoes = int((tFinal - tInicial)/deltaT)
    for i in range(numIteracoes):
        incognitasT.append((i+1)*deltaT)
        K1 = deltaT/2 * funcaoDiferencial2Ordem(incognitasT[i], incognitasX[i], xLinhaAtual)
        Q = deltaT/2 * (xLinhaAtual + K1/2)
        K2 = deltaT/2 * funcaoDiferencial2Ordem(incognitasT[i] + deltaT/2, incognitasX[i] + Q, xLinhaAtual + K1)
        K3 = deltaT/2 * funcaoDiferencial2Ordem(incognitasT[i] + deltaT/2, incognitasX[i] + Q, xLinhaAtual + K2)
        L = deltaT * (xLinhaAtual + K3)
        K4 = deltaT/2 * funcaoDiferencial2Ordem(incognitasT[i] + deltaT, incognitasX[i] + L, xLinhaAtual + 2*K3)
        incognitasX.append(incognitasX[i] + deltaT*(xLinhaAtual + (K1 + K2 + K3)/3))
        xLinhaAtual = xLinhaAtual + (K1 + 2*K2 + 2*K3 + K4)/3

    return incognitasX

def funcaoDiferencial(t, funcao):
    return t + funcao

def funcaoDiferencial2Ordem(t, funcao, funcaoDiferencial):
    return (-1)*9.80665 - funcaoDiferencial*abs(funcaoDiferencial)

#metodoEuler(funcaoDiferencial,0,1,0.1,0)
#print(rungeKutta2(funcaoDiferencial, 0, 1, 0.1, 0))
#print(rungeKutta4(funcaoDiferencial, 0, 1, 0.1, 0))
#print(taylor2Ordem(funcaoDiferencial2Ordem, 0, 1, 0.1, 0, 0))
print(rungeKuttaNystron(funcaoDiferencial2Ordem, 0, 1, 0.1, 0, 0))