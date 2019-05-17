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

def funcaoDiferencial(t, funcao):
    return t + funcao

#metodoEuler(funcaoDiferencial,0,1,0.1,0)
#print(rungeKutta2(funcaoDiferencial, 0, 1, 0.1, 0))
print(rungeKutta4(funcaoDiferencial, 0, 1, 0.1, 0))