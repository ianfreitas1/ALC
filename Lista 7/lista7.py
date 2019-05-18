import math

def derivadaPassoFrente(funcao, x, deltaX):
    top = funcao(x + deltaX) - funcao(x)
    bottom = deltaX
    return top/bottom

def derivadaPassoTras(funcao, x, deltaX):
    top = funcao(x) - funcao(x - deltaX)
    bottom = deltaX
    return top/bottom

def derivadaCentral(funcao, x, deltaX):
    top = funcao(x + deltaX) - funcao(x - deltaX)
    bottom = 2*deltaX
    return top/bottom

def derivadaPassoFrenteRichard(funcao, x, deltaX, p):
    d1 = derivadaPassoFrente(funcao, x, deltaX)
    deltaX2 = deltaX/2
    d2 = derivadaPassoFrente(funcao, x, deltaX2)
    q = deltaX/deltaX2
    resultado = d1 + (d1 - d2)/(q**(-p)-1)
    return resultado

def derivadaPassoTrasRichard(funcao, x, deltaX, p):
    d1 = derivadaPassoTras(funcao, x, deltaX)
    deltaX2 = deltaX/2
    d2 = derivadaPassoTras(funcao, x, deltaX2)
    q = deltaX/deltaX2
    resultado = d1 + (d1 - d2)/(q**(-p)-1)
    return resultado

def derivadaCentralRichard(funcao, x, deltaX, p):
    d1 = derivadaCentral(funcao, x, deltaX)
    deltaX2 = deltaX/2
    d2 = derivadaCentral(funcao, x, deltaX2)
    q = deltaX/deltaX2
    resultado = d1 + (d1 - d2)/(q**(-p)-1)
    return resultado

def funcao(x):
    return math.sin(x)

print(derivadaPassoFrenteRichard(funcao, 1, 0.5, 1))
print(derivadaPassoTrasRichard(funcao, 1, 0.5, 1))
print(derivadaCentralRichard(funcao, 1, 0.5, 1))