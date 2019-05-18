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

def funcao(x):
    return math.sin(x)

print(derivadaPassoFrente(funcao, 1, 0.5))
print(derivadaPassoTras(funcao, 1, 0.5))
print(derivadaCentral(funcao, 1, 0.5))