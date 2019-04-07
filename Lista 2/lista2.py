from math import fabs

# Calculo do autovalor e autovetor correspondente utilizando o Power Method
def powerMethod(A):
    n = len(A)
    X0 = [1.0 for i in range(n)]
    Y = multMV(A, X0)
    X1 = [0.0 for i in range(n)]
    y0 = 1
    y1 = Y[0]
    k = 1
    for i in range(n):
        Y[i] = Y[i]/y1
        X0 = Y
    R = fabs(y1 - y0)/y1
    while (R>0.001):
        y0 = y1
        Y = multMV(A, X0)
        y1 = Y[0]
        for i in range(n):
            Y[i] = Y[i]/y1
            X0 = Y
        R = fabs(y1-y0)/y1
        k+=1
    print "Autovalor: " + str(y1)
    print "Autovetor: " + str(X0)
    print "Iteracoes: " + str(k)

# Retorno de multiplicacao matriz x vetor. Funcao auxiliar do Power Method
def multMV(m, v):
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
    
powerMethod([[1,0.2,0],[0.2,1,0.5],[0,0.5,1]])
