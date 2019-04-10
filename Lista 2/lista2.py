import math
import numpy as np

# Calculo do autovalor e autovetor correspondente utilizando o Power Method
def powerMethod(A):
    n = len(A)
    X0 = [1.0 for i in range(n)]
    Y = multMV(A, X0)
    y0 = 1
    y1 = Y[0]
    k = 1
    for i in range(n):
        Y[i] = Y[i]/y1
        X0 = Y
    R = math.fabs(y1 - y0)/y1
    while (R>0.001):
        y0 = y1
        Y = multMV(A, X0)
        y1 = Y[0]
        for i in range(n):
            Y[i] = Y[i]/y1
            X0 = Y
        R = math.fabs(y1-y0)/y1
        k+=1
    print("Autovalor: " + str(y1))
    print("Autovetor: " + str(X0))
    print("Iteracoes: " + str(k))

# Calculo dos autovalores e autovetores correspondentes utilizando o metodo de Jacobi    
def jacobiNumerico(A):
    if (np.allclose(A, transposta(A), 1e-8) == 0):
        raise ValueError("Erro: Matriz nao simetrica")
    n = len(A)
    X = [[float(i==j) for j in range(n)] for i in range(n)]
    tol = 10**(-3)
    maxi = maxElem(A)
    while (A[maxi[0]][maxi[1]] > tol):
        P = matrizP(A, maxi)
        transP = transposta(P)
        A = multiplyMatrices(transP, multiplyMatrices(A, P))
        X = multiplyMatrices(X, P)
        maxi = maxElem(A)
    autovalores = []
    for i in range(n):
        autovalores.append(A[i][i])
    print("Autovalores: " + str(autovalores))
    print("Autovetores: " + str(X))
        
        
# Metodo de Jacobi para resolucao de sistemas lineares
def jacobiIterativo(A, B):
    if (diagDominante(A) == False):
        raise ValueError("Erro: Matriz nao e diagonal dominante")
    n = len(A)
    X0 = [1.0 for i in range(n)]
    tol = 10**(-2)
    X1 = [0.0 for i in range(n)]
    R = tol + 1
    i = 0
    numerador = denominador = 0
    while (R > tol):
        for j in range(n):
            c = 0
            for k in range(n):
                if (j!=k):
                    c += (A[j][k] * X0[k])
                X1[j] = (B[j]-c)/A[j][j]
        for z in range(n):
            numerador += (X1[z]-X0[z])**2
            denominador += X1[z]**2
        R = float(math.sqrt(numerador))/math.sqrt(denominador)
        X0 = X1
        i += 1
    print("X1: ", X1)
    print("Iteracoes: " + str(i))
    
    
# Metodo de Gauss-Seidel para resolucao de sistemas lineares
def gaussSeidel(A,B):
    if (diagDominante(A) == False):
        raise ValueError("Erro: Matriz nao e diagonal dominante")
    n = len(A)
    X0 = [1.0 for i in range(n)]
    tol = 10**(-2)
    X1 = [0.0 for i in range(n)]
    R = tol + 1
    i = 0
    numerador = 0
    denominador = 0
    c = 0
    d = 0
    while (R > tol):
        for j in range(n):
            c = mult(A[j][:j], X1[:j])
            d = mult(A[j][j+1:], X0[j+1:])
            X1[j] = (B[j] - c - d)/A[j][j]
        for z in range(n):
            numerador += (X1[z]-X0[z])**2
            denominador += X1[z]**2
        R = float(math.sqrt(numerador))/math.sqrt(denominador)
        X0 = X1
        i += 1
    print("X1: ", X1)
    print("R: ", R)
    print("Iteracoes: " + str(i))


def diagDominante(A):
    n = len(A)
    for i in range(n):
        for j in range(n):
            if (i != j and math.fabs(A[i][i]) >= math.fabs(A[i][j]) and math.fabs(A[i][i]) >= math.fabs(A[j][i])):
                return True
    return False
                
# Funcao auxiliar de Jacobi que retorna o maior elemento fora da diagonal da matriz
def maxElem(A):
    n = len(A)
    s = -100
    for i in range(n):
        for j in range(n):
            if (i != j and A[i][j] > s):
                s = A[i][j]
                index = (i, j)
    return index

# Retorna o valor de phi para o metodo de Jacobi
def phi(A, maxi):
    if(A[maxi[0]][maxi[0]] == A[maxi[1]][maxi[1]]):
        return math.pi/4
    else:
        return math.atan2(2*A[maxi[0]][maxi[1]], (A[maxi[0]][maxi[0]] - A[maxi[1]][maxi[1]]))/2

# Retorna a matriz P utilizada no metodo de Jacobi
def matrizP(A, maxi):
    n = len(A)
    P = [[0.0 for i in range(n)] for j in range(n)]
    for i in range(n):
        P[i][i] = 1
        
    phi1 = phi(A, maxi)
    P[maxi[0]][maxi[0]] = math.cos(phi1)
    P[maxi[1]][maxi[1]] = math.cos(phi1)
    P[maxi[0]][maxi[1]] = -math.sin(phi1)
    P[maxi[1]][maxi[0]] = math.sin(phi1)
    
    return P

# Funcao auxiliar de Gauss-Seidel para multiplicacao de vetores
def mult(A, B):
    s = 0
    for i in range(len(A)):
        for j in range(len(B)):
            if (i==j):
                s += A[i]*B[i]
    return s
    
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

# Multiplica duas matrizes
def multiplyMatrices(X, Y):
    n = len(X)
    result = [[0.0 for j in range(n)] for i in range(n)]
    for i in range(len(X)):
        for j in range(len(Y[0])):
            for k in range(len(Y)):
                result[i][j] += X[i][k] * Y[k][j]
    return result

# Retorna a transposta de uma matriz
def transposta(matriz):
    n = len(matriz)
    resp = [[0.0]*n for i in range(n)]
    for i in range(n):
        for j in range(n):
            resp[j][i] = matriz[i][j]
    return resp
    
#powerMethod([[3,2,0],[2,3,-1],[0,-1,3]])
#gaussSeidel([[3,-1,-1],[-1,3,-1],[-1,-1,3]],[1,2,1])
#jacobiIterativo([[3,-1,-1],[-1,3,-1],[-1,-1,3]],[1,2,1])
#print(jacobiNumerico([[3,0.4,5],[0.4,4,0.1],[5,0.1,-2]])
