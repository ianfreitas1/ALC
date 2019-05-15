import math
import numpy as np

# Funcao que retorna a propria matriz apos a decomposicao LU
def decomposicaoLU(matriz):
    det = np.linalg.det(matriz)
    if (np.allclose([det], [0], 1e-8)):
        raise ValueError("Erro: Matriz singular")
    n = len(matriz)
    for k in range(n):
        for i in range(k+1, n):
            matriz[i][k] = float(matriz[i][k]/matriz[k][k])
        for j in range(k+1, n):
            for i in range(k+1, n):
                matriz[i][j] = matriz[i][j]-matriz[i][k]*matriz[k][j]
    return matriz

# Funcao que retorna a triangular inferior apos a decomposicao de Cholesky
def Cholesky(A):
    if (np.allclose(A, transposta(A), 1e-8) == 0):
        raise ValueError("Erro: Matriz nao simetrica")
    n = len(A)
    L = [[0.0] * n for i in range(n)]
    for i in range(n):
        for j in range(i+1):
            s = sum(L[i][k]*L[j][k] for k in range(j))
            L[i][j] = math.sqrt(A[i][i]-s) if (i==j) else (1.0/L[j][j]*(A[i][j]-s))
    return L

# Funcao que realiza uma substituicao para frente para decomposicao LU
def subFrente(L, B):
    n = len(L)
    Y = []
    for i in range(n):
        Y.append(0)
    Y[0] = B[0]
    for i in range(1,n):
        s = 0
        for j in range(i):
            s = s + L[i][j]*Y[j]
        Y[i] = B[i] - s
    return Y

# Funcao que realiza uma substituicao para frente para decomposicao de Cholesky
def subFrente2(L,B):
    n = len(L)
    Y = []
    for i in range(n):
        Y.append(0)
    for i in range(n):
        a = B[i]
        for j in range(i):
            a = a - L[i][j]*Y[j]
        Y[i] = float(a)/L[i][i]
    return Y

# Funcao que realiza uma retro-substituicao
def subTras(U, Y):
    n = len(U)
    X = []
    for i in range(n):
        X.append(0)
    X[n-1] = float(Y[n-1])/U[n-1][n-1]
    for i in range(n-2, -1, -1):
        s = Y[i]
        for j in range(i+1, n):
            s = s - U[i][j]*X[j]
        X[i] = float(s)/U[i][i]
    return X

# Funcao que retorna a transposta da matriz passada como parametro
def transposta(matriz):
    n = len(matriz)
    resp = [[0.0]*n for i in range(n)]
    for i in range(n):
        for j in range(n):
            resp[j][i] = matriz[i][j]
    return resp