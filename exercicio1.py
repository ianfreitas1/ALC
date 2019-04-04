import sys
from math import sqrt


A = [[16, 9, 8, 7, 6, 5, 4, 3, 2, 1],
     [9, 17, 9, 8, 7, 6, 5, 4, 3, 2],
     [8, 9, 18, 9, 8, 7, 6, 5, 4, 3],
     [7, 8, 9, 19, 9, 8, 7, 6, 5, 4],
     [6, 7, 8, 9, 18, 9, 8, 7, 6, 5],
     [5, 6, 7, 8, 9, 17, 9, 8, 7, 6],
     [4, 5, 6, 7, 8, 9, 16, 9, 8, 7],
     [3, 4, 5, 6, 7, 8, 9, 15, 9, 8],
     [2, 3, 4, 5, 6, 7, 8, 9, 14, 9],
     [1, 2, 3, 4, 5, 6, 7, 8, 9, 13],]

B = [[1,2,2],
     [4,4,2],
     [4,6,4]]
C = [[1,0.2,0.4],
     [0.2,1,0.5],
     [0.4,0.5,1]]

def decomposicaoLU(matriz):
    n = len(matriz)
    for k in range(n):
        for i in range(k+1, n):
            matriz[i][k] = float(matriz[i][k]/matriz[k][k])
        for j in range(k+1, n):
            for i in range(k+1, n):
                matriz[i][j] = matriz[i][j]-matriz[i][k]*matriz[k][j]
    return matriz

def decomposicaoCholesky(A):
    n = len(A)
    for i in range(n):
        A[i][i] = sqrt(A[i][i])
        for j in range(i+1, n):
            A[j][i] = float(A[j][i]/A[i][i])
        for k in range(i+1, n):
            for j in range(k, n):
                A[j][k] = A[j][k] - A[j][i]*A[k][i]
    for a in range(n):
        for b in range(n):
            print A[a][b]
    return A

#decomposicaoLU(B)
decomposicaoCholesky(C)