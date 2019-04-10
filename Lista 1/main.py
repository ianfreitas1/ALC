from lista1 import *

if __name__ == "__main__":
    print("1: Decomposicao LU")
    print("2: Decomposicao Cholesky")
    print("3: Resolver AX = B utilizando LU")
    print("4: Resolver AX = B utilizando Cholesky")
    print("5: Resolver Exercicio 3 da lista 1")
    resposta = int(input("Escolha um numero: "))
    if resposta == 1:
        tam = int(input("Escolha o tamanho da matriz quadrada: "))
        matriz = [[0.0]*tam for i in range(tam)]
        for i in range(tam):
            for j in range(tam):
                a = float(input("A" + str(i) + str(j) + ": "))
                matriz[i][j] = a
        print("LU: " + str(decomposicaoLU(matriz)))
    elif resposta == 2:
        tam = int(input("Escolha o tamanho da matriz quadrada: "))
        matriz = [[0.0]*tam for i in range(tam)]
        for i in range(tam):
            for j in range(tam):
                a = float(input("A" + str(i) + str(j) + ": "))
                matriz[i][j] = a
        print("Triangular Inferior (L): " + str(Cholesky(matriz)))
        print("Triangular Superior (U): " + str(transposta(Cholesky(matriz))))
    elif resposta == 3:
        tam = int(input("Escolha o tamanho da matriz quadrada: "))
        matriz = [[0.0]*tam for i in range(tam)]
        for i in range(tam):
            for j in range(tam):
                a = float(input("A" + str(i) + str(j) + ": "))
                matriz[i][j] = a
        B = []
        for i in range(tam):
            b = float(input("B" + str(i) + ": "))
            B.append(b)
        A = decomposicaoLU(matriz)
        print("Y: " + str(subFrente(A, B)))
        print("X: " + str(subTras(A, subFrente(A, B))))
    elif resposta == 4:
        tam = int(input("Escolha o tamanho da matriz quadrada: "))
        matriz = [[0.0]*tam for i in range(tam)]
        for i in range(tam):
            for j in range(tam):
                a = float(input("A" + str(i) + str(j) + ": "))
                matriz[i][j] = a
        B = []
        for i in range(tam):
            b = float(input("B" + str(i) + ": "))
            B.append(b)
        L = Cholesky(matriz)
        Y = subFrente2(L, B)
        X = subTras(transposta(L), Y)
        print("L = " + str(L))
        print("Y = " + str(Y))
        print("X = " + str(X))
    elif resposta == 5:
        A = [[16, 9, 8, 7, 6, 5, 4, 3, 2, 1],
             [9, 17, 9, 8, 7, 6, 5, 4, 3, 2],
             [8, 9, 18, 9, 8, 7, 6, 5, 4, 3],
             [7, 8, 9, 19, 9, 8, 7, 6, 5, 4],
             [6, 7, 8, 9, 18, 9, 8, 7, 6, 5],
             [5, 6, 7, 8, 9, 17, 9, 8, 7, 6],
             [4, 5, 6, 7, 8, 9, 16, 9, 8, 7],
             [3, 4, 5, 6, 7, 8, 9, 15, 9, 8],
             [2, 3, 4, 5, 6, 7, 8, 9, 14, 9],
             [1, 2, 3, 4, 5, 6, 7, 8, 9, 13]]
        B = [4.0, 0.0, 8.0, 0.0, 12.0, 0.0, 8.0, 0.0, 4.0, 0.0]
        print("Resolvendo AX=B utilizando decomposicao LU")
        A2 = decomposicaoLU(A)
        print("Y: " + str(subFrente(A2, B)))
        print("X: " + str(subTras(A2, subFrente(A2, B))))
        print("----------------------------------------------")
        print("Resolvendo AX=B utilizando decomposicao de Cholesky")
        A = [[16, 9, 8, 7, 6, 5, 4, 3, 2, 1],
             [9, 17, 9, 8, 7, 6, 5, 4, 3, 2],
             [8, 9, 18, 9, 8, 7, 6, 5, 4, 3],
             [7, 8, 9, 19, 9, 8, 7, 6, 5, 4],
             [6, 7, 8, 9, 18, 9, 8, 7, 6, 5],
             [5, 6, 7, 8, 9, 17, 9, 8, 7, 6],
             [4, 5, 6, 7, 8, 9, 16, 9, 8, 7],
             [3, 4, 5, 6, 7, 8, 9, 15, 9, 8],
             [2, 3, 4, 5, 6, 7, 8, 9, 14, 9],
             [1, 2, 3, 4, 5, 6, 7, 8, 9, 13]]
        B = [4.0, 0.0, 8.0, 0.0, 12.0, 0.0, 8.0, 0.0, 4.0, 0.0]
        A3 = Cholesky(A)
        Y = subFrente2(A3, B)
        X = subTras(transposta(A3), Y)
        print("Y = " + str(Y))
        print("X = " + str(X))