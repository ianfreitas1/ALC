from lista2 import *

if __name__ == "__main__":
    print("1: Calcular autovalores e autovetores com Power Method")
    print("2: Calcular autovalores e autovetores com Jacobi numerico")
    print("3: Resolver AX=B utilizando Gauss-Seidel")
    print("4: Resolver AX=B utilizando Jacobi Iterativo")
    resposta = int(input("Escolha um numero: "))
    if resposta == 1:
        tam = int(input("Escolha o tamanho da matriz quadrada: "))
        matriz = [[0.0]*tam for i in range(tam)]
        for i in range(tam):
            for j in range(tam):
                a = float(input("A" + str(i) + str(j) + ": "))
                matriz[i][j] = a
        powerMethod(matriz)
    elif resposta == 2:
        tam = int(input("Escolha o tamanho da matriz quadrada: "))
        matriz = [[0.0]*tam for i in range(tam)]
        for i in range(tam):
            for j in range(tam):
                a = float(input("A" + str(i) + str(j) + ": "))
                matriz[i][j] = a
        jacobiNumerico(matriz)
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
        gaussSeidel(matriz, B)
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
        jacobiIterativo(matriz, B)