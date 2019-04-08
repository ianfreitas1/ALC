from lista3 import *

if __name__ == "__main__":
    print("Executando metodo da Regressao Linear")
    num = int(input("Escolha o numero de pontos: "))
    X = []
    Y = []
    for i in range(num):
        X.append(float(input(("X"+str(i)+": "))))
        Y.append(float(input(("Y"+str(i)+": "))))
    print("Coeficientes: "+ str(regressaoLinear(X, Y)[0]) + ", " + str(regressaoLinear(X, Y)[1]))
    