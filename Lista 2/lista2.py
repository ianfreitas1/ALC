from math import fabs, sqrt

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

# Metodo de Jacobi para resolucao de sistemas lineares
def jacobiIterativo(A, B):
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
        R = float(sqrt(numerador))/sqrt(denominador)
        X0 = X1
        i += 1
    print X1
    print("Iteracoes: " + str(i))
    
    
# Metodo de Gauss-Seidel para resolucao de sistemas lineares
def gaussSeidel(A,B):
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
        R = float(sqrt(numerador))/sqrt(denominador)
        X0 = X1
        i += 1
    print(X1)
    print(R)
    print("Iteracoes: " + str(i))

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
    
#powerMethod([[3,2,0],[2,3,-1],[0,-1,3]])
#gaussSeidel([[3,-1,-1],[-1,3,-1],[-1,-1,3]],[1,2,1])
#jacobiIterativo([[3,-1,-1],[-1,3,-1],[-1,-1,3]],[1,2,1])

