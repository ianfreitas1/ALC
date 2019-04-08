from math import sqrt

def gaussSeidel(A,B):
    n = len(A)
    X0 = [1.0 for i in range(n)]
    tol = 10**(-5)
    X1 = [0.0 for i in range(n)]
    R = tol + 1
    i = 0
    numerador = 0
    denominador = 0
    while (R>tol):
        c = 0
        d = 0   
        for y in range(n):
            for j in range(0, y-1):
                print(y,j)
                c += A[y][j] * X1[j]
            for k in range(y+1, n):
                #print(y, k)
                d += A[y][k] * X0[k]
            X1[y] = float(B[y] - c - d)/A[y][y]
        for z in range(n):
            numerador += (X1[z]-X0[z])**2
            denominador += X1[z]**2
        R = float(sqrt(numerador))/sqrt(denominador)
        X0 = X1
        i += 1
            
    print(X1)
    print(R)
gaussSeidel([[3,-1,-1],[-1,3,-1],[-1,-1,3]],[1,2,1])