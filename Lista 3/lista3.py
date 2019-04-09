# Obtem os parametros a e b de uma reta a partir do metodo da regressao linear
def regressaoLinear(x, y):
    n = len(x)
    a11 = float(n)
    a12 = 0.0
    a21 = 0.0
    a22 = 0.0
    for i in range(n):
        a12 += x[i]
        a21 += x[i]
        a22 += x[i]**2
    detA = (a11*a22) - (a12*a21)
    invA = [[a22/detA, -a12/detA],[-a21/detA, a11/detA]]
    c1 = 0
    c2 = 0
    for i in range(n):
        c1 += y[i]
        c2 += x[i]*y[i]
    C = [c1, c2]
    B = multMV(invA, C)
    return B
    
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

#print regressaoLinear([1,2,3],[2,3.5,6.5])