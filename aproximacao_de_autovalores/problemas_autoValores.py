import numpy as np

# Metodos para aproximação de raizes de polinomios

def MetodoHorne(n = 2, x0 = 0, coef = []):
    if coef == None: return (None, None)
    y = coef[n]
    z = coef[n]
    for j in range(n-1, 0, -1):
        y = x0*y + coef[j]
        z = x0*z + y
    y = x0*y + coef[0]
    return (y, z)

def MetodoNewton(coef, p0 = 1, tol=10**(-3), N=300):
    n = len(coef) - 1
    i = 1
    while i <= N:
        f, f_derivada = MetodoHorne(n, p0, coef)
        print("f: ", f)
        print("f_derivada: ", f_derivada)
        p = p0 - (f / f_derivada)
        if abs(p - p0) < tol:
            return p
        print("aproximaçao p: ", p,"\n")
        i += 1
        p0 = p
    return None

def Busca(x):
    p = 0
    max = x[p]
    for i in range(len(x)):
        if abs(x[i]) > max:
            p = i
    return p

def MetodoDaPotencia(n, matriz, x, tol, max_it):
    k = 1
    x_aux = x.copy()
    p = Busca(x_aux)
    x_aux = x_aux / x_aux[p]
    while k <= max_it:
        y = matriz.dot(x_aux)
        mu = y[p]
        p = Busca(y)
        if y[p] == 0:
            return (0, x_aux)
        ERR = max(x_aux - (y / y[p]))
        x_aux = y / y[p]
        print(x_aux)
        print(mu)
        print("")
        if ERR < tol:
            return (mu, x_aux)
        k = k + 1
    return (None, None)

'''Teste metodo de Newton'''    
# coef = [-4, 3, -3, 0, 2]
# x0 = -2
# n = 4
# x = MetodoNewton(coef, -2, 10**(-5))


'''Teste metodo da potencia'''    
matriz = np.array([(-2,-3),(6,7)])
dim_matriz = len(matriz)
x = np.array([1,1])
tol = 10**(-3)
max_it = 300
auto_valor, auto_vetor = MetodoDaPotencia(dim_matriz, matriz, x, tol, max_it)