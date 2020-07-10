import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix
import time
import sys

def GeraVetor_b(sparse_matrix):
    '''Gera um vetor b somando todos os elementos de cada linha da matriz esparsa'''
    b = np.zeros(sparse_matrix.shape[0])
    for i in range(len(sparse_matrix)):
        b[i] += sum(sparse_matrix[i])
    return b

def main(matriz, b):
    tolerancia = 10 ** (-3)
    max_int = 300
    w = 1.5
    # print("\nMetodo de Jacobi\n")
    # sol_Jacobi, num_interacoes_Jacobi = MetodoJacobi(matriz, b, tolerancia, max_int)
    # print("\nMetodo de Gauss-Seidel\n")
    # sol_GS, num_interacoes_GS = MetodoGaussSiedel(matriz, b, tolerancia, max_int)
    print("\nMetodo de relaxação SOR\n")
    sol_SOR, num_interacoes_SOR = SOR(matriz, b, tolerancia, max_int, w)
        
def MetodoJacobi(matriz, b, tol=10^-3, max_int=300):
    '''Aplica o metodo de Jacobi sobre uma matriz A, para calcular a solucao x do sistema Ax=b'''
    k = 0
    dim_matriz = len(matriz)
    x = np.zeros(dim_matriz)
    x_0_aux = np.zeros(len(matriz))
    while k <= max_int:
        for i in range(dim_matriz):
            sum_aux = 0
            for j in range(dim_matriz):
                if j == i: continue
                if matriz[i][j] == 0: continue
                print(j)
                sum_aux += matriz[i][j]*x_0_aux[j]
            x[i] = (-sum_aux + b[i]) / matriz[i][i]
        k += 1
        # print("k = ", k,"\nx = ",x,"\n")
        if max(np.abs(x - x_0_aux)) / max(abs(x)) < tol: return (x, k)
        for l in range(dim_matriz):
            x_0_aux[l] = x[l]
    return (None, None)

def MetodoGaussSiedel(matriz, b, tol=10^-3, max_int=300):
    '''Aplica o metodo de Gauss-Siedel sobre uma matriz A, para calcular a solucao x do sistema Ax=b'''
    k = 0
    dim_matriz = len(matriz)
    x = np.zeros(dim_matriz)
    x_0_aux = np.zeros(len(matriz))
    while k <= max_int:
        for i in range(dim_matriz):
            sum_aux1 = 0
            sum_aux2 = 0
            for j in range(dim_matriz):
                if j == i: continue
                if matriz[i][j] == 0: continue
                if j <= i-1: sum_aux1 += matriz[i][j]*x[j]
                else: sum_aux2 += matriz[i][j]*x_0_aux[j]
            x[i] = (- sum_aux1 - sum_aux2 + b[i]) / matriz[i][i]
        k += 1
        # print("k = ", k,"\nx = ",x,"\n")
        if max(np.abs(x - x_0_aux)) / max(abs(x)) < tol: return (x, k)
        for l in range(dim_matriz):
            x_0_aux[l] = x[l]
    return (None, None)

def SOR(matriz, b, tol=10^-3, max_int=300, w=1):
    '''Aplica uma refinaçao w, sobre o metodo de Gauss-Siedel sobre uma matriz A, para calcular a solucao x do sistema Ax=b'''
    k = 0    
    dim_matriz = len(matriz)
    x = np.zeros(dim_matriz)
    x_0_aux = np.zeros(dim_matriz)
    while k <= max_int:
        for i in range(dim_matriz):
            sum_aux1 = 0
            sum_aux2 = 0
            for j in range(dim_matriz):
                if j == i: continue
                if matriz[i][j] == 0: continue
                if j <= i-1: sum_aux1 += matriz[i][j]*x[j]
                else: sum_aux2 += matriz[i][j]*x_0_aux[j]
            x[i] = (1-w)*x_0_aux[i] + ((w*(- sum_aux1 - sum_aux2 + b[i])) / matriz[i][i])
        k += 1
        # print("k = ", k,"\nx = ",x,"\n")
        if max(np.abs(x - x_0_aux)) / max(abs(x)) < tol: return (x, k)
        for l in range(dim_matriz):
            x_0_aux[l] = x[l]
    return (x, k)

def GradienteConjugadoNormal(matriz, b, dim_matriz, tol=10**(-3)):
    x = np.zeros(len(matriz))
    r = b - matriz.dot(x)
    if r.dot(r) < tol: return (x,0)
    v = r
    alfa = v.dot(v)
    k = 0
    while k <= dim_matriz:
        u = matriz.dot(v)
        t = alfa / (v.dot(u))
        x = x + t*v
        r = r - t*u
        beta = r.dot(r)
        if beta < tol:
            return (x,k)
        s = beta / alfa
        v = r + s*v
        alfa = beta
        k += 1
    return (x, k)

def SelecionaMatrizTeste(teste):
    if teste == 1:
        # Converge para: Gauss-Siedel e Jacobi
        matriz = np.array([(10,-1,2,0),(-1,11,-1,3),(2,-1,10,-1),(0,3,-1,8)])
        b = np.array([6,25,-11,15])
        return (matriz, b)
    if teste == 2:
        # Converge para: Gauss-Siedel
        matriz = np.array([(2,-1,1),(2,2,2),(-1,-1,2)])
        b = np.array([-1,4,-5])
        return (matriz, b)
    if teste == 3:
        # Converge para: Gauss-Siedel e Jacobi
        matriz = np.array([(4,1,-1),(-1,3,1),(2,2,5)])
        b = np.array([5,-4,1])
        return (matriz, b)
    if teste == 4:
        # Converge para: Gauss-Siedel e Jacobi
        matriz = np.array([(-2,1,1/2),(1,-2,-1/2),(0,1,2)])
        b = np.array([4,-4,0])
        return (matriz, b)
    if teste == 5:
        # Converge para: Gauss-Siedel e Jacobi
        matriz = np.array([(4,1,-1,1),(1,4,-1,-1),(-1,-1,5,1),(1,-1,1,3)])
        b = np.array([-2,-1,0,1])
        return (matriz, b)
    if teste == 6:
        # Converge para: Gauss-Siedel e Jacobi
        matriz = np.array([(4,-1,0,-1,0,0),(-1,4,-1,0,-1,0),(0,-1,4,0,0,-1),(-1,0,0,4,-1,0),(0,-1,0,-1,4,-1),(0,0,-1,0,-1,4)])
        b = np.array([0,5,0,6,-2,6])
        return (matriz, b)
    if teste == 7:
        # Converge para: Gauss-Siedel e Jacobi (Jacobi em quase 300 interações)
        matriz = np.array([(1,0,-1),(-1/2,1,-1/4),(1,-1/2,1)])
        b = np.array([0.2, -1.425, 2])
        return (matriz, b)
    if teste == 8:
        # Converge para: ninguem
        matriz = np.array([(1,0,-2),(-1/2,1,-1/4),(1,-1/2,1)])
        b = np.array([0.2, -1.425, 2])
        return (matriz, b)
    if teste == 9:
        sparse_matrix_coo = csr_matrix(mmread("./LFAT5/LFAT5.mtx").toarray())
        # Gera a forma completa da matriz esparsa
        sparse_matrix = sparse_matrix_coo.toarray()
        # Gera o vetor b com a matriz esparsa
        b = GeraVetor_b(sparse_matrix)
        return (sparse_matrix, b)

if __name__ == "__main__":
    teste = 9
    matriz, b = SelecionaMatrizTeste(teste)
    main(matriz, b)