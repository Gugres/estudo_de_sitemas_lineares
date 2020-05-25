import numpy as np
import pandas as pd
import time
import sys


def main(matriz, b):
    tolerancia = 10 ** (-3)
    max_int = 300
    w = 0.8
    # print("\nMetodo de Jacobi\n")
    # sol_Jacobi, num_interacoes_Jacobi = MetodoJacobi(matriz, b, tolerancia, max_int)
    print("\nMetodo de Gauss-Seidel\n")
    sol_GS, num_interacoes_GS = MetodoGaussSiedel(matriz, b, tolerancia, max_int)
    print("\nMetodo de relaxação SOR\n")
    sol_SOR, num_interacoes_SOR = SOR(matriz, b, tolerancia, max_int, w)

def MetodoJacobi(matriz, b, tol=10^-3, max_int=300):
    '''Aplica o metodo de Jacobi sobre uma matriz A, para calcular a solucao x do sistema Ax=b'''
    dim_matriz = len(matriz)
    k = 0
    x_0_aux = np.zeros(len(matriz))
    x = np.zeros(dim_matriz)
    while k <= max_int:
        for i in range(dim_matriz):
            sum_aux = 0
            for j in range(dim_matriz):
                if j == i: continue
                if matriz[i][j] == 0: continue
                sum_aux += matriz[i][j]*x_0_aux[j]
            x[i] = (-sum_aux + b[i]) / matriz[i][i]
        k += 1
        print("k = ", k,"\nx = ",x,"\n")
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
        print("k = ", k,"\nx = ",x,"\n")
        if max(np.abs(x - x_0_aux)) / max(abs(x)) < tol: return (x, k)
        for l in range(dim_matriz):
            x_0_aux[l] = x[l]
    return (None, None)

def SOR(matriz, b, tol=10^-3, max_int=300, w=1):
    '''Aplica uma refinaçao w, sobre o metodo de Gauss-Siedel sobre uma matriz A, para calcular a solucao x do sistema Ax=b'''
    dim_matriz = len(matriz)
    k = 0
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
            x[i] = (1-w)*x_0_aux[i] + ((w*(- sum_aux1 - sum_aux2 + b[i])) / matriz[i][i])
        k += 1
        print("k = ", k,"\nx = ",x,"\n")
        if max(np.abs(x - x_0_aux)) / max(abs(x)) < tol: return (x, k)
        for l in range(dim_matriz):
            x_0_aux[l] = x[l]
    return (None, None)

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

if __name__ == "__main__":
    teste = 2
    matriz, b = SelecionaMatrizTeste(teste)
    main(matriz, b)