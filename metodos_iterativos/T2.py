import numpy as np
import pandas as pd
from scipy.io import mmread
import time
import sys

def GeraVetor_b(sparse_matrix):
    '''Gera um vetor b somando todos os elementos de cada linha da matriz esparsa'''
    b = np.zeros(sparse_matrix.shape[0])
    for i in range(len(sparse_matrix)):
        b[i] += sum(sparse_matrix[i])
    return b

def OrdenaCOO(sparse_matrix):
    rows = []
    cols = []
    for i in range(sparse_matrix.shape[0]):
        for j in range(sparse_matrix.shape[1]):
            if sparse_matrix[i][j] != 0.0:
                rows.append(i)
                cols.append(j)
    return (rows, cols)

# Considerações:
# A matriz no formato COO esta ordenada por linhas (i=0/j=1,2,5,...,n; i=1/j=3,7,8,...,n)

def ProdutoVetorMatrizCOO(matriz_COO, rows, cols, vetor):
    '''Aplica a multiplicação de uma matriz no formato COO por um vetor (A*x)'''
    result = np.zeros(len(vetor))
    l = 0
    for i in range(len(result)):
        while l < len(rows):
            if rows[l] == i:
                result[i] += matriz_COO[i][cols[l]] * vetor[cols[l]]
            if rows[l] > i:
                break
            l += 1
    return result

def main():
    # Le o arquivo e ordena a representação COO da matriz esparsa por linha
    sparse_matrix_coo = mmread("./LFAT5/LFAT5.mtx")
    # Gera a forma completa da matriz esparsa
    sparse_matrix = sparse_matrix_coo.toarray()
    # Gera o vetor b com a matriz esparsa
    b = GeraVetor_b(sparse_matrix)
    # Gera um array de todas as posições "i" e "j" com elementos não nulos
    rows, cols = OrdenaCOO(sparse_matrix)
    x_result_SOR, num_int_SOR = SOR_SparseMatrix(sparse_matrix, b, rows, cols)

# Considerações:
# - A matriz possui uma diagonal principal com todos os elementos não nulos
# - A matriz é sempre quadrada
# - A representação COO da matriz esta ordenada por linha

def SOR_SparseMatrix(matriz, b, rows, cols, w=1, tol=10**(-3), max_int=300):
    '''Aplica uma refinaçao w, sobre o metodo de Gauss-Siedel, para calcular a solucao x do sistema Ax=b'''
    k = 0    
    dim_matriz = len(matriz)
    x = np.zeros(dim_matriz)
    x_0_aux = np.zeros(len(matriz))
    while k < max_int:
        i_index = 0
        for i in range(dim_matriz):
            sum_aux1 = 0
            sum_aux2 = 0
            while i_index < len(rows):
                if rows[i_index] == i and cols[i_index] == i: 
                    i_index += 1
                    continue
                if rows[i_index] == i and cols[i_index] <= i-1: 
                    sum_aux1 += matriz[i][cols[i_index]]*x[cols[i_index]]
                    i_index += 1
                    continue
                if rows[i_index] == i and cols[i_index] > i: 
                    sum_aux2 += matriz[i][cols[i_index]]*x_0_aux[cols[i_index]]
                    i_index += 1
                    continue
                else: break
            x[i] = (1-w)*x_0_aux[i] + ((w*(- sum_aux1 - sum_aux2 + b[i])) / matriz[i][i])
        k += 1
#         print("k = ", k,"\nx = ",x,"\n")
        if max(np.abs(x - x_0_aux)) / max(abs(x)) < tol: return (x, k)
        for l in range(dim_matriz):
            x_0_aux[l] = x[l]
    return (x, k)

def GradienteConjugadoNormalCOO(matriz, rows, cols, b, dim_matriz, tol=10**(-3)):
    x = np.zeros(len(matriz))
    r = b - ProdutoVetorMatrizCOO(matriz, rows, cols, x)
    if r.dot(r) < tol: return (x,0)
    v = r
    alfa = v.dot(v)
    k = 0
    while k <= dim_matriz:
        u = ProdutoVetorMatrizCOO(matriz, rows, cols, v)
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

if __name__ == "__main__":
    main()