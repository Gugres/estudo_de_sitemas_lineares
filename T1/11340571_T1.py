import numpy as np

def GeraMatrizHilbert(n):
    matrizHilbert = np.zeros((n,n+1))
    for i in range(n):
        for j in range(n + 1):
            if j != n:
                matrizHilbert[i][j] = 1/((i+1) + (j+1) - 1)
            else:
                matrizHilbert[i][j] = np.sum(matrizHilbert[i]) 
    return matrizHilbert

def main():
    # dim_matriz = input(int("Entre com a dimensão (n) da matriz A, do sistema Ax=b: "))
    dim_matriz = 2
    matrizHilbert = GeraMatrizHilbert(dim_matriz)
    # solucao, determinante = EliminacaoSemPivotamento(matrizHilbert, dim_matriz, dim_matriz + 1)
    EliminacaoSemPivotamento(matrizHilbert, dim_matriz, dim_matriz + 1)

def EliminacaoSemPivotamento(matriz, i, j):
    # print(matriz)
    trocasDeLinha = 0
    p = 0
    linha = 0
    x = [0] * i
    while linha < (i - 1):
        if matriz[p][linha] == 0:
            p += 1
            if p == i: return None
            continue
        if p != linha:
            matriz[(linha,p),:] = matriz[(p,linha),:]
        linha_final = linha + 1
        while linha_final < i:
            m = matriz[linha_final][linha] / matriz[linha][linha]
            matriz[linha_final] = matriz[linha_final] - m*matriz[linha]
            linha_final += 1
        linha += 1
        p = linha
    # x = CalculaResultado(matriz, x, i-1)
    determinante = CalculaDeterminante(matriz, i, trocasDeLinha)
    print(determinante)

def EliminacaoComPivotamento(matriz, i, j):
    pass

def CalculaDeterminante(matriz, n, trocasDeLinha):
    determinante = 1
    while 0 < n:
        determinante *= matriz[n-1][n-1]
        n -= 1
    determinante *= pow(-1,trocasDeLinha)
    return determinante

def CalculaResultado(matriz, x, n):
    if matriz[n][n] == 0: return None
    x[n] = matriz[n][n+1] / matriz[n][n]
    n -= 1
    while 0 <= n:
        bakward = 0


if __name__ == "__main__":
    main()