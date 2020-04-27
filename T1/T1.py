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

def main(dim_matriz):
    # dim_matriz = input(int("Entre com a dimensão (n) da matriz A, do sistema Ax=b: "))
    # dim_matriz = 50
    matrizHilbert = GeraMatrizHilbert(dim_matriz)
    # print("Sistema Ax = b dado por: \nA = ", matrizHilbert[:,:(dim_matriz)].tolist(), 
    #                 "\nb = ", matrizHilbert[:,dim_matriz:].tolist())
    resultado_semPivo, determinante_semPivo = EliminacaoSemPivotamento(matrizHilbert, dim_matriz, dim_matriz + 1)
    resultado_comPivo, determinante_comPivo = EliminacaoComPivotamento(matrizHilbert, dim_matriz, dim_matriz + 1)
    print("Norma^2 sem pivotamento: ", CalculaNormaDaDiferenca(resultado_semPivo, dim_matriz),
            "\nNorma^2 com pivotamento: ", CalculaNormaDaDiferenca(resultado_comPivo, dim_matriz))
    print("Determinante sem pivotamento: ", determinante_semPivo,"\nDeterminante com pivotamento: ", determinante_comPivo)

def EliminacaoSemPivotamento(matriz, i, j):
    matriz_aux = matriz.copy()
    trocasDeLinha = 0
    p = 0
    linha = 0
    x = [0] * i
    while linha < (i-1):
        if matriz_aux[p][linha] == 0:
            p += 1
            if p == i: return None
            continue
        if p != linha:
            matriz_aux[(linha,p),:] = matriz_aux[(p,linha),:]
            trocasDeLinha += 1
        linha_final = linha + 1
        while linha_final < i:
            m = matriz_aux[linha_final][linha] / matriz_aux[linha][linha]
            matriz_aux[linha_final] = matriz_aux[linha_final] - m*matriz_aux[linha]
            linha_final += 1
        linha += 1
        p = linha
    x = CalculaResultado(matriz_aux, x, i-1)
    determinante = CalculaDeterminante(matriz_aux, i, trocasDeLinha)
    return (x, determinante)

def EliminacaoComPivotamento(matriz, i, j):
    matriz_aux = matriz.copy()
    trocasDeLinha = 0
    linha = 0
    x = [0] * i
    while linha < (i-1):
        linha_maior = linha
        p = linha + 1
        max = matriz_aux[linha_maior][linha]
        while p < i:
            if max < matriz_aux[p][linha]: 
                max = matriz_aux[p][linha]
                linha_maior = p
            p += 1 
        if linha_maior != linha:
            matriz_aux[(linha,linha_maior),:] = matriz_aux[(linha_maior,linha),:]
            trocasDeLinha += 1
        linha_final = linha + 1
        while linha_final < i:
            m = matriz_aux[linha_final][linha] / matriz_aux[linha][linha]
            matriz_aux[linha_final] = matriz_aux[linha_final] - m*matriz_aux[linha]
            linha_final += 1
        linha += 1
    x = CalculaResultado(matriz_aux, x, i-1)
    determinante = CalculaDeterminante(matriz_aux, i, trocasDeLinha)
    return (x, determinante)

def CalculaDeterminante(matriz, n, trocasDeLinha):
    determinante = 1
    while 0 < n:
        determinante *= matriz[n-1][n-1]
        n -= 1
    determinante *= pow(-1,trocasDeLinha)
    return determinante

def CalculaResultado(matriz, x, n):
    linha = n
    coluna = n + 1
    if matriz[n][n] == 0: return None
    x[n] = matriz[linha][coluna] / matriz[linha][linha]
    linha -= 1
    while 0 <= linha:
        coef_aux = 0
        linha_aux = linha + 1
        while linha_aux <= n:
            coef_aux += matriz[linha][linha_aux]*x[linha_aux]
            linha_aux += 1
        x[linha] = (matriz[linha][coluna] - coef_aux) / matriz[linha][linha]
        linha -= 1
    return x

def CalculaNormaDaDiferenca(resultado, n):
    resultado_correto = np.array([1] * n)
    resultado = np.array(resultado)
    norma_2 = sum(resultado - resultado_correto)
    return norma_2

def TestesMocados():
    matriz = np.array([(1.0,2.0,3.0,6.0),(2.0,3.0,5.0,10.0),(1.0,3.0,1.0,5.0)])    
    resultado_semPivo, determinante_semPivo = EliminacaoSemPivotamento(matriz, 3, 3 + 1)
    resultado_comPivo, determinante_comPivo = EliminacaoComPivotamento(matriz, 3, 3 + 1)
    print("Norma^2 sem pivotamento: ", CalculaNormaDaDiferenca(resultado_semPivo, 3),
            "\nNorma^2 com pivotamento: ", CalculaNormaDaDiferenca(resultado_comPivo, 3))
    print("Determinante sem pivotamento: ", determinante_semPivo,"\nDeterminante com pivotamento: ", determinante_comPivo)

if __name__ == "__main__":
    main(7)
    # TestesMocados()