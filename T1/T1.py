import numpy as np
import pandas as pd

def GeraMatrizHilbert(n):
    '''Gera uma matriz de Hilbert, onde cada termo a_ij = 1 / (i + j - 1)'''
    matrizHilbert = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if j != n:
                matrizHilbert[i][j] = 1/((i+1) + (j+1) - 1)
    return matrizHilbert

def GeraMatrizAleatoria(n):
    matriz = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            matriz[i][j] = np.random.randint(1,100,1)
    return matriz

def GeraVetor(matriz):
    '''Gera um vetor b considerando que b_i = somatoria da linha i da matriz'''
    b = np.array([0.0] * len(matriz))
    for i in range(len(matriz)):
        b[i] = np.sum(matriz[i])
    return b

def main():
    parte = 2
    if parte == 1:
        Parte1()
    if parte == 2:
        Parte2(False)

def Parte1():
    limite_dim = 100
    resultados = pd.DataFrame({'dim_matriz':([0]*(limite_dim-2)),'norma2 sem Pivo':([0]*(limite_dim-2)),'determinante sem Pivo':([0]*(limite_dim-2)),
                                'norma2 com Pivo':([0]*(limite_dim-2)),'determinante com Pivo':([0]*(limite_dim-2)), 'norma2 solucao Numpy': ([0]*(limite_dim-2)), 'determinante solucao Numpy': ([0]*(limite_dim-2))})
    norma2_semPivo = []
    norma2_comPivo = []
    norma2_solucaoNumpy = []
    determinantes_semPivo = []
    determinantes_comPivo = []
    determinantes_solucaoNumpy = []
    dim_matriz = []
    for i in range(2,100):    
        matriz = GeraMatrizHilbert(i)
        b = GeraVetor(matriz)
        matriz_semPivo, resultado_semPivo, determinante_semPivo = EliminacaoSemPivotamento(matriz, b, i)
        matriz_comPivo, resultado_comPivo, determinante_comPivo = EliminacaoComPivotamento(matriz, b, i)
        norma2_semPivo.append(CalculaNormaDaDiferenca(resultado_semPivo, i))
        norma2_comPivo.append(CalculaNormaDaDiferenca(resultado_comPivo, i))
        determinantes_semPivo.append(determinante_semPivo)
        determinantes_comPivo.append(determinante_comPivo)
        dim_matriz.append(i)
        norma2_solucaoNumpy.append(CalculaNormaDaDiferenca(np.linalg.solve(matriz,b), i))
        determinantes_solucaoNumpy.append(np.linalg.det(matriz))
        # if i >= 10:
        #     print(np.linalg.solve(matriz,b),"\n")
        # print("Matriz de Hilbert:\n\n", matriz, "\n")
        # print("Vetor 'b' do sistema: ", b, "\n")
        # print("Matriz escalonada pelo metodo de solução sem pivotamento:\n\n", matriz_semPivo, "\n")
        # print("Matriz escalonada pelo metodo de solução com pivotamento:\n\n", matriz_comPivo, "\n")
        # print("Solucao encontrada pelo metodo de solucao sem pivotamento: ", resultado_semPivo, "\n")
        # print("Solucao encontrada pelo metodo de solucao com pivotamento: ", resultado_comPivo, "\n")
        # return
    resultados['dim_matriz'] = dim_matriz
    resultados['norma2 sem Pivo'] = norma2_semPivo
    resultados['determinante sem Pivo'] = determinantes_semPivo
    resultados['norma2 com Pivo'] = norma2_comPivo
    resultados['determinante com Pivo'] = determinantes_comPivo
    resultados['norma2 solucao Numpy'] = norma2_solucaoNumpy
    resultados['determinante solucao Numpy'] = determinantes_solucaoNumpy
    resultados.to_excel("resultados_matrizHilbert.xlsx", index=False)

def Parte2(teste):
    if teste == False:
        limite_dim = 71
        resultados = pd.DataFrame({'dim_matriz':([0]*(limite_dim-2)),'norma2 Cholesky':([0]*(limite_dim-2)),'determinante Cholesky':([0]*(limite_dim-2)),
                                    'norma2 sem Pivo':([0]*(limite_dim-2)),'determinante sem Pivo':([0]*(limite_dim-2)), 'norma2 solucao Numpy': ([0]*(limite_dim-2)), 'determinante solucao Numpy': ([0]*(limite_dim-2))})
        norma2_Cholesky = []
        norma2_semPivo = []
        norma2_solucaoNumpy = []
        determinantes_Cholesky = []
        determinantes_semPivo = []
        determinantes_solucaoNumpy = []
        dim_matriz = []
        i = 2
        while i < 71:
            matriz = GeraMatrizAleatoria(i)
            b = GeraVetor(matriz)
            matriz_semPivo, resultado_semPivo, determinante_semPivo = EliminacaoSemPivotamento(matriz, b, i)
            if determinante_semPivo == 0:
                continue 
            matriz_simetrica = np.dot(np.transpose(matriz), matriz)
            matrizL, matrizL_trans, resultado_cholesky, determinante_cholesky = DecomposicaoDeCholesky(matriz_simetrica, GeraVetor(matriz_simetrica), i)
            matriz_semPivo, resultado_semPivo, determinante_semPivo = EliminacaoSemPivotamento(matriz_simetrica, b, i)
            norma2_Cholesky.append(CalculaNormaDaDiferenca(resultado_cholesky, i))
            norma2_semPivo.append(CalculaNormaDaDiferenca(resultado_semPivo, i))
            determinantes_Cholesky.append(determinante_cholesky)
            determinantes_semPivo.append(determinante_semPivo)
            dim_matriz.append(i)
            norma2_solucaoNumpy.append(CalculaNormaDaDiferenca(np.linalg.solve(matriz_simetrica,b), i))
            determinantes_solucaoNumpy.append(np.linalg.det(matriz_simetrica))
            # print(matriz)
            # print(np.transpose(matriz))
            # return
            i += 1
            # print("\nMatrizes triangulares Cholesky:\n")
            # print(matrizL, "\n")
            # print(matrizL_trans, "\n")
        resultados['dim_matriz'] = dim_matriz
        resultados['norma2 Cholesky'] = norma2_Cholesky
        resultados['determinante Cholesky'] = determinantes_Cholesky
        resultados['norma2 sem Pivo'] = norma2_semPivo
        resultados['determinante sem Pivo'] = determinantes_semPivo
        resultados['norma2 solucao Numpy'] = norma2_solucaoNumpy
        resultados['determinante solucao Numpy'] = determinantes_solucaoNumpy
        resultados.to_excel("resultados_Cholesky.xlsx", index=False)
    
    else:
        matriz_teste = np.array([(4,2,2),(2,6,2),(2,2,5)])
        b_teste = np.array([8,10,9])
        matrizL, matrizL_trans, resultado_cholesky, determinante_cholesky = DecomposicaoDeCholesky(matriz_teste, b_teste, len(matriz_teste))
        print(matrizL,"\n")
        print(matrizL_trans,"\n\n")
        print(resultado_cholesky)

def EliminacaoSemPivotamento(matriz, b, i):
    '''Recebe uma matriz e um vetor 'b' e realiza a eliminação de Gauss sem pivotamento.
        Retorna o vetor solução 'x' e o determinante da matriz'''
    matriz_aux = np.array(matriz.copy())
    b_aux = b.copy()
    trocasDeLinha = 0
    p = 0
    linha = 0
    try:
        while linha < (i-1):
            if matriz_aux[p][linha] == 0:
                if p == i: return (None, None, None)
                p += 1
                continue
            if p != linha:
                matriz_aux[(linha,p),:] = matriz_aux[(p,linha),:]
                b_aux[[linha,p]] = b_aux[[p,linha]]
                trocasDeLinha += 1
            linha_final = linha + 1
            while linha_final < i:
                m = matriz_aux[linha_final][linha] / matriz_aux[linha][linha]
                matriz_aux[linha_final] = matriz_aux[linha_final] - m*matriz_aux[linha]
                b_aux[linha_final] = b_aux[linha_final] - m*b_aux[linha]
                linha_final += 1
            linha += 1
            p = linha
    except:
        print("Erro no metodo sem pivotamento, na matriz de ordem: ", i, "\n")
        return (None, None, None)
    x = CalculaResultadoMatrizSuperior(matriz_aux, b_aux, i-1)
    determinante = CalculaDeterminante(matriz_aux, i-1, trocasDeLinha)
    return (matriz_aux, x, determinante)

def EliminacaoComPivotamento(matriz, b, i):
    '''Recebe uma matriz e um vetor 'b' e realiza a eliminação de Gauss com pivotamento.
    Retorna o vetor solução 'x' e o determinante da matriz'''
    matriz_aux = matriz.copy()
    b_aux = b.copy()
    trocasDeLinha = 0
    linha = 0
    try:
        while linha < (i-1):
            linha_maior = linha
            p = linha + 1
            max = matriz_aux[linha_maior][linha]
            while p < i:
                if max < matriz_aux[p][linha] and matriz_aux[p][linha] != 0.0: 
                    max = matriz_aux[p][linha]
                    linha_maior = p
                p += 1
            if matriz_aux[linha_maior][linha] == 0: return (None, None, None) 
            if linha_maior != linha:
                matriz_aux[(linha,linha_maior),:] = matriz_aux[(linha_maior,linha),:]
                b_aux[[linha,linha_maior]] = b_aux[[linha_maior,linha]]
                trocasDeLinha += 1
            linha_final = linha + 1
            while linha_final < i:
                m = matriz_aux[linha_final][linha] / matriz_aux[linha][linha]
                matriz_aux[linha_final] = matriz_aux[linha_final] - m*matriz_aux[linha]
                b_aux[linha_final] = b_aux[linha_final] - m*b_aux[linha]
                linha_final += 1
            linha += 1
    except:
        print("Erro no metodo com pivotamento, na matriz de ordem: ", i, "\n")
        return (None, None, None)
    x = CalculaResultadoMatrizSuperior(matriz_aux, b_aux, i-1)
    determinante = CalculaDeterminante(matriz_aux, i-1, trocasDeLinha)
    return (matriz_aux, x, determinante)

def DecomposicaoDeCholesky(matriz, b, n):
    matriz_aux = matriz.copy()
    matrizL = np.zeros((n,n))
    matrizL_trans = np.zeros((n,n))
    try:
        matrizL[0][0] = matriz_aux[0][0] ** (1/2)
        matrizL_trans[0][0] = matriz_aux[0][0] ** (1/2)
        for j in range(1,n):
            matrizL[j][0] = matriz_aux[j][0] / matrizL[0][0]
            matrizL_trans[0][j] = matriz_aux[j][0] / matrizL[0][0]
        for i in range(1, n-1):
            aux = 0
            k = 0 
            while k <= i-1:
                aux += ((matrizL[i][k]) ** 2)
                k += 1
            matrizL[i][i] = (matriz_aux[i][i] - aux) ** (1/2)
            matrizL_trans[i][i] = (matriz_aux[i][i] - aux) ** (1/2)
            for j in range(i+1, n):
                aux2 = 0
                l = 0
                while l <= i-1:
                    aux2 += (matrizL[j][l] * matrizL[i][l])
                    l += 1
                matrizL[j][i] = (matriz_aux[j][i] - aux2) / matrizL[i][i]
                matrizL_trans[i][j] = (matriz_aux[j][i] - aux2) / matrizL[i][i]
        aux3 = 0
        for q in range(0, n-1):
            aux3 += ((matrizL[n-1][q]) ** 2)
        matrizL[n-1][n-1] = (matriz_aux[n-1][n-1] - aux3) ** (1/2)
        matrizL_trans[n-1][n-1] = (matriz_aux[n-1][n-1] - aux3) ** (1/2)
    except:
        print("Erro no metodo de Cholesky, na matriz de ordem: ", n, "\n")
        return (None, None, None, None)
    y = CalculaResultadoMatrizInferior(matrizL, b, n-1)
    x = CalculaResultadoMatrizSuperior(matrizL_trans, y, n-1)
    determinante = (CalculaDeterminante(matrizL, n-1, 0)) * (CalculaDeterminante(matrizL_trans, n-1, 0))
    return (matrizL, matrizL_trans, x, determinante)

def CalculaDeterminante(matriz, n, trocasDeLinha):
    '''Recebe uma matriz e retorna o determinante. Como a matriz passada é uma triângular superior
        calcula o determinante como o produto da diagonal principal * (-1)^numero de trocas de linha'''
    determinante = 1
    while 0 <= n:
        determinante *= matriz[n][n]
        n -= 1
    determinante *= pow(-1,trocasDeLinha)
    return determinante

def CalculaResultadoMatrizSuperior(matriz, b, n):
    '''Recebe uma matriz triangular superior e um vetor 'b', e retorna o vetor 'x' do sistema Ax = b'''
    x = np.array([0.0] * (n+1))
    linha = n
    if matriz[n][n] == 0: return None
    x[n] = b[linha] / matriz[linha][linha]
    linha -= 1
    while 0 <= linha:
        coef_aux = 0
        linha_aux = linha + 1
        while linha_aux <= n:
            coef_aux += matriz[linha][linha_aux] * x[linha_aux]
            linha_aux += 1
        x[linha] = (b[linha] - coef_aux) / matriz[linha][linha]
        linha -= 1
    return x

def CalculaResultadoMatrizInferior(matriz, b, n):
    '''Recebe uma matriz triangular inferior e um vetor 'b', e retorna o vetor 'x' do sistema Ax = b'''
    x = np.array([0.0] * (n+1))
    linha = 0
    if matriz[0][0] == 0: return None
    x[0] = b[linha] / matriz[linha][linha]
    linha =+ 1
    while linha <= n:
        coef_aux = 0
        linha_aux = linha - 1
        while 0 <= linha_aux:
            coef_aux += matriz[linha][linha_aux] * x[linha_aux]
            linha_aux -= 1
        x[linha] = (b[linha] - coef_aux) / matriz[linha][linha]
        linha += 1
    return x

def CalculaNormaDaDiferenca(resultado, n):
    '''Calcula a norma da diferença do vetor 'x' calculado e o vetor [1] * n'''
    resultado_aux = resultado.copy()
    resultado_correto = np.array([1.0] * n)
    norma_2_vetor = (resultado_aux - resultado_correto)
    norma_2 = np.sqrt(sum(abs(norma_2_vetor)**2))
    return norma_2

if __name__ == "__main__":
    main()