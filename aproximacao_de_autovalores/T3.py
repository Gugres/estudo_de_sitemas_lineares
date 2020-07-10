import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as plt
from scipy.io import mmread
from scipy import sparse as scs
import time
import sys

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

def Busca(x):
    p = 0
    max = abs(x[p])
    for i in range(len(x)):
        if abs(x[i]) > max:
            p = i
    return p

def AproximacaoAutoValorDominante(matriz, rows, cols):
    '''Aproxima o auto-valor dominante usando o circulo de Gersgorin, usando o formato COO da matriz'''
    q = 0
    limite_inferior = 0
    limite_superior = 0
    somatoria = 0
    a,i,j = (0,0,0)
    while i < len(rows):
        if rows[i] == j and rows[i] != cols[i]:
            somatoria += abs(matriz[rows[i]][cols[i]])
            i += 1
            continue
        if rows[i] == j and rows[i] == cols[i]:
            a = matriz[rows[i]][cols[i]]
            i += 1
            continue
        limite_inferior = a - somatoria
        limite_superior = somatoria + a
        if limite_superior > q: q = limite_superior
        j += 1
        somatoria = 0
    return q + 1

def main():
    matrizes = ["bcsstm20/bcsstm20.mtx", "bcsstk10/bcsstk10.mtx", "bcsstk05/bcsstk05.mtx", "mesh1em6/mesh1em6.mtx",
                "mesh1e1/mesh1e1.mtx", "mesh2e1/mesh2e1.mtx", "Trefethen_20/Trefethen_20.mtx", "Trefethen_700/Trefethen_700.mtx",
                "nos4/nos4.mtx", "nasa1824/nasa1824.mtx"]
    nao_nulos = []
    esparsidade = []
    tamanho = []
    tempos_mps = []
    tempos_mpi = []
    tempos_qr = []
    inter_mps = []
    inter_mpi = []
    raio_esp_mps = []
    raio_esp_mpi = []
    raio_esp_qr = []
    menor_autoValor_qr = []
    for nome_matriz in matrizes:
        # Le o arquivo e ordena a representação COO da matriz esparsa por linha
        sparse_matrix_coo = mmread("./" + nome_matriz)
        # Gera a forma completa da matriz esparsa
        matriz = sparse_matrix_coo.toarray()
        # Gera um array de todas as posições "j" com elementos não nulos
        rows, cols = OrdenaCOO(matriz)
        # Vetor inicial
        x = np.ones(len(matriz))
        print("\nmatriz: ", nome_matriz)
        print("Não nulos:", len(rows))
        print("Esparsidade: ", 1 - (len(rows) / (matriz.shape[0]*matriz.shape[1])))
        print("Tamanho da matriz: ", len(matriz))
        nao_nulos.append(len(rows))
        esparsidade.append(1 - (len(rows) / (matriz.shape[0]*matriz.shape[1])))
        tamanho.append(len(matriz))
        # Configurações
        tol = 10**(-3)
        max_it = 300
        # Aplicando o metodo da potência simétrico
        tempo_1_mps = time.time()
        mu_mps, x_mps, k_mps = MetodoDaPotenciaSimetrico(matriz, rows, cols, x, tol, max_it)
        tempo_2_mps = time.time()
        delta_mps = tempo_2_mps - tempo_1_mps
        tempos_mps.append(delta_mps)
        inter_mps.append(k_mps)
        raio_esp_mps.append(mu_mps)
        # Metodo da potencia inversa
        tempo_1_mpi = time.time()
        mu_mpi, x_mpi, k_mpi = MetodoDaPotenciaInversaNorma2(matriz, rows, cols, x, tol, max_it)
        tempo_2_mpi = time.time()
        delta_mpi = tempo_2_mpi - tempo_1_mpi
        tempos_mpi.append(delta_mpi)
        inter_mpi.append(k_mpi)
        raio_esp_mpi.append(mu_mpi)
        # Decomposicao QR para encontrar o espectro completo
        tempo_1_qr = time.time()
        T, U = DecomposiçãoQR(matriz)
        tempo_2_qr = time.time()
        delta_qr = tempo_2_qr - tempo_1_qr
        # Encontrando o maior e o menor auto-valor
        raio_espectral = 0
        menor_auto_valor = T[0][0]
        for i in range(len(T)):
            if abs(T[i][i]) > raio_espectral:
                raio_espectral = T[i][i]
            if abs(T[i][i]) < menor_auto_valor:
                menor_auto_valor = T[i][i]
        tempos_qr.append(delta_qr)
        raio_esp_qr.append(raio_espectral)
        menor_autoValor_qr.append(menor_auto_valor)
        print("\nRaio espectral: ", raio_espectral)
        print("Menor auto-valor: ", raio_espectral)
    df = pd.DataFrame({"matrizes": matrizes, "dimensao": tamanho, "nao_nulos": nao_nulos, "esparsidade": esparsidade,
                        "tempos_metPotSim": tempos_mps, "interacoes_metPotSim": inter_mps, "raio_esp_metPotSim": raio_esp_mps, 
                        "tempos_metPotInv": tempos_mpi, "interacoes_metPotInv": inter_mpi, "raio_esp_metPotInv": raio_esp_mpi,
                        "tempos_QR": tempos_qr, "raio_esp_QR": raio_esp_qr, "menor_autoValor_QR": menor_autoValor_qr})
    df.to_csv("analise_resultados.csv")

def MetodoDaPotencia(matriz, rows, cols, x, tol, max_it):
    k = 1
    x_aux = x.copy()
    p = Busca(x_aux)
    x_aux = x_aux / x_aux[p]
    while k <= max_it:
        y = ProdutoVetorMatrizCOO(matriz, rows, cols, x_aux)
        mu = y[p]
        p = Busca(y)
        if y[p] == 0:
            return (0, x_aux, k)
        ERR = max(x_aux - (y / y[p]))
        x_aux = y / y[p]
        if ERR < tol:
            return (mu, x_aux, k)
        k = k + 1
    return (mu, x_aux, k)

def MetodoDaPotenciaModificado(matriz, rows, cols, x, tol, max_it):
    k = 1
    mu_0 = 0
    mu_1 = 0
    x_aux = x.copy()
    p = Busca(x_aux)
    x_aux = x_aux / x_aux[p]
    while k <= max_it:
        y = ProdutoVetorMatrizCOO(matriz, rows, cols, x_aux)
        mu_2 = y[p]
        mu = mu_0 - ((mu_1 - mu_0) ** 2) / (mu_2 - 2*mu_1 + mu_0)
        p = Busca(y)
        if y[p] == 0:
            return (0, x_aux, k)
        ERR = max(x_aux - (y / y[p]))
        x_aux = y / y[p]
        if ERR < tol:
            return (mu, x_aux, k)
        k = k + 1
        mu_0 = mu_1
        mu_1 = mu_2
    return (mu, x_aux, k)

def MetodoDaPotenciaSimetrico(matriz, rows, cols, x, tol, max_it):
    k = 1
    mu_0 = 0
    mu_1 = 0
    x_aux = x.copy()
    x_aux = x_aux / ((sum(x_aux**2))**(1/2))
    while k <= max_it:
        y = ProdutoVetorMatrizCOO(matriz, rows, cols, x_aux)
        mu_2 = x_aux.dot(y)
        mu = mu_0 - ((mu_1 - mu_0) ** 2) / (mu_2 - 2*mu_1 + mu_0)
        if (sum(y**2))**(1/2) == 0:
            return (0, x_aux, k)
        ERR_aux = (x_aux - (y / ((sum(y**2))**(1/2))))
        ERR = (sum(ERR_aux**2))**(1/2)
        x_aux = y / ((sum(y**2))**(1/2))
        if ERR < tol:
            return (mu, x_aux, k)
        k = k + 1
        mu_0 = mu_1
        mu_1 = mu_2
    return (mu, x_aux, k)

def MetodoDaPotenciaInversaNormaInfinito(matriz, rows, cols, x, tol, max_it):
    mu_0 = 0
    mu_1 = 0
    n = len(matriz)
    x_aux = x.copy()
    matriz_aux = matriz.copy()
    # Usando o circulo de Gersgorin
    q = AproximacaoAutoValorDominante(matriz, rows, cols)
    k = 1
    p = Busca(x_aux)
    x_aux = x_aux / x_aux[p]
    B = (matriz_aux - q*np.identity(n))
    while k <= max_it:
        y = np.linalg.solve(B, x_aux)
        if np.array(y).any == None:
            return (mu, x_aux, k)
        mu_2 = y[p]
        mu = mu_0 - ((mu_1 - mu_0) ** 2) / (mu_2 - 2*mu_1 + mu_0)
        p = Busca(y)
        ERR = max(x_aux - (y/y[p]))
        x_aux = y / y[p]
        if ERR < tol:
            mu = (1/mu) + q
            return (mu, x_aux, k)
        k += 1
        mu_0 = mu_1
        mu_1 = mu_2
    return (mu, x_aux, k)

def MetodoDaPotenciaInversaNorma2(matriz, rows, cols, x, tol, max_it):
    mu_0 = 0
    mu_1 = 0
    n = len(matriz)
    x_aux = x.copy()
    matriz_aux = matriz.copy()
    k = 1
    # Usando o circulo de Gersgorin
    q = AproximacaoAutoValorDominante(matriz, rows, cols)
    x_aux = x_aux / ((sum(x_aux**2))**(1/2))
#     B = (matriz_aux - q*np.identity(n))
    B = scs.csr_matrix(matriz_aux - q*np.identity(n))
    while k <= max_it:
#         y = np.linalg.solve(B, x_aux)
        y = scs.linalg.spsolve(B,x_aux)
        if np.array(y).any == None:
            return (mu, x_aux, k)
        mu_2 = y.dot(x_aux)
        mu = mu_0 - ((mu_1 - mu_0) ** 2) / (mu_2 - 2*mu_1 + mu_0)
        ERR_aux = x_aux + (y / ((sum(y**2))**(1/2)))
        ERR = (sum(ERR_aux**2))**(1/2)
        x_aux = y / ((sum(y**2))**(1/2))
#         print("\ny: ", y)
#         print("x: ", x_aux)
        if ERR < tol:
            mu = (1/mu) + q
            return (mu, x_aux, k)
        k += 1
        mu_0 = mu_1
        mu_1 = mu_2
    return (mu, x_aux, k)

def DecomposiçãoQR(matriz):
    T = matriz.copy()
    U = np.identity(len(matriz))
#     U = scs.identity(len(matriz))
    i = 0
    while i < 60:
        Q,R = np.linalg.qr(T)
#         Q = scs.csr_matrix(Q)
#         R = scs.csr_matrix(R)
#         T = R.dot(Q).toarray()
        T = R.dot(Q)
        U = U.dot(Q)
        i += 1
    return (T,U)

if __name__ == "__main__":
    main()