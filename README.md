# estudo_de_sitemas_lineares

  Objetiva realizar um estudo de soluções de sistemas linares, através de alguns métodos para diminuição do erro numérico e aproveitamento da forma da matriz.
  
  * Eliminição de Gauss sem pivotamento
  * Eliminição de Gauss com pivotamento parcial
  * Eliminição de Gauss com pivotamento parcial escalado
  * Fatoração de Cholesky

# Requerimentos

  São utilizadas quatro bibliotecas para análise dos resultados e geração das matrizes:
  
  * numpy
  * pandas
  * time
  * sys

  Para instalação das bibliotecas que não são internas no python, utilize os seguintes comandos:
  
```
pip install numpy
pip install pandas
```

# Testando a aplicação

  Para testar a aplicação, basta utilizar um dos quatro comandos a seguir dentro do diretorio:
 
```
python T1.py 1 0 # Gera uma base de dados com 100 resoluções dos metodos de eliminação de Gauss com matrizes diferentes de Hilbert
python T1.py 1 1 # Testa os metodos de eliminação de Gauss com apenas uma matriz mocada
python T1.py 2 0 # Gera uma base de dados com 71 resolucoes de matrizes aleatorias pela fatoração de Cholesky e metodo de eliminacao de Gauss
python T1.py 2 1 # Testa a fatoração de Cholesky com apenas uma matriz mocada
```

# Bibliografia

  Todos os códigos apresentados nos algoritmos, bem como suas devidas explicações, encontram-se no livro ["_Numerical Analyses_"](https://www.google.com/search?sxsrf=ALeKk02kn-IlANTHneIM5UBjtV09U8a8og%3A1588687297722&ei=wXGxXt2zK7zO5OUPzdGWsAo&q=numerical+analysis+burden&oq=numerical+analyses+bur&gs_lcp=CgZwc3ktYWIQAxgAMgQIABANMgQIABANMgQIABANMggIABAWEAoQHjoECAAQR1ClG1iyIGC_K2gAcAJ4AIAB_AKIAdILkgEDMy00mAEAoAEBqgEHZ3dzLXdpeg&sclient=psy-ab), Burden e Faires, nas páginas 358-379 e 412-419.
