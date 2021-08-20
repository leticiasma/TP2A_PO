import numpy as np
from simplex import Simplex
from numpy.core.numeric import isclose

n,m = input().split() #n vertices, m arcos
n = int(n)
m = int(m)
capacidades = np.array(input().split(), dtype='float32') #vetor de capacidades c, com m inteiros

N =  np.zeros((n, m))

b = np.zeros(n-2)
b = np.concatenate((b, capacidades), axis=0)

for i in range(n):
    linha = np.array(input().split())
    N[i] = linha

c_PL = N[0]
zeros = np.zeros(m)
c_PL = np.concatenate((c_PL, zeros), axis=0)

for i in range(2*m):
    if(not np.isclose(c_PL[i],0)):
        c_PL[i] *= -1

identidade = np.eye(m)
identidades = np.concatenate((identidade, identidade), axis=1)

#Deletando os vértices "s" e "t" da matriz
M = np.delete(N, 0, 0)
M = np.delete(M, -1, 0)
matriz_zeros = np.zeros((n-2,m))
M_zeros = np.concatenate((M, matriz_zeros), axis=1)

matriz_restricoes = np.concatenate((M_zeros, identidades), axis=0)

simplex = Simplex(n,m,c_PL,matriz_restricoes,b)

simplex.resolve_PL()