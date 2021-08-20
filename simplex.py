import numpy as np
from numpy.core.numeric import isclose
#import functools
import math

class Simplex:

    def __init__(self, n, m, c, A, b):

        self.n = A.shape[0]
        self.m = A.shape[1]
        self.c = c
        self.A = A
        self.b = b

        self.numVertices = n
        self.numArcos = m

    def print_tableau(self, tipo):

        if(tipo == "simples"):
            print(self.certif_otimalidade_PL," | ",self.c_PL," | ",self.val_obj_PL)
            for i in range (self.n):
                print(self.matriz_transform[i],"|",self.A_tableau[i]," | ",self.b_tableau[i])
            print("\n\n")
        elif(tipo == "extendido"):
            print(self.certif_otimalidade_PL," | ",self.c_PL," | ",self.val_obj_PL)
            print(self.certif_otimalidade_AUX," | ",self.c_AUX," | ",self.val_obj_AUX)
            for i in range (self.n):
                print(self.matriz_transform[i]," | ",self.A_tableau[i]," | ",self.b_tableau[i])
            print("\n\n")
            
    def encontra_sol(self):
        self.bases_PL = []

        solucao = np.zeros(self.numArcos, dtype = "int32")

        canonicos = np.eye(self.n)

        for i in range(self.numArcos): #ISSO PARECE MEIO ERRADO... NÃO SEI SE RESOLVE TODOS OS CASOS
            if (np.isclose(self.c_PL[i], 0)):
                vetor_aux = []

                achou_base = False

                for linha in range(self.n):
                    vetor_aux.append(self.A_tableau[linha][i])
                    
                for j in range(self.n):
                    #if functools.reduce(lambda x, y : x and y, map(lambda p, q: p == q,canonicos[:][j],vetor_aux), True):                     
                    #if ((canonicos[:][j] == vetor_aux).all()):
                    if (canonicos[:][j] == vetor_aux).all():

                        achou_base = True
                        
                        base = i
                        indice_do_1 = j
                        
                        par = (i, indice_do_1)
                        self.bases_PL.append(par)

                        solucao[i] = self.b_tableau[indice_do_1]

                if(not achou_base):
                    solucao[i] = 0 

            else:
                solucao[i] = 0

        solucao_string = ""
        for i in range(self.numArcos):
            if (i != self.numArcos-1):
                solucao_string += str(solucao[i]) + " "
            else:
                solucao_string += str(solucao[i])
        
        print(solucao_string)

    def print_otima(self):

        print(str(int(self.val_obj_PL)))

        self.encontra_sol()
        
        certificado = np.zeros(self.numVertices, dtype = "int32")
        certificado[0] = 1

        #certificado = "1.0000000 " #ACHO Q TUDO VAI SER INTEIRO SEM ISSO DE PRECISAO DE CASAS DECIMAIS
        for i in range(self.numVertices-2):
            certificado[i+1] = self.certif_otimalidade_PL[i]

        certificado_string = ""
        for i in range (self.numVertices):
            if(i != self.numVertices-1):
                certificado_string += str(certificado[i]) + " "
            else:
                certificado_string += str(certificado[i])

        print(certificado_string)

    def resolve_PL(self):
        self.monta_tableau_extendido()
        self.tableau()

    def monta_tableau_extendido(self):        
        self.A_tableau = self.A
        self.b_tableau = self.b
        self.c_PL = self.c

        for i in range(self.m):
            if(not np.isclose(self.c_PL[i],0)):
                self.c_PL[i] *= -1

        self.val_obj_PL = 0

        self.certif_otimalidade_PL = np.zeros(self.n)

        self.matriz_transform = np.eye(self.n)

        for i in range(self.n):
            if(self.b_tableau[i] < 0):
                self.b_tableau[i] *= -1
                self.A_tableau[i] *= -1
                self.matriz_transform[i] *= -1              
    
    def busca_ci_negativo(self, c):
        for j in range (c.shape[0]): #j é coluna
            if(c[j] < 0 and not np.isclose(c[j], 0)):
                return True, j   
        return False, -1

    def tableau(self):
        tem_ci_negativo, j = self.busca_ci_negativo(self.c_PL)

        while(tem_ci_negativo):
            #Achar a linha na coluna j com a menor razao b/A
            razao = math.inf #2**10
            linha_pivo = 0
            coluna_pivo = 0

            for i in range(self.n): #i é linha
                if(self.A_tableau[i][j] > 0 and self.b_tableau[i]/self.A_tableau[i][j] < razao):
                    linha_pivo = i
                    coluna_pivo = j
                    razao = self.b_tableau[i]/self.A_tableau[i][j] 
   
            #Pivotear: Eliminacao Gaussiana de forma que apenas A[linha_pivo][coluna_pivo] seja 1 e o restante da coluna seja 0
            pivo = self.A_tableau[linha_pivo][coluna_pivo]

            self.matriz_transform[linha_pivo] /= pivo
            self.A_tableau[linha_pivo] /= pivo
            self.b_tableau[linha_pivo] /= pivo

            valor_op_PL = -1*self.c_PL[coluna_pivo]

            self.certif_otimalidade_PL += valor_op_PL*self.matriz_transform[linha_pivo]                      
 
            self.c_PL += valor_op_PL*self.A_tableau[linha_pivo] 

            self.val_obj_PL += valor_op_PL*self.b_tableau[linha_pivo]

            for i in range(self.n):
                if(i != linha_pivo):
                    elemento = self.A_tableau[i][coluna_pivo]

                    if(np.sign(elemento) == 1 or np.sign(elemento) == -1):
                        valor_op = -1*elemento
                                
                        self.matriz_transform[i] += valor_op*self.matriz_transform[linha_pivo]
                        self.A_tableau[i] += valor_op*self.A_tableau[linha_pivo]
                        self.b_tableau[i] += valor_op*self.b_tableau[linha_pivo]
                    else:
                        continue

            tem_ci_negativo, j = self.busca_ci_negativo(self.c_PL)

        self.print_otima()