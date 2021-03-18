# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

import numpy as np
from math import *
import matplotlib.pyplot as plt
import time
import random


print("ALGORITHME DE GAUSS")
print()
print()



def ResolutionSystTriInf(Taug):
   n,m=  Taug.shape
   if m != n+1:
       print("ce n'est pas une matrice augmentée")
       return
   x = np.zeros(n)
   for i in range(n):
       somme = 0
       for k in range(i):
           somme = somme + x[k]*Taug[i,k]
       x[i] = (Taug[i,-1]-somme)/Taug[i,i]
   return x



def ResolutionSystTriSup(Taug):
   n,m = Taug.shape
   if m != n+1:
       print("ce n'est pas une matrice augmentée")
       return
   x = np.zeros(n)
   for i in range(n-1,-1,-1):
       somme = 0
       for k in range(i,n):
           somme = somme + x[k]*Taug[i,k]
       x[i] = (Taug[i,-1]-somme)/Taug[i,i]
   return x




def ReductionGauss(Aaug):
    n,m = Aaug.shape
    for i in range(0,n-1):
        if Aaug[i,i] != 0:
           for k in range (i+1,n):
               g = Aaug[k,i]/Aaug[i,i]
               for j in range (0,n+1):
                   Aaug[k,j] = Aaug[k,j]-g*Aaug[i,j]
        else:
            print("Pivot nul, veuillez choisir une autre matrice")
            return
    return(Aaug)
    
Aaug=np.array([[2,5,6,7],[4,11,9,12],[-2,-8,7,3]])
Taug = ReductionGauss(Aaug)



print("\n")
print("Résolution avec la fonction Gauss(A,B) \n")


def Gauss(A,B):

    Aaug=np.concatenate((A,B.T),axis=1)
    #print('Matrice augmenté :')
    #print(Aaug)
    Taug=ReductionGauss(Aaug)
    x=ResolutionSystTriSup(Taug)
    return x
    

A=np.array([[2,5,6],[4,11,9],[-2,-8,7]])
B=np.array([[7,12,3]])
#print ('solution : \n',Gauss(A,B))







print()
print()
print("############### Décomposition LU ###############")
print()
print()



def DecompostionLU(A):
    n,m = A.shape
    U = np.zeros((n,m))      
    L = np.eye(n)
    for i in range (0,n):
        for k in range (i+1,n):
            pivot = A[i][i]
            pivot = A[k][i]/pivot
            L[k][i] = pivot
            for j in range (i,n):
                A[k][j]= A[k][j]-(pivot*A[i][j])
        U = A
    return (L,U)
    

A =np.array ([[2,5,6],[4,11,9],[-2,-8,7]])
B =np.array([[7],[12],[3]])


def ResolutionLU(A,B):
    n,m = A.shape
    L,U = DecompostionLU(A)
    
    #print("L= \n",L)
    #print("U=\n",U)
    

    Y = ResolutionSystTriInf(np.concatenate((L,B),axis=1))
    Y = np.asarray(Y).reshape(n,1)
    X = ResolutionSystTriSup(np.concatenate((U,Y),axis=1))
    
    #print("Y=\n",Y)
    #print("X=\n",X)
    return(X)
    
    
X = ResolutionLU(A,B)





   

 

#print("X=\n",ResolutionLU(A,b))
    

print()
print()
print("############### Gauss Choix Pivot Partiel ###############")
print()
print()


A=np.array([[1,2,2],[2,1,2],[2,2,1]])
B=np.array([[1],[2],[1]])

def GaussChoixPivotPartiel(A,B):
    X=np.hstack([A,B])  
    n,m = X.shape
    if m!= n+1:
        print("erreur de format")
    else:
        for k in range(n-1):
            for i in range(k+1,n):
                for var in range(i,n):
                    if abs(X[var,k]) > abs(X[k,k]):
                        L0 = np.copy(X[k,:])
                        X[k,:] = X[var,:]
                        X[var,:] = L0
                if X[k,k] == 0:
                    print('erreur pivot nul')
                g = X[i,k]/X[k,k]
                X[i,:] = X[i,:] - g*X[k,:]
    
    X=ResolutionSystTriSup(X)

    return(X)

#print(GaussChoixPivotPartiel(A,B))

print()
print()
print("############### Gauss Choix Pivot Total###############")
print()
print()


def GaussChoixPivotTotal(A,B):
    X=np.hstack([A,B])  
    n,m = X.shape
    if m!= n+1:
        print("erreur de format")
    else:
        for k in range(m-1):
            for i in range(k+1,m-1):
                for var in range(i,m-1):
                    if abs(X[var,k]) > abs(X[k,k]):
                        L0 = np.copy(X[k,:])
                        X[k,:] = X[var,:]
                        X[var,:] = L0
                if X[k,k] == 0:
                    print('erreur pivot nul')
                g = X[i,k]/X[k,k]
                X[i,:] = X[i,:] - g*X[k,:]
    
    X=ResolutionSystTriSup(X)

    return(X)

#print(GaussChoixPivotTotal(A,B))









# Courbes du temps en fonction de la dimension de la matrice 

Temps = []
Indices = []
TempsLU = []
IndicesLU = []
TempsPGP = []
IndicesPGP = []
TempsPT=[]
IndicesPT=[]
TempsLIN = []
IndicesLIN = []


'''for n in range(50,550,50):
    try:
        A = np.random.randint(low = 1,high = n,size = (n,n))
        B = np.random.randint(low = 1,high = n,size = (1,n))
        A = np.array(A, dtype = float)
        B = np.array(B, dtype = float)
        t1 = time.time()
        Gauss(A,B)
        t2 = time.time()
        T1 = t2 - t1
        Temps.append(T1)
        Indices.append(n)

    except:
        print('')
            
for n in range(50,550,50):
    try:
        A= np.random.randint(low = 1,high = n,size = (n,n))
        B= np.random.randint(low = 1,high = n,size = (n,1))
        A = np.array(A, dtype = float)
        B = np.array(B, dtype = float)
        t3 = time.time()
        ResolutionLU(A,B)
        t4 = time.time()
        T2 = t4 - t3
        TempsLU.append(T2)
        IndicesLU.append(n)

    except:
        print('')

for n in range(50,550,50):
    try:
        A= np.random.randint(low = 1,high = n,size = (n,n))
        B= np.random.randint(low = 1,high = n,size = (n,1))
        A = np.array(A, dtype = float)
        B = np.array(B, dtype = float)
        t5 = time.time()
        GaussChoixPivotPartiel(A,B)
        t6 = time.time()
        T3 = t6 - t5
        TempsPGP.append(T3)
        IndicesPGP.append(n)
        

    except:
        print('')
        
for n in range(50,550,50):
    try:
        A= np.random.randint(low = 1,high = n,size = (n,n))
        B= np.random.randint(low = 1,high = n,size = (n,1))
        A = np.array(A, dtype = float)
        B = np.array(B, dtype = float)
        t7 = time.time()
        GaussChoixPivotTotal(A,B)
        t8 = time.time()
        T4 = t8 - t7
        TempsPT.append(T4)
        IndicesPT.append(n)
        

    except:
        print('')'''
        
for n in range(50,550,50):
    try:
        A= np.random.randint(low = 1,high = n,size = (n,n))
        B= np.random.randint(low = 1,high = n,size = (n,1))
        A = np.array(A, dtype = float)
        B = np.array(B, dtype = float)
        t9 = time.time()
        np.linalg.solve(A,B)
        t10 = time.time()
        T5 = t10 - t9
        TempsLIN.append(T5)
        IndicesLIN.append(n)

    except:
        print('')




x1 = Indices
y1 = Temps
x2 = IndicesLU
y2 = TempsLU
x3 = IndicesPGP
y3 = TempsPGP
x4 = IndicesPT
y4 = TempsPT
x5 = IndicesLIN
y5 = TempsLIN

#plt.plot(x1,y1,color='green', label="Gauss")
#plt.plot(x2,y2,color='red', label="DécompostionLU")
#plt.plot(x3,y3,color='blue', label="PivotPartiel")
#plt.plot(x4,y4,color='black', label="PivotTotal")
plt.plot(x5,y5,color='orange', label="Linalg.solve")
plt.legend()
plt.xlabel('dimension matrice')
plt.ylabel('temps en secondes')
plt.title("Graphique présentant les courbes du temps en fonction de n ")
plt.show()

# Courbes des erreurs en fonction de la dimension de la matrice 


'''Erreurs = []

for n in range(0,6,2):
    A= np.random.randint(low = 1,high = n,size = (n,n))
    B= np.random.randint(low = 1,high = n,size = (1,n))
    X = Gauss(A,B)
    L = np.linalg.norm(A.dot(X)-B)
    Erreurs.append(L)
    Indices.append(n)

x6 = Erreurs
y6 = Indices
plt.plot(x6,y6,color='green')
plt.title("Représentation erreurs de calculs")
plt.xlabel('Dimension Matrice')
plt.ylabel('Erreurs')
plt.show()'''





 



