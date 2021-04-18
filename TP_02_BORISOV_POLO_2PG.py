"""
BORISOV Alexander 

POLO Mateo 

2PG

"""

import numpy as np
from math import *
import matplotlib.pyplot as plt
import time
import random

print()
print()
print("############### Décomposition de Cholesky ###############")
print()
print()
def Cholesky(A):
    n,m=A.shape
    L=np.zeros((n,m))
    for i in range(0,n):
        for j in range(i,n):
            if(i==j):
                S=0
                for k in range(0,i):
                    S=S+L[i][k]**2
                L[i][i]=sqrt(A[i][i]-S)
            else:
                S=0
                for k in range(0,i):                
                    S=S+(L[i][k])*(L[j][k])
                L[j][i]=(A[i][j]-S)/(L[i][i])
    return(L)
A=np.array([[4,-2,-4],[-2,10,5],[-4,5,6]])


print()
L=Cholesky(A)
print("L=",L)
print()


print()
print()
print("############### Résolution de système ###############")
print()
print()

def ResolutionSystTriInf(Taug):
    n,m=Taug.shape
    if m != n+1:
        print("Ce n'est pas une matrice colonne")
        return
    x=np.zeros(n)
    for i in range(n):
        somme=0
        for k in range(i):
            somme+=x[k]*Taug[i,k]
        x[i]=(Taug[i,-1]-somme)/Taug[i,i]
    
    return(x)



def ResolutionSystTriSup(Taug):
    n,m=Taug.shape
    if m != n+1:
        print("Ce n'est pas une matrice colonne")
        return
    x=np.zeros(n)
    for i in range(n-1,-1,-1):
        somme=0
        for k in range(i,n):
            somme+=x[k]*Taug[i,k]
        x[i]=(Taug[i,-1]-somme)/Taug[i,i]
    
    return(x)


def ResolCholesky(A,B):
    L = Cholesky(A)
    L_T = np.transpose(L)
    
    Aaug = np.hstack([L,B])
    Y = ResolutionSystTriInf(Aaug)
    Y = Y.reshape(-1,1)
    
    Baug = np.hstack([L_T,Y])
    X = ResolutionSystTriSup(Baug)
    return X


B = np.array([[6],[-9],[-7]])
ResolCholesky(A,B)

# Cholesky avec np.linalg.cholesky       
def ResolCholeskyMACHINE(A,B):
    L = np.linalg.cholesky(A)
    L_T = np.transpose(L)
    
    Aaug = np.hstack([L,B])
    Y = ResolutionSystTriInf(Aaug)
    Y = Y.reshape(-1,1)
    
    Baug = np.hstack([L_T,Y])
    X = ResolutionSystTriSup(Baug)
    return X




print ("La solution est : X = ", ResolCholesky(A,B))

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


#Courbes du temps 

TempsCHO = []
IndicesCHO = []
ErreurCHO = []
TempsLU = []
IndicesLU = []
ErreurLU = []
TempsCHO_2 = []
IndicesCHO_2 = []
ErreurCHO_2 = []
TempsLIN = []
IndicesLIN = []
ErreurLIN = []

for n in range(50,550,50):
    try:
        M= np.random.randint(low = 1,high = 10,size = (n,n))
        while np.linalg.det(M)==0 :
            M= np.random.randint(low=1,high=10,size=(n,n))
        A = M.dot(M.T)
        B= np.random.randint(low = 1,high = 10,size = (n,1))
        A = np.array(A, dtype = float)
        B = np.array(B, dtype = float)
        t1 = time.time()
        X = ResolCholesky(A,B)
        t2 = time.time()
        T1 = t2 - t1
        erreur = np.linalg.norm(A.dot(X)-B)
        TempsCHO.append(T1)
        IndicesCHO.append(n)
        ErreurCHO.append(erreur)
        

    except:
        print('')
        
for n in range(50,550,50):
    try:
        A= np.random.randint(low = 1,high = n,size = (n,n))
        B= np.random.randint(low = 1,high = n,size = (n,1))
        A = np.array(A, dtype = float)
        B = np.array(B, dtype = float)
        t3 = time.time()
        X = ResolutionLU(A,B)
        t4 = time.time()
        T2 = t4 - t3
        erreur2 = np.linalg.norm(A.dot(X)-B)
        TempsLU.append(T2)
        IndicesLU.append(n)
        ErreurLU.append(erreur2)


    except:
        print('')
        

for n in range(50,550,50):
    try:
        M = np.random.randint(low = 1,high = 10,size = (n,n))
        while np.linalg.det(M)==0 :
            M= np.random.randint(low=1,high=10,size=(n,n))
        A = M.dot(M.T)
        B= np.random.randint(low = 1,high = 10,size = (n,1))
        A = np.array(A, dtype = float)
        B = np.array(B, dtype = float)
        t5 = time.time()
        X = ResolCholeskyMACHINE(A,B)
        t6 = time.time()
        T3 = t6 - t5
        erreur3 = np.linalg.norm(A.dot(X)-B)
        TempsCHO_2.append(T3)
        IndicesCHO_2.append(n)
        ErreurCHO_2.append(erreur3)

        
    except:
        print('')
        
for n in range(50,550,50):
    try:
        A= np.random.randint(low = 1,high = n,size = (n,n))
        B= np.random.randint(low = 1,high = n,size = (n,1))
        A = np.array(A, dtype = float)
        B = np.array(B, dtype = float)
        t7 = time.time()
        X = np.linalg.solve(A,B)
        t8 = time.time()
        T4 = t8 - t7
        erreur4 = np.linalg.norm(A.dot(X)-B)
        TempsLIN.append(T4)
        IndicesLIN.append(n)
        ErreurLIN.append(erreur4)


    except:
        print('')
        


x7 = IndicesCHO
y7 = TempsCHO
x2 = IndicesLU
y2 = TempsLU
x3 = IndicesCHO_2
y3 = TempsCHO_2
x4 = IndicesLIN
y4 = TempsLIN



plt.plot(x7,y7,color='green', label="Cholesky")
plt.plot(x2,y2,color='red', label="DécompostionLU")
plt.plot(x3,y3,color='blue', label="Cholesky_Machine")
plt.plot(x4,y4,color='orange', label="Linalg.solve")
plt.legend()
plt.xlabel('Dimension Matrice')
plt.ylabel('Temps en seconde')
plt.title("Graphique présentant les courbes logarithmiques du temps en fonction de n ")
plt.show()


# Courbes d'erreurs 
'''x7 = IndicesCHO
z7 = ErreurCHO
x2 = IndicesLU
z2 = ErreurLU
x3 = IndicesCHO_2
z3 = ErreurCHO_2
x4 = IndicesLIN
z4 = ErreurLIN

#plt.plot(x7,z7,color='green', label="Cholesky")
plt.plot(x2,z2,color='red', label="DécompostionLU")
#plt.plot(x3,z3,color='blue', label="Cholesky_Machine")
#plt.plot(x4,z4,color='orange', label="Linalg.solve")
plt.legend()
plt.xlabel('Dimension Matrice')
plt.ylabel('Erreur || AX - B ||')
plt.title("Graphique présentant les courbes d'erreur en fonction de n ")
plt.show()'''