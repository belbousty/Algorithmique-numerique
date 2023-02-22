import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import block_diag
import math

def M_generator(N:int):
    M = np.zeros((N,N))
    for i in range (0,N-1):
        M[i,i] = -4
        M[i,i+1] = 1
        M[i+1,i] = 1
    M[N-1,N-1]= -4
    return M

def A_generator(N): # generateur de la matrice A de l'équation AX=B
    M =M_generator(N)
    A_prime = block_diag(*([M]*N))   
    for i in range(0,N**2-N):
        A_prime[i+N,i] = 1
        A_prime[i,i+N] = 1
        
    return A_prime
    
def transform_Vec_Mat(X, N): #transformation d'un vecteur de taille N**2 en une matrice N*N
    M = np.zeros((N,N))
    for i in range(0,N):
        for j in range(0,N):
            M[i,j] = X[i*N+j,0]
    return M

def transform_Mat_Vec(M, N): #transformation d'une matrice N*N en un vecteur de taille N**2
    X = np.zeros((N**2,1))
    for j in range(0,N):
        for i in range(0,N):
            X[i*N+j,0] = M[i,j]
    return X

def B_generator_center(N): #generator of f values //center

    Z = np.zeros((N,N))
    for i in range(2,N-1):
        for j in range(i,N-i):
            Z[i,j] = i/N**5
            Z[j,i] = i/N**5
            Z[N-1-i, j] = i/N**5
            Z[j, N-1-i] = i/N**5
    return transform_Mat_Vec(Z,N)


def factorisation_incomplete_cholesky(A): 
  (ligne, colonne) = np.shape(A)
  T = np.zeros((ligne, colonne))

  for j in range(colonne):
    for i in range(j+1):
      sigma = sum(T[i][k]*T[j][k] for k in range(i))
      if(A[j][i] == 0 ):
        T[j][i] = 0
      elif(A[i][i] == 0):
        T[i][i] = 0
      else :
        if (i == j):
          T[i][i] = math.sqrt(abs(A[i][i] - sigma))
        else:
          T[j][i] = (A[i][j] - sigma)/T[i][i]
      
  return T    

 ################################## solution en utilisant la biblio numpy
def solution_center(N): # solution pour un radiateur au centre
    A = A_generator(N)
    B = np.zeros((N, N))
    B[N//2, N//2] = -1
    B = transform_Mat_Vec(B,N)
    Z = np.linalg.solve(A,B)

    return transform_Vec_Mat(Z,N)

plt.imshow(solution_center(13), cmap='hot', interpolation='bicubic')
plt.title("cas d'un radiateur placé au centre du carré [numpy]")
plt.show()


def solution_north(N): # solution pour un mur chaud au nord
    A = A_generator(N)
    B = np.zeros((N, N))   
    B[N//4, N//2] = -1
    B = transform_Mat_Vec(B,N)
    Z = np.linalg.solve(A,B)
    return transform_Vec_Mat(Z,N)

plt.imshow(solution_north(13), cmap='hot', interpolation='quadric')
plt.title("cas d'un mur chaud placé au nord du carré [numpy]")
plt.show()

#################################

def factorisation_cholesky(A):
  (ligne, colonne) = np.shape(A)
  T = np.zeros((ligne, colonne))

  for j in range(colonne):
    for i in range(j+1):
      sigma = sum(T[i][k]*T[j][k] for k in range(i))
      if (i == j):
        T[i][i] = math.sqrt(abs(A[i][i] - sigma))
      else:
        T[j][i] = (A[i][j] - sigma)/T[i][i]
  return T 

################################# utilisation de la partie 1 

def cholesky_method_center(N):
    A = A_generator(N)
    T = factorisation_cholesky(-1*A)
    B = np.zeros((N, N))   
    B[N//2, N//2] = 1
    B = transform_Mat_Vec(B,N)
    Y = np.linalg.solve(T, B)
    X = np.linalg.solve(T.T,Y)
    return transform_Vec_Mat(X,N)

plt.imshow(cholesky_method_center(13), cmap='hot', interpolation='quadric')
plt.title("cas d'un radiateur placé au centre du carré")
plt.show()

def cholesky_method_north(N):
    A = A_generator(N)
    T = factorisation_cholesky(-1*A)
    B = np.zeros((N, N))   
    B[N//4, N//2] = 1
    B = transform_Mat_Vec(B,N)
    Y = np.linalg.solve(T, B)
    X = np.linalg.solve(T.T,Y)
    return transform_Vec_Mat(X,N)

plt.imshow(cholesky_method_north(13), cmap='hot', interpolation='quadric')
plt.title("cas d'un mur chaud placé au nord du carré")
plt.show()

################################################ Erreur 
N=2
A = A_generator(N)
T = factorisation_cholesky(-1*A)
A_facto = -1*np.dot(T,T.T)
z=0
x= np.linspace(0,200,N**4)
erreur_A = np.zeros((1,N**4))
for i in range (0,N**2):
  for j in range (0,N**2):
    if (A[i,j] != 0):
      erreur_A[0,i*N+j] = (abs((A_facto[i,j])-(A[i,j])/abs(A[i,j])))
    if erreur_A[0,i*N+j] != 0 :
      z+=1
print(z)

plt.plot(x, erreur_A[0],label="erreur de la décomposition")
plt.xscale("log")
plt.show()