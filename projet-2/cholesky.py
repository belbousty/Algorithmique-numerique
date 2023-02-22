import numpy as np
import math
import random
import matplotlib.pyplot as plt


# Question 1 
def factorisation_cholesky(A):
  """
  param : A matrix
  returns : T Lower triangulaire matrix

  complexity : n**(3) 
  """
  (ligne, colonne) = np.shape(A)
  T = np.zeros((ligne, colonne))

  for j in range(colonne):
    for i in range(j+1):
      sigma = sum(T[i][k]*T[j][k] for k in range(i))
      if (i == j):
        T[i][i] = math.sqrt((A[i][i] - sigma))
      else:
        T[j][i] = (A[i][j] - sigma)/T[i][i]
  return T 

#Question 3
def generate_matrix_sym(n) :
    M = np.random.rand(n,n)
    A = np.dot(M,M.T)
    for i in range(n):
      for j in range(n):
        if A[i][j] == 0:
          generate_matrix_sym(n)
    return A
    

def matrice_copy(A):
  (n,n) = np.shape(A)
  M = np.zeros([n,n])
  for i in range (n):
    for j in range(n):
      M[i][j] = A[i][j]
  return M

def definite_matrix(n,m) :
  """
  param : n, the size of the matrix
          m, number of non zeros in the matrix
  returns : Definite matrix
    """
  B = generate_matrix_sym(n)
  M = matrice_copy(B)
  count = 0
  while count < n**2 - n - m :
      i = np.random.randint(0,n)
      j = np.random.randint(0,n)
      if i != j and M[i][j] != 0 :
          M[i][j] = 0
          M[j][i] = 0
          count += 2
  for i in np.linalg.eigvals(M) :
      if i <= 0 :
          matrice_creuse(B,n,m)
          break
  return M

def matrice_creuse(B,n,m) :
  """
  param : n, the size of the matrix
          m, number of non zeros in the matrix
  returns : Definite matrix
    """
  M = matrice_copy(B)
  count = 0
  while count < n**2 - n - m :
      i = np.random.randint(0,n)
      j = np.random.randint(0,n)
      if i != j and M[i][j] != 0 :
          M[i][j] = 0
          M[j][i] = 0
          count += 2
  for i in np.linalg.eigvals(M) :
      if i <= 0 :
          matrice_creuse(B,n,m)
          break
  return M

#Question 4
def factorisation_incomplete_cholesky(A):
  """
  param : A matrix
  returns : T Lower triangulaire matrix

  complexity : n**(3) 
  """
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
          T[i][i] = math.sqrt((A[i][i] - sigma))
        else:
          T[j][i] = (A[i][j] - sigma)/T[i][i]
      
  return T 

#Question 4
def factorisation_incomplete_cholesky_abs(A):
  """
  param : A matrix
  returns : T Lower triangulaire matrix

  complexity : n**(3) 
  This function includes an abs to prevent math error while calculating the square root; nevertheless it doesn't give
  the good results, so we just use it in the 3rd part of the project (heat equation) where it solves our problem
  """
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

#Question 5

def plot_cond(A):
  A1 = factorisation_cholesky(A)

  condA = np.linalg.cond(A)
  cond1 = np.linalg.cond(np.dot(A1,A1.T))
  pre_cond1 = np.linalg.cond(np.dot(np.linalg.inv(np.dot(A1,A1.T)),A))

  cond_name = ['A','tr(T).T','inv(tr(T).T).A']
  cond_name_creuse = ['A','tr(T).T','inv(tr(T).T).A']
  cond = [condA]
  pre_cond = [pre_cond1]
  cond_1 = [cond1]

  barWidth = 0.01

  fig = plt.subplots(figsize = (12,8))

  br1 = np.arange(len(cond))
  br2 = [x + barWidth for x in br1]
  br3 = [x + barWidth for x in br2]

  plt.bar(br1,condA,color='r', width = barWidth,label = "Conditionnement de A")
  plt.bar(br2,cond1,color='g', width = barWidth, label = "Conditionnement de tr(T)*-1.T*-1")
  plt.bar(br3,pre_cond1,color='b', width = barWidth, label = "Conditionnement de M**-1.A")
  plt.xticks([r + barWidth for r in range(len(cond))],["Matrice symétrique défini positive"])
  plt.yscale('log')

  plt.legend()
  plt.show()

#Question 5

def plot_cond_creuse(B):
  B1 = factorisation_incomplete_cholesky(B)

  condB = np.linalg.cond(B)
  cond2 = np.linalg.cond(np.dot(B1,B1.T))
  pre_cond2 = np.linalg.cond(np.dot(np.linalg.inv(np.dot(B1,B1.T)),B))

  cond_name_creuse = ['A','tr(T).T','inv(tr(T).T).A']
  cond = [condB]
  pre_cond = [pre_cond2]
  cond_1 = [cond2]

  barWidth = 0.01

  fig = plt.subplots(figsize = (12,8))

  br1 = np.arange(len(cond))
  br2 = [x + barWidth for x in br1]
  br3 = [x + barWidth for x in br2]

  plt.bar(br1,cond,color='r', width = barWidth,label = "Conditionnement de A")
  plt.bar(br2,cond_1,color='g', width = barWidth, label = "Conditionnement de tr(T)*-1.T*-1")
  plt.bar(br3,pre_cond,color='b', width = barWidth, label = "Conditionnement de M**-1.A")
  plt.xticks([r + barWidth for r in range(len(cond))],["Matrice symétrique défini positive creuse"])
  plt.yscale('log')

  plt.legend()
  plt.show()

fig,ax1 = plt.subplots(2)

def list_cond():
  L = [definite_matrix(4,0) for i in range(10)]
  L_ch =[factorisation_cholesky(L[i]) for i in range (10)]
  A1 = [np.linalg.cond(L[i]) for i in range (10)]
  A2 = [np.linalg.cond(np.dot(L_ch[i],L_ch[i].T)) for i in range (10)]
  A3 = [np.linalg.cond(np.dot(np.linalg.inv(np.dot(L_ch[i],L_ch[i].T)),L[i])) for i in range (10)]
  return A1,A2,A3

def list_cond_creuse():
  L = [definite_matrix(4,2) for i in range(10)]
  L_ch =[factorisation_incomplete_cholesky(L[i]) for i in range (10)]
  A1 = [np.linalg.cond(L[i]) for i in range (10)]
  A2 = [np.linalg.cond(np.dot(L_ch[i],L_ch[i].T)) for i in range (10)]
  A3 = [np.linalg.cond(np.dot(np.linalg.inv(np.dot(L_ch[i],L_ch[i].T)),L[i])) for i in range (10)]
  return A1,A2,A3

def plot_list_cond():
  x = np.linspace(0,10,10)
  y1 = list_cond()[0]
  y2 = list_cond()[1]
  y3 = list_cond()[2]

  ax1[0].plot(x,y1,'b--',label="Conditionnement de A")
  ax1[0].plot(x,y2,'r', label = "Conditionnement de tr(T)*-1.T*-1")
  ax1[0].plot(x,y3,'g*', label = "Conditionnement de tr(T)*-1.T*-1.A")
  ax1[0].set_title('Matrice définie positive')
  ax1[0].legend()

## matrice creuse
  x1 = np.linspace(0,10,10)
  y_1 = list_cond_creuse()[0]
  y_2 = list_cond_creuse()[1]
  y_3 = list_cond_creuse()[2]

  ax1[1].plot(x,y_1,'b--',label="Conditionnement de A")
  ax1[1].plot(x,y_2,'r', label = "Conditionnement de tr(T)*-1.T*-1")
  ax1[1].plot(x,y_3,'g*', label = "Conditionnement de tr(T)*-1.T*-1.A")
  ax1[1].set_title('Matrice creuse définie positive')
  ax1[1].legend()
  plt.show()

A2 = np.array([[1,1,1,1],[1,5,5,5],[1,5,14,14],[1,5,14,15]])
A21 = factorisation_cholesky(A2)
#B = generate_matrix_sym(4)
#B1 = matrice_creuse(B,4,6)
#A1 = matrice_creuse(A2,4,10)
#T = factorisation_incomplete_cholesky(B1)



if __name__ == "__main__":
  plot_list_cond()