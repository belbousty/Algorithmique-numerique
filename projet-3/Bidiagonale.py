import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import block_diag

#householder algorithm not optimized
def Householder(X,Y):
    Xc = np.array(X).reshape(-1,1)
    Yc = np.array(Y).reshape(-1,1)
    n = np.size(Xc)
    I = np.eye(n)
    U = Xc - Yc
    norm = 0
    for i in range(n):
        norm += (Xc[i]-Yc[i])**2
    norm = np.sqrt(norm)
    U = U/norm
    H = I - 2*( U @ np.transpose(U))
    return H

def Bidiagonal_Form(A, n): #  this function returns the bidiagonal form of A
    Q_left = np.eye(n)
    Q_right = np.eye(n)
    BD = A
    for i in range (0,n-2):
        Vl = np.zeros((n-i,1))
        print()
        Vl[0,0] = np.linalg.norm(BD[i:,i])
        Ql = Householder(BD[i:n,i], Vl)
        if (i!= 0): 
            Ql = block_diag(np.eye(i), Ql)
        Q_left = np.dot(Q_left,Ql)
        BD = np.dot(Ql,BD)
        if (i < n-2):
            Vr = np.zeros((1,n-i-1))
            Vr[0,0] = np.linalg.norm(BD[i,i+1:])
            Qr = Householder(BD[i, i+1:n], Vr)
            Qr = block_diag(np.eye(n-np.size(Qr[0])), Qr)
            Q_right  = np.dot(Qr,Q_right)
            BD = np.dot(BD,Qr)
    return (BD, Q_right, Q_left)

A = np.array([[4,1,-2, 2]
            ,[1,2,0,1]
            ,[-2,0,3,-2]
            ,[2,1,-2,-1]])
print(Bidiagonal_Form(A, 4))

Bidiagonal_Form_A = Bidiagonal_Form(A, 4)[0]
B =  np.array([[5, 4.5, 0, 0],
            [0, 1.76, 1, 0], 
            [0, 0, -1.6, 1.03],
            [0, 0, 0, 2]])

# B is the result of the computed theorical bidiagonal form 
x= np.linspace(0,10,4**2)
erreur_A = np.zeros((1,4**2))
for i in range (0,4):
  for j in range (0,4):
    if (B[i,j] != 0):
      erreur_A[0,i*4+j] = (abs(Bidiagonal_Form_A[i,j]-(B[i,j]))/abs(B[i,j]))
plt.plot(x, erreur_A[0],label="Erreur relative")
plt.show()

        