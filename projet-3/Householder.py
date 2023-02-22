import numpy as np
from numba import jit,prange


#Takes in argument two vectors of the same size (LIST)
#Returns the Householder matrix of the vectors given in entry
@jit(nopython=True,fastmath = True)
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

#Non optimized product:
#Takes in argument a Householder Matrix(ARRAY) and a vector(LIST)
#Returns their product (COLUMN VECTOR)

def HouseholderVProduct(H,Z):
    n = np.size(X)
    Xc=np.array(X).reshape(-1,1)
    Y = np.array([0]*n)
    for i in range(n):
        s=0
        for j in range(n):
            s+=H[i,j]*X[j]
        Y[i]=s
    return np.array(Y).reshape(-1,1)

#Takes in argument a Householder matrix H
#Returns the vector U such as : H =  Id - 2*U*tU

def HouseholderextractU(H):
    n = np.size(H[0])
    I = np.eye(n)
    W = (I-H)/2
    U = np.zeros(n)
    U[0] = np.sqrt(W[0,0])
    for k in range(1,n):
        if W[k,0]>=0:
            U[k] = np.sqrt(W[k,k])
        else:
            U[k] = - np.sqrt(W[k,k])
    return U.reshape(-1,1)

#Takes in argument a Householder matrix, it's vectors (X,Y) and a vector
#Returns the product H*Z (With less operations)

def HouseholderVOptimum(X,Y,H,Z):
    Zc = np.array(Z).reshape(-1,1)
    U = np.array(X)-np.array(Y)
    n = np.size(Zc)
    norm =0
    for i in range(n):
        norm+= U[i]**2
    norm=np.sqrt(norm)
    U = U/norm
    cte = 0
    for i in range(n):
        cte+=Zc[i][0]*U[i]
    output = Zc-2*cte*(np.array(U).reshape(-1,1))
    return output

#Takes in argument a Householder matrix and a square matrix (H,M)
#Returns the product :  H*M

def HouseholderMOptimumLeft(X,Y,H,M):
    n = np.size(H[0])
    O = np.eye(n)
    for i in range(n):
        V = HouseholderVOptimum(X,Y,H, np.array(M[:,i]))
        O[:,i] = np.array(V).reshape(1,-1)
    return O


# A = Householder([3,4,0], [0,0,5])
# B = Householder([7,1,1], [12,0,-9])
# print(A)
# print(B)
# print(HouseholderMOptimumLeft(A,B))


#Takes in argument a Householder matrix and a square matrix (H,M)
#Returns the product : M*H

def HouseholderMOptimumRight(X,Y,H,M):
    Out = HouseholderMOptimumLeft(X,Y,np.array(H).T, np.array(M).T)
    return np.array(Out).T


X = [3,4,0]
Y = [0,0,5]
H = Householder(X,Y)
print(HouseholderMOptimumRight(X,Y,H,np.eye(np.size(H[0]))))

