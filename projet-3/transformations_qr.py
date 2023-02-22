
import numpy as np
import matplotlib.pyplot as plt
import time
import matplotlib.patches as mpatches
import random
from numba import jit,prange
from Householder import Householder


@jit(nopython=True,fastmath = True)
def test_bidiag(A):
    (L,l) = np.shape(A)
    for i in range (0,L):
        for j in range (0,l):
            if (i != j and i + 1 != j):
                if (A[i,j] > 0.00001):
                    return False
    return True


@jit(nopython=True,fastmath = True)
def err_rel(theo,appr):
    if (theo == 0):
        return 0
    return abs(appr - theo)/abs(theo)

# def decomp_qr(NMax,BD):
#     n = np.shape(BD)[0]
#     U = np.eye(n)
#     V = np.eye(n)
#     S = BD
#     for i in range (0, NMax):
#         (Q1, R1) = decomp_qr(np.transpose(S))
#         (Q2, R2) = decomp_qr(np.transpose(R1))
#         S = R2
#         U = np.dot(U,Q2)
#         V = np.dot(np.transpose(Q1),V)
#     return (U,S,V)


@jit(nopython=True,fastmath = True)
def eigenvalue(A):
    (L,l) = np.shape(A)
    x = np.zeros((l,1))
    for i in range (0,l):
        x[i,0] = 1
    for j in range(20):
        x = np.dot(A,x)/(np.linalg.norm(np.dot(A,x)))
    l = np.dot(np.dot(np.transpose(x),A),x)/np.linalg.norm(x)
    return l[0,0]




@jit(nopython=True,fastmath = True)
def decomp_qr_preli(NMax,BD):
    n = np.shape(BD)[0]
    U = np.eye(n)
    V = np.eye(n)
    S = BD
    for i in range (0, NMax):
        (Q1, R1) = np.linalg.qr(np.transpose(S))
        (Q2, R2) = np.linalg.qr(np.transpose(R1))
        S = R2
        # Test if the matrixs R1 and R2 are bidiagonal
        if (not test_bidiag(R1) or not test_bidiag(R2)):
            print("Not Bidiag")
        U = np.dot(U,Q2)
        V = np.dot(np.transpose(Q1),V)
    return (U,S,V)


@jit(nopython=True,fastmath = True)
def diag(n):
    BD = np.zeros((n,n))
    for i in range (n):
        BD[i,i] = random.randrange(0,100)
        if (i + 1 < n):
            BD[i,i + 1] = random.randrange(0,100)
    return BD


@jit(nopython=True,fastmath = True)
def vect(n):
    vect = np.zeros((1,n))
    for i in range (n):
        vect[1,i] = random.randrange(0,100)
    return vect


@jit(nopython=True,fastmath = True)
def column_matrix(A,val):
    (L,l) = np.shape(A)
    vect = np.zeros((1,L))
    for i in range (L):
        vect[1,i] = A[val,i]
    return vect


@jit(nopython=True,fastmath = True)
def householder_qr(NMax,A):
    (L,l) = np.shape(A)
    res = np.copy(A)
    for k in range (0, NMax):
        vect_1 = column_matrix(A,k)
        vect_2 = np.zeros((1,L))
        Q = Householder(vect_1, vect_2)
        res = np.dot(Q,res)
    return res


@jit(nopython=True,fastmath = True)
def cmp_decomp_qr_preli(NMax,BD):
    (L,l) = np.shape(BD)
    
    x = []

    y_nump = []
    y_hh = []

    n = np.shape(BD)[0]
    U_nump = np.eye(n)
    V_nump = np.eye(n)
    S_nump = np.copy(BD)

    U_hh = np.eye(n)
    V_hh = np.eye(n)
    S_hh = np.copy(BD)

    for k in range (NMax):
        (Q1_nump, R1_nump) = np.linalg.qr(np.transpose(S_nump))
        (Q2_nump, R2_nump) = np.linalg.qr(np.transpose(R1_nump))
        S_nump = np.copy(R2_nump)
        # Test if the matrixs R1 and R2 are bidiagonal
        if (not test_bidiag(R1_nump) or not test_bidiag(R2_nump)):
            print("Not Bidiag")
        U_nump = np.dot(U_nump,Q2_nump)
        V_nump = np.dot(np.transpose(Q1_nump),V_nump)
        res_nump = np.dot(U_nump,S_nump)
        res_nump = np.dot(res_nump,V_nump)
        sum_err_rel_nump = 0
        for i in range(L):
            for j in range(l):
                sum_err_rel_nump = sum_err_rel_nump + err_rel(BD[i,j],res_nump[i,j])
        sum_err_rel_nump = sum_err_rel_nump/(L*l)



        (Q1_hh, R1_hh) = householder_qr(k,np.transpose(S_hh))
        (Q2_hh, R2_hh) = householder_qr(k,np.transpose(R1_hh))
        S_hh = np.copy(R2_hh)
        # Test if the matrixs R1 and R2 are bidiagonal
        if (not test_bidiag(R1_hh) or not test_bidiag(R2_hh)):
            print("Not Bidiag")
        U_hh = np.dot(U_hh,Q2_hh)
        V_hh = np.dot(np.transpose(Q1_hh),V_hh)
        res_hh = np.dot(U_hh,S_hh)
        res_hh = np.dot(res_hh,V_hh)
        sum_err_rel_hh = 0
        for i in range(L):
            for j in range(l):
                sum_err_rel_hh = sum_err_rel_hh + err_rel(BD[i,j],res_hh[i,j])
        sum_err_rel_hh = sum_err_rel_hh/(L*l)


        x.append(k)
        y_nump.append(sum_err_rel_nump)
        y_hh.append(sum_err_rel_hh)
    return (x,y_nump,y_hh)
    


@jit(nopython=True,fastmath = True)
def cmp_decomp_qr(NMax,BD):
    (L,l) = np.shape(BD)
    x = []

    y_nump = []
    y_hh = []
    
    n = np.shape(BD)[0]

    U_nump = np.eye(n)
    V_nump = np.eye(n)
    S_nump = np.copy(BD)

    U_hh = np.eye(n)
    V_hh = np.eye(n)
    S_hh = np.copy(BD)

    for k in range (1,NMax):
        (Q1_nump, R1_nump) = np.linalg.qr(np.transpose(S_nump))
        (Q2_nump, R2_nump) = np.linalg.qr(np.transpose(R1_nump))
        S_nump = np.copy(R2_nump)
        # Test if the matrixs R1 and R2 are bidiagonal
        if (not test_bidiag(R1_nump) or not test_bidiag(R2_nump)):
            print("Not Bidiag")
        U_nump = np.dot(U_nump,Q2_nump)
        V_nump = np.dot(np.transpose(Q1_nump),V_nump)
        res_nump = np.dot(U_nump,V_nump)
        sum_err_rel_nump = 0
        for i in range(L):
            for j in range(l):
                sum_err_rel_nump = sum_err_rel_nump + abs(res_nump[i,j])
        sum_err_rel_nump = sum_err_rel_nump/(L*l)


        (Q1_hh, R1_hh) = householder_qr(k,np.transpose(S_hh))
        (Q2_hh, R2_hh) = householder_qr(k,np.transpose(R1_hh))
        S_hh = np.copy(R2_hh)
        # Test if the matrixs R1 and R2 are bidiagonal
        if (not test_bidiag(R1_hh) or not test_bidiag(R2_hh)):
            print("Not Bidiag")
        U_hh = np.dot(U_hh,Q2_hh)
        V_hh = np.dot(np.transpose(Q1_hh),V_hh)
        res_hh = np.dot(U_hh,V_hh)
        sum_err_rel_hh = 0
        for i in range(L):
            for j in range(l):
                sum_err_rel_hh = sum_err_rel_hh + abs(res_hh[i,j])
        sum_err_rel_hh = sum_err_rel_hh/(L*l)


        x.append(k)
        y_hh.append(sum_err_rel_hh)
        y_hh.append(sum_err_rel_hh)
    return (x,y_nump,y_hh)
    
    
def plot_test_1(x,y1,y2):
    plt.plot(x,y1,y2,color='b',label='Err_rel')
    plt.xlabel('Nombre ditérations',fontsize=20)
    plt.ylabel('Erreur relative de U*S*V = BD',fontsize=20)
    plt.tick_params(axis='x', labelsize=20)
    plt.tick_params(axis='y', labelsize=20)
    blue_patch = mpatches.Patch(color='blue', label='Numpy')
    red_patch = mpatches.Patch(color='red', label='Householder')
    plt.legend(handles=[blue_patch,red_patch],fontsize=20)
    plt.yscale('log')
    plt.show()


def plot_test_2(x,y1,y2):
    plt.plot(x,y1,y2,color='b',label='Err_rel')
    plt.xlabel('Nombre ditérations',fontsize=20)
    plt.ylabel('Erreur absolue du produit U*V',fontsize=20)
    plt.tick_params(axis='x', labelsize=20)
    plt.tick_params(axis='y', labelsize=20)
    blue_patch = mpatches.Patch(color='blue', label='Numpy')
    red_patch = mpatches.Patch(color='red', label='Householder')
    plt.legend(handles=[blue_patch,red_patch],fontsize=20)
    plt.yscale('log')
    plt.show()


BD = diag(300)

(x,y1,y2) = cmp_decomp_qr_preli(1000,BD)
plot_test_1(x,y1,y2)

# (x,y1,y2) = cmp_decomp_qr(1000,BD)
# plot_test_2(x,y1,y2)

# A = diag(4)
# l = eigenvalue(A)

# print(l)
# print(A)