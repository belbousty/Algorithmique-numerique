from tkinter import N
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as sli
import time
import matplotlib.patches as mpatches
from numba import jit,prange


img_full = plt.imread("p3_takeoff_base.png")


@jit(nopython=True,fastmath = True)
def err_rel(theo,appr):
    if (theo < 0.001):
        theo = 0.001
    if (appr < 0.001):
        appr = 0.001
    return abs(appr - theo)/abs(theo)


@jit(nopython=True,fastmath = True)
def black_and_white(img):
    (L,l,n) = np.shape(img)
    res = np.zeros((L,l))
    for i in range (L):
        for j in range (l):
            res[i,j] = (img[i,j,0] + img[i,j,1] + img[i,j,2])/3
    return res

@jit(nopython=True,fastmath = True)
def compress(S, k):
    L = np.shape(S)[0]
    S_bis = np.copy(S)
    for i in range (k, L):
            S_bis[i] = 0
    return S_bis


@jit(nopython=True,fastmath = True)
def decompose_color(img):
    (L,l,n) = np.shape(img)
    img_R = np.zeros((L,l))
    img_G = np.zeros((L,l))
    img_B = np.zeros((L,l))
    for i in range (L):
        for j in range (l):
            img_R[i,j] = img[i,j,0]
            img_G[i,j] = img[i,j,1]
            img_B[i,j] = img[i,j,2]
    return (img_R,img_G,img_B)


@jit(nopython=True,fastmath = True)
def vector_to_diag(S):
    L = np.shape(S)[0]
    S_bis = np.zeros((L,L))
    for i in range (L):
        S_bis[i,i] = S[i]
    return S_bis


@jit(nopython=True,fastmath = True)
def vector_to_diag_3(S):
    (L,l) = np.shape(S)
    S_bis = np.zeros((L,L,3))
    for i in range (0,L):
            S_bis[i,i,0] = S[i,0]
            S_bis[i,i,1] = S[i,1]
            S_bis[i,i,2] = S[i,2]
    print(S_bis)
    return S_bis

@jit(nopython=True,fastmath = True)
def recompose(U,S,V):
    S = vector_to_diag(S)
    A = np.dot(U,S)
    A = np.dot(A,V)
    return A


@jit(nopython=True,fastmath = True)
def img_recompose(R,G,B):
    (L,l) = np.shape(R)
    img = np.zeros((L,l,3))
    for i in range (L):
        for j in range (l):
            img[i,j,0] = R[i,j]
            img[i,j,1] = G[i,j]
            img[i,j,2] = B[i,j]
    return img


@jit(nopython=True,parallel = True,fastmath = True)
def clipping_check(img):
    (L,l,n) = np.shape(img)
    for i in prange (L):
        for j in prange (l):
            for k in prange (n):
                if (img[i,j,k] > 1):
                    img[i,j,k] = 1
                if (img[i,j,k] < 0):
                    img[i,j,k] = 0


@jit(nopython=True,fastmath = True)
def cmp_compress(N,img):
    (R,G,B) = decompose_color(img)
    (U_R, S_R, V_R) = np.linalg.svd(R, full_matrices=False)
    (U_G, S_G, V_G) = np.linalg.svd(G, full_matrices=False)
    (U_B, S_B, V_B) = np.linalg.svd(B, full_matrices=False)
    (L,l,n) = np.shape(img)
    x = []
    y = []
    for ite in range (1,N):
        S_R_cmp = compress(S_R,ite)
        S_G_cmp = compress(S_G,ite)
        S_B_cmp = compress(S_B,ite)

        A_R = recompose(U_R,S_R_cmp,V_R)
        A_G = recompose(U_G,S_G_cmp,V_G)
        A_B = recompose(U_B,S_B_cmp,V_B)

        res = img_recompose(A_R,A_G,A_B)
        sum_err_rel = 0
        for i in range(L):
            for j in range(l):
                for k in range(n):
                    err = err_rel(img[i,j,k],res[i,j,k])
                    sum_err_rel = sum_err_rel + err
        sum_err_rel = sum_err_rel/(L*l*n)
        x.append(ite)
        y.append(sum_err_rel)
    return (x,y)

def plot_show(x,y):
    plt.plot(x,y,color='b',label='Err_rel')
    plt.title('Erreur relative de la compression')
    plt.xlabel('Ordre de compression')
    plt.ylabel('Erreur relative')
    blue_patch = mpatches.Patch(color='blue', label='Numpy')
    plt.legend(handles=[blue_patch])
    plt.yscale('log')
    plt.show()


<<<<<<< HEAD
(x,y) = cmp_compress(400,img_full)
plot_show(x,y)
=======
def compress_img(N,img_full):
    (R,G,B) = decompose_color(img_full)
    (U_R, S_R, V_R) = np.linalg.svd(R, full_matrices=False)
    (U_G, S_G, V_G) = np.linalg.svd(G, full_matrices=False)
    (U_B, S_B, V_B) = np.linalg.svd(B, full_matrices=False)

    S_R_cmp = compress(S_R,N)
    S_G_cmp = compress(S_G,N)
    S_B_cmp = compress(S_B,N)

    A_R = recompose(U_R,S_R_cmp,V_R)
    A_G = recompose(U_G,S_G_cmp,V_G)
    A_B = recompose(U_B,S_B_cmp,V_B)

    res = img_recompose(A_R,A_G,A_B)

    clipping_check(res)

    plt.imshow(res)
    plt.show()

compress_img(100,img_full)


# (x,y) = cmp_compress(400,img_full)
# plot_show(x,y)
>>>>>>> 301ada699f6e8a62f114182705a356897e5e918e
