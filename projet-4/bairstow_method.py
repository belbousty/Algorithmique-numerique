import numpy as np
import matplotlib.pyplot as plt
from math import *
import random as rand
import time
import cmath
from numpy.polynomial import Polynomial


def div_eucli(Q, Z):
    q, s = np.polydiv(Q, Z)
    if len(s)== 2:
        return q, s[0], s[1]
    else :
        return q, 0, s[0]



def jacob(P, U):
    """
    returns the Jacobienne matrix of Euclidean division
    """
    Z = np.array([1, U[0], U[1]])
    (q, R, S)=div_eucli(P, Z)
    q1, R1, S1=div_eucli(q, Z)
    l1 = [U[0]*R1-S1, -R1]
    l2 = [U[1]*R1, -S1]
    ja = [l1, l2]
    return ja


def racines( B, C):
    delta = B*B - 4 * C
    if delta > 0:
        return [(-B + np.sqrt(delta)) / 2, (-B - np.sqrt(delta)) / 2]
    elif delta == 0:
        return [-B / 2 , -B / 2]
    else:
#            sol = np.array( [-B + np.sqrt(-delta)j) / 2 ),( -B - np.sqrt(-delta)j ) / 2 ) ])
        return [0,0]


def fu(P, U):
    """
    returns the rest of euclidean division
    """
    a = np.array([1, U[0], U[1]])
    q,r,s = div_eucli(P, a)
    return [r,s]

def subtraction(X,Y):
    W = np.zeros(len(X))
    for i in range(len(X)):
        W[i] = X[i]-Y[i]
    return W


def newton_raphson(f, P, J,U0,Nmax,epsilon):
    U = []
    U.append(U0[0])
    U.append(U0[1])
    Y = f(P, U)
    k = 0
    while(np.linalg.norm(Y) > epsilon and k < Nmax):
        V = np.linalg.lstsq(J(P,U),Y,rcond=None)[0]
        U = subtraction(U,V)
        Y = f(P, U)
        k += 1
    return U

def unicité(L):
    V = []
    for e in L:
        if e not in V:
            V.append(e)
    return V        

def bairstow(P, b, c, N, eps):
    root=[]
    l= len(P)
    U0=[b, c]
    U=U0
    q = P
    for i in range (l//2):
        U = np.array(newton_raphson(fu, q, jacob, U0, N, eps))
        (q, r, s) = div_eucli(q, [1, U[0], U[1]]) 
        root.append(racines( U[0], U[1])[0])
        root.append(racines( U[0], U[1])[1])
    if l%2 == 1 :
        root.append(racines( U[0], U[1])[0])
        root.append(racines( U[0], U[1])[1])
    return unicité(root)

def poly_gen(n):
    s=[]
    for i in range(n):
        s.append(rand.randint(0,20))
    return s

def complexite_numpy(P, n) :
    p = list(reversed(P))
    debut = time.time()
    np.roots(P)
    fin = time.time()
    return (fin - debut)

def complexite_bairstow(P, n):
    b = rand.randint(0,10)
    c = rand.randint(0,10)
    debut = time.time()
    bairstow(P, b, c, 100, 0.01)
    fin = time.time()
    return (fin - debut)

def complexite():
    x=[]
    y=[]
    z=[]
    for i in range(1, 30) :
        p = poly_gen(i)
        x.append(complexite_bairstow(p, i))
        y.append(complexite_numpy(p, i))
        z.append(i)
    p2 = plt.plot(z, x, label="bairstow")
    p3 = plt.plot(z, y, label="numpy")
    plt.legend(loc="best", shadow=True)
    plt.xlabel("degre n of the polynomial")
    plt.ylabel("Time in ms")
    plt.title("Time execution between numpy and bairstow algorithme")
    plt.show()

print(Polynomial([42,-83,53,-13,1]).roots())
print(bairstow([1,-13,53,-83,42],0.5,0.5,100,0.01))
#print(newton_raphson(fu, [1,-13,53,-83,42], jacob, [0.5,0.5], 100, 0.01))
#print(racines(5,1))
complexite()