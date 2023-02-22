from cProfile import label
from re import M
import numpy as np
from scipy.special import legendre 
import matplotlib.pyplot as plt

Nm=1000
Xm = np.linspace(-999/1000, 999/1000, Nm)
N = 9
X = np.linspace(-99/100, 99/100, N)
Jacob = np.zeros((N, N))
for i in range(0, N):
    for j in range (0, N):
        if (i != j): 
            Jacob[i, j] = -1/((X[i]-X[j])**2)
        else :
            sum1 = 0
            sum2 = 0
            for k in range(0, N) :
                if (k != i):
                    sum2+=1/((X[i]-X[k]))**2
            Jacob[i, j] = -1/((X[i]+1)**2)-1/((X[i]-1)**2)-sum2

def Computef(U) :
    f = np.empty(N)
    for i in range(0, N):
        sum1 = 0
        for k in range(0, N) :
            if (k != i):
                sum1+=1/(U[i]-U[k])    
        f[i] = 1/(U[i]+1)+1/(U[i]-1)+sum1
    return f
                

def newton_raphson(f,J,U0,Nmax,epsilon):
    U = U0
    Y = f(U)    
    k = 0
    while(np.linalg.norm(Y) > epsilon and k < Nmax):
        H = Jacob
        V = np.linalg.lstsq(H,-Y,rcond=None)[0]
        alpha = 1
        
        while(np.linalg.norm(f(U + alpha * V)) > np.linalg.norm(Y)):
            alpha /= 2
            if (alpha < 10**(-3)):
                print(k)    
                print("Error : Algorithm Failed")
        U = U + alpha * V
        Y = f(U)
        k += 1

    return U

axe = np.zeros((1,N))
U = newton_raphson(Computef, Jacob, X, 19,10**(-1000))
plt.plot(U, axe[0], 'o', label="roots")

#########################################################
legen2 = legendre(2)

dlegen2 = legen2.deriv()
vdlegen2 = np.zeros((1,Nm))
for i in range (0,Nm): 
    vdlegen2[0, i] = dlegen2(Xm[i])
roots2 = np.roots(dlegen2)
plt.plot(Xm, vdlegen2[0], label="p = 2")
########################################################################""
legen3 = legendre(3)

dlegen3 = legen3.deriv()
vdlegen = np.zeros((1,N))
roots3 = np.roots(dlegen3)

vdlegen3 = np.zeros((1,Nm))
for i in range (0,Nm): 
    vdlegen3[0, i] = dlegen3(Xm[i])
plt.plot(Xm, vdlegen3[0], label="p = 3")
#############################################################################

legen4 = legendre(4)
dlegen4 = legen4.deriv()
vdlegen = np.zeros((1,N))
roots4 = np.roots(dlegen4)
vdlegen4 = np.zeros((1,Nm))
for i in range (0,Nm): 
    vdlegen4[0, i] = dlegen4(Xm[i])
plt.plot(Xm, vdlegen4[0], label="p = 4")

##############################################################################
legen5 = legendre(5)
dlegen5 = legen5.deriv()
vdlegen = np.zeros((1,N))
roots5 = np.roots(dlegen5)
vdlegen5 = np.zeros((1,Nm))
for i in range (0,Nm): 
    vdlegen5[0, i] = dlegen5(Xm[i])
plt.plot(Xm, vdlegen5[0], label="p = 5")

plt.legend()
plt.show()
