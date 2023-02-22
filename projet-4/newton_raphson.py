import numpy as np
import matplotlib.pyplot as plt
from math import exp
from math import sqrt
from math import pi

def newton_raphson(f,J,U0,Nmax,epsilon):
    U = U0
    Y = f(U)
    k = 0
    while(np.linalg.norm(Y) > epsilon and k < Nmax):
        H = J(f,U)
        V = np.linalg.lstsq(H,-Y,rcond=None)[0]
        alpha = 1

        while(np.linalg.norm(f(U + alpha * V)) > np.linalg.norm(Y)):
            alpha /= 2
            if (alpha < 1e-3):
                print("Error : Algorithm Failed")
                return None
        U = U + alpha * V
        Y = f(U)
        k += 1
    return U

def Jacobienne(f : "function",U : "Vector") :
    """
    returns the Jacobienne matrix of the function applied on U
    """
    t = 10**(-10)
    J = np.zeros((len(f(U)),len(U)))
    for i in range(len(f(U))):
        for j in range(len(U)):
            L = [0 for i in range(len(U))]
            L[j] = t
            J[i][j] = (f([U[k] + L[k] for k in range(len(U))])[i] - f(U)[i]) / t
    return J

f = lambda x : np.array([x[0]**2 - 2])


#print(newton_raphson(f,Jacobienne,np.array([2]),100,10**(-10)))


gaussienne = lambda x : np.array( [ exp(-(x[0]**2)/2)/(sqrt(2*pi)) -0.3 ] )
#print(newton_raphson(gaussienne,Jacobienne,np.array([4]),100,10**(-10)))