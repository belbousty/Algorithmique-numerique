import numpy as np
import matplotlib.pyplot as plt
from newton_raphson import *
from math import *
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox

import matplotlib.pyplot as plt
import matplotlib.image as mpimg



def Jacobienne(f,U) :
    t = 10**(-10)
    J = np.zeros((len(f(U)),len(U)))
    for i in range(len(f(U))):
        
        for j in range(len(U)):
            L = [0 for i in range(len(U))]
            L[j] = t
            J[i][j] = (f([U[k] + L[k] for k in range(len(U))])[i] - f(U)[i]) / t
    return J



elastric_force = lambda k,x : -k*x

centrifigal_force = lambda k,X,X0 : np.array([-k*(X[0] - X0[0]),-k*(X[1] - X0[1])])

def gravitational_force(k,X0):
    return lambda X: np.array([-k*(X[0] - X0[0])/((X[0] - X0[0])**2 + (X[1] - X0[1])**2)**(3/2)
                                            ,-k*(X[1] - X0[1])/((X[0] - X0[0])**2 + (X[1] - X0[1])**2)**(3/2)]) 


G = 6.67e-11
Me = 5.97e24
Ms = 1.9892e30
Mj = 1.898e27
k0 = G*Ms
k1 = G*Me
k0_r = 1
k1_r = k0/k1
k2 = G*Mj
Rs = 6.9639e8
Re = 6.4e6
Re_r = Rs/Re
Rs_s = 1 - Re_r
Rj = 6.9911e7
w = k0_r + k1_r
Xx0 = np.array([0,0])
Xx1 = np.array([1,0])
Xx2 = np.array([k1_r/(k0_r+k1_r),0])

def distance(X,Y):
    return sqrt((X[0] - Y[0])**2 + (X[1] - Y[1])**2)

def Lagragian_pts(K0,X0,K1,X1,K2,X2):
    return lambda X: np.array([ K2*(X[0] - X2[0]) - K0*(X[0] - X0[0])/(distance(X,X0)**3) - K1*(X[0] - X1[0])/(distance(X,X1))**(3)
                                            ,K2*(X[1] - X2[1]) - K0*(X[1] - X0[1])/(distance(X,X0)**3) - K1*(X[1] - X1[1])/(distance(X,X1))**(3)]) 



L1 = newton_raphson(Lagragian_pts(1,np.array([0,0]),0.01,np.array([1,0]),1,np.array([0.01,0])),Jacobienne,np.array([1.5,0]),100,10e-5)


L2 = newton_raphson(Lagragian_pts(1,np.array([0,0]),0.01,np.array([1,0]),1,np.array([0.01,0])),Jacobienne,np.array([0.5,0]),100,10e-5)


L3 = newton_raphson(Lagragian_pts(1,np.array([0,0]),0.01,np.array([1,0]),1,np.array([0.01,0])),Jacobienne,np.array([-1,0]),100,10e-5)


L4 = newton_raphson(Lagragian_pts(1,np.array([0,0]),0.01,np.array([1,0]),1,np.array([0.01,0])),Jacobienne,np.array([0.8,0.5]),100,10e-5)

L5 = newton_raphson(Lagragian_pts(1,np.array([0,0]),0.01,np.array([1,0]),1,np.array([0.01,0])),Jacobienne,np.array([0.8,-0.5]),100,10e-5)

X = [L1[0],L2[0],L3[0],L4[0],L5[0]]
Y = [L1[1],L2[1],L3[1],L4[1],L5[1]]


ax = plt.subplots()[1]


ax.set_xlim(-2, 2)
ax.set_ylim(-1, 1)

soleil = mpimg.imread('soleil.png')
jupiter = mpimg.imread('jupiter.png')
satellite = mpimg.imread('satellite.png')

imageSo = OffsetImage(soleil, zoom=0.06)
imageJu = OffsetImage(jupiter, zoom=0.05)
imageSa = OffsetImage(satellite, zoom=0.04)

J = AnnotationBbox(imageJu, (1, 0))
so = AnnotationBbox(imageSo, (0,0))
S1 = AnnotationBbox(imageSa,(L1[0]+0.01,L1[1]))
S2 = AnnotationBbox(imageSa,(L2[0]-0.02,L2[1]))
S3 = AnnotationBbox(imageSa,(L3[0],L3[1]))
S4 = AnnotationBbox(imageSa,(L4[0],L4[1]))
S5 = AnnotationBbox(imageSa,(L5[0],L5[1]))

ax.add_artist(J)
ax.add_artist(so)
ax.add_artist(S1)
ax.add_artist(S2)
ax.add_artist(S3)
ax.add_artist(S4)
ax.add_artist(S5)


plt.draw()
plt.savefig('add_picture_matplotlib_figure.png',bbox_inches='tight')
plt.show()
