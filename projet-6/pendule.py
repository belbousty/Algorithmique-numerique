import numpy as np
import matplotlib.pyplot as plt
from inspect import getframeinfo
from scipy.integrate import odeint 


g = 9.81
L=1
M=1


# PENDULE SIMPLE

angles = linspace(-pi/2,pi/2,100)
temps = linspace(0.5*pi, 10, 100000)


def eq_pendule_simple(u,t):
    return (u[1],-g/L*sin(u[0]))


def calcul_freq(y):
     L = []
     if y[0] >= 0:
          signe = 1
     else:
          signe = -1
     for k in range(len(y)):
          if signe*y[k] < 0:
               L.append(t[k])
               signe = signe*(-1)
          if len(L) == 3:
               break
     return 1/(L[2]-L[0])  # f = 1/T
               
          
def obtenir_freq(i):
    courbe = odeint(eq_pendule_simple,[angles[i],0],t)
    y = courbe[:,0]
    return calcul_freq(y)


courbe = odeint(eq_pendule_simple,(pi/2,0),t)
plt.plot(t, courbe[:,0])
plt.xlabel('temps(s)')
plt.ylabel(''r'$\theta$(rad)')
plt.show()


F = []

for i in range(len(angles)):
     F.append(obtenir_freq(i)*2*pi)  #Je transforme les fréquences en pulsation

plt.plot(angles,F)
plt.xlabel(''r'$\theta$(rad)')
plt.ylabel('Pulsation propre(rad/s)')
plt.show()




# PENDULE DOUBLE

def fonctions(init,temps):
    angle1,angle2,dangle1dt,dangle2dt = init
    ddangle1ddt = (-g*(2*M+M)*np.sin(angle1) - M*g*np.sin(angle1-2*angle2) -2*np.sin(angle1-angle2)*M*((dangle2dt**2)*L+(dangle1dt**2)*L*np.cos(angle1-angle2)))/(L*(2*M+M-M*np.cos(2*angle1-2*angle2)))
    ddangle2ddt =(2*np.sin(angle1-angle2)*((dangle1dt**2)*L*(M+M)+g*(M+M)*np.cos(angle1) + (dangle2dt**2)*L*M*np.cos(angle1-angle2)))/(L*(2*M+M-M*np.cos(2*angle1-2*angle2)))
    return dangle1dt, dangle2dt, ddangle1ddt, ddangle2ddt

temps  = np.linspace(0,200,10000)


def calcul(angles_init,temps):
    return scipy.integrate.odeint(fonctions, angles_init, temps)



def plot_position():
    angles_init = [np.pi/2, np.pi/2, 0, 0]

    valeurs = calcul(angles_init,temps)

    x1 = L*np.sin(valeurs[:,0])   
    y1 = -L*np.cos(valeurs[:,0])   

    x2 = x1 + L*np.sin(valeurs[:,1])
    y2 = y1 - L*np.cos(valeurs[:,1])
    
    plt.plot(x2,y2)
    plt.xlabel('x en mètres')
    plt.ylabel('y en mètres')
    plt.show()

plot_position()





