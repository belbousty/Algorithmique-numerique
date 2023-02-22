from array import array
from cProfile import label
from calendar import c
from curses import ERR
from re import X
import numpy as np 
import matplotlib.pyplot as plt
import RK4 
from scipy.integrate import odeint

###
X = lambda X,t : [X[1], -X[0]]
Y = lambda Y,t : [Y[1], -Y[0]]
####
x0 = [1,0]
y0 = [0,1]
t0 = 0
h = 10**-2
N = 1000
Tp = np.linspace(0,N*h,1000)
####
def show_X_Y_solution_for_2sys() :
    Xp = RK4.meth_n_step(x0, t0, N, h, X, "RK4")
    Yp = RK4.meth_n_step(y0, t0, N, h, Y, "RK4")
    plt.plot(Tp, Xp, label="x(t)")
    plt.plot(Tp, Yp, label="y(t)")
    plt.legend()
    plt.grid()
    plt.show()
    plt.plot(Xp, Yp)
    plt.legend()
    plt.grid()
    plt.show()

def show_X_Y_solutions_for_3sys():
    X_c = lambda x,t : [x[1], np.cos(t)-2*x[0]]
    Y_c = lambda y,t : [y[1], np.cos(t)-2*y[0]]
    ####
    xc0 = [1,2**0.5]
    yc0 = [0,1]
    tc0 = 0
    hc = 10**-1
    Nc = 1000
    Tpc = np.linspace(0,N*h,1000)
    Xp = RK4.meth_n_step(xc0, tc0, Nc, hc, X_c, "RK4")
    Yp = RK4.meth_n_step(yc0, tc0, Nc, hc, Y_c, "RK4")
    plt.plot(Tpc, Xp, label="x(t)")
    plt.plot(Tpc, Yp, label="y(t)")
    plt.legend()
    plt.grid()
    plt.show()
    plt.plot(Xp, Yp)
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":
    show_X_Y_solution_for_2sys()
    show_X_Y_solutions_for_3sys()