import numpy as np

def step_RK4( y, t, h, f):
    k1 = f(y,t)

    k2 = f([x + h/2*k for x,k in zip(y,k1) ], t+h/2) 
    k3 = f([x + h/2*k for x,k in zip(y,k2) ], t+h/2)
    k4 = f([x + h/2*k for x,k in zip(y,k3) ], t)
    y = [y+(h/6)*(k1+ 2*k2 +2*k3 + k4) for y,k1,k2,k3,k4 in zip(y,k1,k2,k3,k4)]
    return y

def meth_n_step(y0, t0, N, h, f, meth="RK4"): 
    Y = [y0]
    for i in range (0, N-1):
        Y.append(step_RK4( Y[i], t0+h*i, h, f))
    return [y[0] for y in Y]

 
    
    
