import matplotlib as plt
import numpy as np

airDensity = 1.292 

##the curve representing the upper surface of the wing. Useful in order to calculate the lenhth of the curve for the pressure calculation
def f(x):
    return x**3 + x**2



##For a fixed l in [0,1], this equation defines a curve situated between the upper side of the wing and the maximal altitude beyond which 
# the air is not disturbed by the wing (in [3h_min, 3h_max]).
#Useful in order to calculate the lenhth of the curve for the pressure calculation
def f_l(x, l, f, h_max):
    return (1-l)*f(x)+3*h_max


#y is the height
#h_min (resp. h_max) the minumum (resp. maximum) of the airfoil
def P(y, h_min, h_max, curvesLength, time):
    V = curvesLength/time
    P_d = 0.5 * airDensity * V**2
    return P_d


#function f need to be fixed
#number h need to be fixed
if __name__ == "__main__"


