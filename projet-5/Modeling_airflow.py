
import numpy as np
import database_to_array as db 
import matplotlib.pyplot as plt
import scipy.interpolate 
import airfoil_refinement 
import length_curves

(f_e, f_i) = airfoil_refinement.get_airfoil_curves("b737a.dat")
(hMax, hMin) = airfoil_refinement.get_Max_Min("b737a.dat")
density = 1.29
LAMBDA = np.linspace(0,1, 100)
X = np.linspace(0,1, 60)

# Compute  the function f_lambda for upper curves
def Calcul_flambda_e(lambdaa):
    def calcul(x):
        return ((1-lambdaa)*f_e(x)+3*lambdaa*hMax)
    return calcul

# Compute  the function f_lambda for lower curves
def Calcul_flambda_i(lambdaa):
    def calcul(x):
        return (((1-lambdaa)*f_i(x)+3*lambdaa*hMin))
    return calcul
######################################################

#Compute Pression
def Lambda_Pression_e(lambdaa):
    speed = length_curves.distance(Calcul_flambda_e(lambdaa), length_curves.trapezoid, 0, 1, 100)
    pression = 0.5*density*(speed)**2
    return pression

def Lambda_Pression_i(lambdaa):
    value = length_curves.distance(Calcul_flambda_i(lambdaa), length_curves.trapezoid , 0, 1, 100) 
    pression = 0.5*density*(value)**2
    return pression
######################################################

PRESSION_MATRIX = np.zeros((100,60))

(dim,ex,ey,ix,iy) = db.load_foil("b737a.dat")
def plot_curves() :
    cs = scipy.interpolate.CubicSpline(ex,ey)
    plt.plot(ex, cs(ex), label="S")
    cs1 = scipy.interpolate.CubicSpline(ix,iy)
    lambdaaa = 0
    for j in range(0,30):
        plt.plot(X ,[Calcul_flambda_e(lambdaaa)(x) for x in X])
        plt.plot(X ,[Calcul_flambda_i(lambdaaa)(x) for x in X])
        lambdaaa += 0.033
    plt.plot(ix, cs1(ix), label="S")
    plt.show()

def plot_map():
    lambdaaa = 0
    for j in range(0, 30):
        lambdaaa += 0.033
        for i in range (0,60):
          PRESSION_MATRIX[50-1-int(Calcul_flambda_e(lambdaaa)(X[i])*100), i] = Lambda_Pression_e(lambdaaa)


    for k in range(0, 40):
        lambdaaa -= 0.033
        for m in range (0,60):  
            PRESSION_MATRIX[60-1-int(Calcul_flambda_i(lambdaaa)(X[m])*100), m] = Lambda_Pression_i(lambdaaa)      
    M = np.zeros((100, 60))
    for i in range(0,100):
        M[i] = PRESSION_MATRIX[100-i-1]
    plt.imshow(M, cmap="hot",interpolation="bilinear", vmin=np.min(PRESSION_MATRIX[np.nonzero(PRESSION_MATRIX)]), vmax=np.max(PRESSION_MATRIX),)
    plt.show()

if __name__ == "__main__":
    print("Projecting Map ...")
    plot_map() 
    print("Projecting Curves ...")
    plot_curves()    
    
    
    
    
