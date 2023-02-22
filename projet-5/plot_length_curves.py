from length_curves import *
import matplotlib.pyplot as plt

#METHOD IMPLEMENTATION FOR GRAPHICS

#On a pas besoin de toute la matrice pour le plot, parce que si on veut fair eune matrice de taille 100*100
#ça prend énormément de temps de calcul
def Romberg_plot(f, a, b, n, l, k):
    L = []
    h = (b-a)/2
    R = np.zeros((n, n))
    R[0, 0] = (f(a)+f(b))*(b-a)/2
    for i in range(1, n):
        R[i, 0] = R[i-1, 0]*0.5 + h*sommeRomberg(f, a, i, h)
        for j in range(1, i + 1):
            R[i, j] = ((4**j)*R[i, j-1] - R[i - 1, j - 1])/(4**j - 1)
        h /= 2
    return R[l - 1, k - 1]

def distance_Romberg(f, a, b, n, l, k):    
    return Romberg_plot(sqrt_length(f), a, b, n, l, k)


#Return an array with the value of the distance for each iteration n
def distance_iteration(f, method, a, b, n, val):
    L = []
    for i in range(1, n + 1):
        if(method == Romberg_plot):
            if(distance_Romberg(f, a, b, n, n, i) != 0):
                L.append(distance_Romberg(f, a, b, n, n, i))
        else:
            L.append(distance(f, method, a, b, i))

    return L;

def real_value(x, n):
    L = []
    for i in range(n):
        L.append(x)
    return L
#Return an array containing the difference between the real value and the calculated value by a method for each iteration
def distance_difference(f, method, a, b, n, val):
    L = []
    for i in range(1, n + 1):
        L.append(np.abs(distance(f, method, a, b, i) - val))
    return L;

def relative_error(f, method, a, b, n, val):
    L = []
    for i in range(1, n + 1):
        if(method == Romberg_plot):
            if(distance_Romberg(f, a, b, n, n, i) != 0):
                L.append(np.abs((distance_Romberg(f, a, b, n, n, i) - val)/val))
        else:
            L.append(np.abs((distance(f, method, a, b, i) - val)/val))
    return L;

def absolute_error(f, method, a, b, n, val):
    L = []
    for i in range(1, n + 1):
        if(method == Romberg_plot):
            if(distance_Romberg(f, a, b, n, n, i) != 0):
                L.append((distance_Romberg(f, a, b, n, n, i) - val)/val)
        else:
            L.append((distance(f, method, a, b, i) - val)/val)
    return L;

def trace(which_function, val):

    X = np.linspace(1, 1000, 1000)
    X1 = np.linspace(1, 15, 15)
    Y1 = np.array(which_function(function, rectangle, 0, 2, 1000, val))
    Y2 = np.array(which_function(function, middle_point, 0, 2, 1000, val))
    Y3 = np.array(which_function(function, trapezoid, 0, 2, 1000, val))
    Y4 = np.array(which_function(function, Simpson, 0, 2, 1000, val))
    Y5 = np.array(which_function(function, Romberg, 0, 2, 15, val))
    Y6 = np.array(real_value(val, 1000))

    plt.plot(X, Y1, "_-", label = "Rectangle")
    plt.plot(X, Y2, "_-", label = "Point-Milieu")
    plt.plot(X, Y3, "_-", label = "Trapèze")
    plt.plot(X, Y4, "_-", label = "Simpson")
    plt.plot(X1, Y5, "_-", label = "Romberg")
    if(which_function == distance_iteration):
        plt.plot(X, Y6, "_-", label = "Valeur réelle")

    plt.xscale('log')
    if(which_function == absolute_error or which_function == relative_error):
        plt.yscale('log')

    plt.legend(prop = {'size' : 12})
    #plt.title("Longueur de la courbe de la fonction f en fonction de n")
    plt.grid()

    plt.show()

trace(distance_iteration, 4.64678)
#print(distance_iteration(function, Romberg_plot, 0, 2, 4, 4.64678))