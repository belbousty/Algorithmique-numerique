import numpy as np 

def function(x):
    return x**2 + 1

def f2(x):
    return 1/(1+x**2)

#METHOD IMPLEMENTATION

##rectangle method implementation
def rectangle(f, a, b, n):
    h = (b - a)/n
    somme = 0
    for i in range(n):
        somme += f(a + i*h)
    return h*somme; 


##middle point method implementation
def middle_point(f, a, b, n):
    h = (b - a)/n
    somme = 0
    for i in range(n):
        somme += f(a + i*h + h/2)
    return h*somme


##trapezoidal method implementation
def trapezoid(f, a, b, n):
    h = (b - a)/n
    somme = 0
    for i in range(n):
        somme += f(a + i*h)
    return h*((f(a)+f(b))/2 + somme)


##simpson method implementation
def Simpson(f, a, b, n):
    h = (b - a)/n
    somme1 = 0
    somme2 = 0
    for i in range(n):
        somme1 += f(a + i*h)
        somme2 += f(a + i*h + h/2)
    return h*((f(a)+f(b))/6 + somme1*1/3 + somme2*2/3)



def sommeRomberg(f, a, N, h):
    somme = 0
    for i in range(1, 2**(N-1) + 1):
        somme += f(a + (2*i - 1)*h)
    return somme


##romberg method implementation
def Romberg(f, a, b, n):
    h = (b-a)/2
    R = np.zeros((n, n))
    R[0, 0] = (f(a)+f(b))*(b-a)/2
    for i in range(1, n):
        R[i, 0] = R[i-1, 0]*0.5 + h*sommeRomberg(f, a, i, h)
        for j in range(1, i + 1):
            R[i, j] = ((4**j)*R[i, j-1] - R[i - 1, j - 1])/(4**j - 1)
        h /= 2
    return R[n-1, n-1]


##return the function f'
def derived_f(f):
    def df(x):
        return (f(x + 1e-3) - f(x))/1e-3
    return df


## return the function x -> sqrt(1+f'(x)**2)
def sqrt_length(f):
    df = derived_f(f)
    def sq(x):
        return np.sqrt(1 + df(x)**2)
    return sq


## return the length curve of the function f between a and b 
## according to the method "method" and with n points
def distance(f, method, a, b, n):    
    return method(sqrt_length(f), a, b, n)

if __name__ == "__main__":
    print("Rectangle Method : ", distance(function, rectangle, 0, 2, 1000))
    print("Middle point method : ", distance(function, middle_point, 0, 2, 1000))
    print("Trapeziodal method :", distance(function, trapezoid, 0, 2, 1000))
    print("Simpson method : ", distance(function, Simpson, 0, 2, 1000))
    print("Méthode de Romberg avec function : ", distance(function, Romberg, 0, 2, 15))
    #print("Méthode de Romberg vérif cours : ", Romberg(f2, -1, 2, 4))


