from cmath import exp
from database_to_array import *
from scipy import interpolate
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})



#return cubic interpolation of (X,Y) coordinates with scipy method
def f(x_points, y_points, method = "Not scipy"):
    if method == "Scipy":
        def F(x):
            tck = interpolate.splrep(x_points, y_points)
            return interpolate.splev(x, tck)
        return F
    else:
        return splint(x_points, y_points, spline(x_points, y_points, 0, 0))




## return the second derivate f''(x) from (xi,yi) and f''(x1) and f''(x2)
def spline(x, y, yp1, ypn):
    n = len(x)
    u = np.zeros((n, 1))
    y2 = np.zeros((n, 1))
    if yp1 > 0.99e30:
        y2[0] = 0
        u[0] = 0
    else:
        y2[0] = -0.5
        u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0])/x[1] - x[0] - yp1)

    for i in range(1, n-1):
        sig = (x[i]-x[i-1])/(x[i+1]-x[i-1])
        p = sig*y2[i-1]+2.0
        y2[i]=(sig-1.0)/p
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1])
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p
    
    if ypn > 0.99e30:
        qn = 0
        un = 0
    else:
        qn = 0.5
        un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
    
    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0)
    k = n-2
    while k >= 0:
        y2[k]=y2[k]*y2[k+1]+u[k]
        k -= 1
    
    return y2


##return the cubic spline interpolation from (x,y) and y''
def splint(xa, ya, y2a):

    def f(x):
        klo = 0
        khi = len(xa) - 1
        while (khi - klo) > 1:
            k = (khi + klo) >> 1
            if xa[k] > x:
                khi = k
            else:
                klo = k
        h = xa[khi] - xa[klo]
        if h == 0:
            return None
        a = (xa[khi] - x)/h
        b = (x - xa[klo])/h
        return a*ya[klo] + b*ya[khi] + ((a**3 - a)*y2a[klo] + (b**3 - b)*y2a[khi])*h**2/6
    
    return f




##Compute the cubic spline interpolatio
##Returns (f_e, f_i) two functions of the curve of respectiveley the extrados and the intrados
def get_airfoil_curves(file, method = "not_scipy"):
    (dim,ex,ey,ix,iy) = load_foil(file)
    if (method == "scipy"):
        return (f(ex, ey), f(ix, iy))
    else:
        return (splint(ex, ey, spline(ex, ey, 0, 0)), splint(ix, iy, spline(ix, iy, 0, 0)))
        
    


#get min and max of the airfoil (h_min and h_max)
def get_Max_Min(file):
    (dim,ex,ey,ix,iy) = load_foil(file)
    return (max(ey), min(iy))


def get_Dimension(file):
    (dim,ex,ey,ix,iy) = load_foil(file)
    return dim


#get min and max of the airfoil (h_min and h_max)
def get_Max_Min(file):
    (dim,ex,ey,ix,iy) = load_foil(file)
    return (max(ey), min(iy))


def plot_airfoil(file):
    (dim,ex,ey,ix,iy) = load_foil(file)
    
    ex_pts = np.linspace(ex[0],ex[-1],10000)
    ix_pts = np.linspace(ix[0],ix[-1],10000)

    ecs = [f(ex,ey)(x) for x in ex_pts] 
    ics =  [f(ix,iy)(x) for x in ix_pts]


    fig, ax = plt.subplots()
    ax.plot(ex, ey, 'o', label='data extrados')
    ax.plot(ix, iy, 'o', label ='data intrados')
    ax.plot(ex_pts, ecs, label='extrados')
    ax.plot(ix_pts, ics, label = 'intrados')
    plt.title("Profil aérien d'une aile d'un boeing 737")
    plt.legend()
    plt.show()


def compare_tab(T1, T2, e):
    for i in range(len(T1)):
        if abs(T1[i] - T2[i]) > e:
            return False
    return True


def test_csi(file):
    print("TEST CUBIC SPLINE INTERPOLATION")

    (dim,ex,ey,ix,iy) = load_foil(file)
    ex_pts = np.linspace(ex[0],ex[-1],100)
    ix_pts = np.linspace(ix[0],ix[-1],100)

    F_scipy_e = f(ex,ey)
    F_scipy_i = f(ix,iy)

    F_test_e = splint(ex, ey, spline(ex, ey, 0, 0))
    F_test_i = splint(ix, iy, spline(ix, iy, 0, 0))

    e_scipy = [F_scipy_e(x) for x in ex_pts] 
    i_scipy =  [f(ix,iy)(x) for x in ix_pts]

    e_test = [F_test_e(x) for x in ex_pts]
    i_test = [F_test_i(x) for x in ix_pts]

    assert(compare_tab(e_scipy, e_test, 0.01))
    assert(compare_tab(i_scipy, i_test, 0.01))

    print("SUCCESS")

def plot_test_csi(file):
    (dim,ex,ey,ix,iy) = load_foil(file)
    ex_pts = np.linspace(ex[0],ex[-1],100)
    ix_pts = np.linspace(ix[0],ix[-1],100)

    F_scipy_e = f(ex,ey)
    F_scipy_i = f(ix,iy)

    F_test_e = splint(ex, ey, spline(ex, ey, 0, 0))
    F_test_i = splint(ix, iy, spline(ix, iy, 0, 0))

    e_scipy = [F_scipy_e(x) for x in ex_pts] 
    i_scipy =  [f(ix,iy)(x) for x in ix_pts]

    e_test = [F_test_e(x) for x in ex_pts]
    i_test = [F_test_i(x) for x in ix_pts]

    fig, ax = plt.subplots()
    ax.plot(ex, ey, 'o', label='data extrados')
    ax.plot(ix, iy, 'o', label ='data intrados')
    ax.plot(ex_pts, e_scipy, label='true_extrados')
    ax.plot(ix_pts, i_scipy, label = 'true_intrados')

    ax.plot(ex_pts, e_test, label='test_extrados')
    ax.plot(ix_pts, i_test, label = 'test_intrados')
    plt.show()



if __name__ == "__main__":
    if len(sys.argv) == 1:
        file = "b737a.dat"
    else:
        file = sys.argv[1]

    plot_airfoil(file)
    plot_test_csi(file)
    test_csi(file)