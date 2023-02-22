import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import time

def step_euler(y, t, h, f):
    """
    Computes a single step by using Euler's method
    :param y: initial
    :param t: initial
    :param h: time step
    :param f: function
    :return: y, t
    """
    y += h * f(y, t)
    t += h

    return y, t

def step_middle_point(y, t, h, f):
    """
    Computes a single step by using the middle point method
    :param y: initial
    :param t: initial
    :param h: time step
    :param f: function
    :return: y, t
    """
    y_middle = y + h * f(y, t)

    y += h * f(y_middle, t + h)
    t += h

    return y, t

def step_heun(y, t, h, f):
    """
    Computes a single step by using the Heun method
    :param y: initial
    :param t: initial
    :param h: time step
    :param f: function
    :return: y, t
    """
    y_intermediate = y + h * f(y, t)

    y = y_intermediate + (h/2) * (f(y, t) + f(y_intermediate, t + h))
    t += h

    return y, t

def step_rk4(y, t, h, f):
    """
    Computes a single step by using the Runge-Kutta method
    :param y: initial
    :param t: initial
    :param h: time step
    :param f: function
    :return: y, t
    """
    k1 = f(y, t)
    k2 = f(y + h/2 * k1, t + h/2)
    k3 = f(y + h/2 * k2, t + h/2)
    k4 = f(y + h * k3, t + h)

    y_plus = y + (h/6) * (k1 + 2 * k2 + 2 * k3 + k4)
    t += h

    return y_plus, t

def meth_n_step(y0, t0, N, h, f, meth):
    """
    Computes N steps by using the chosen method by using steps of size h
    :param y0: initial
    :param t0: initial
    :param N: number of steps
    :param h: time step
    :param f: function
    :param meth: method
    :return: y, t
    """
    y = y0
    t = t0

    for _ in range(N):
        y, t = meth(y, t, h, f)

    return y, t

def meth_epsilon(y0, t0, tf, eps, f, meth):
    """
    Computes N steps by using the chosen method by using steps of size h
    Uses an error parameter epsilon to stop the loop
    :param y0: initial
    :param t0: initial
    :param tf: final
    :param eps: error parameter
    :param f: function
    :param meth: method
    :return: y, t
    """
    y = y0
    t = t0

    # Start by defining the step size h
    h = (tf - t0) / 100

    while abs(y - tf) > eps:
        y, t = meth(y, t, h, f)

    return y, t

def point_milieu(f, y0, t):
    """
    Computes the solution of the ODE using the median point method
    """
    n = len(t)
    sol = [y0]
    for k in range(n-1):
        h = t[k+1] - t[k]
        y_mid = sol[k] + (h/2) * f(sol[k], t[k])
        yp_1 = f(y_mid, t[k] + (h/2))
        y_1 = sol[k] + h * yp_1
        sol.append(y_1)

    return sol

def euler(f, y0, t):
    """
    Computes the solution of the ODE using the Euler method
    """
    n = len(t)
    sol = [y0]
    for k in range(n-1):
        h = t[k+1] - t[k]
        y_1 = sol[k] + h * f(sol[k], t[k])
        sol.append(y_1)

    return sol

def heun(f, x0, t):
    """
    Computes the solution of the ODE using the Heun method
    """
    n = len(t)
    sol = [x0]
    for k in range(0, n-1):
        h = t[k+1] - t[k]
        p1 = f(sol[k], t[k])
        p2 = f(sol[k] + h * p1, t[k+1])
        sol.append(sol[k] + h * (p1 + p2) / 2)
    return sol

def rk4(f, y0, t):
    """
    Computes the solution of the ODE using the Runge-Kutta method
    """
    n = len(t)
    sol = [y0]
    for k in range(0, n-1):
        h = t[k+1] - t[k]
        pn1 = f(sol[k], t[k])
        yn2 = sol[k] + (1/2) * h * pn1
        pn2 = f(yn2, t[k] + (1/2) * h)
        yn3 = sol[k] + (1/2) * h * pn2
        pn3 = f(yn3, t[k] + (1/2) * h)
        yn4 = sol[k] + h * pn3
        pn4 = f(yn4, t[k] + h)
        sol.append(sol[k] + (h/6) * (pn1 + 2 * pn2 + 2 * pn3 + pn4))

    return sol

def main():
    # Implement the first differential equation
    # Then plot the solution
    # case y(0) = 1
    # case y'(t) = y(t) / (1 + t^2)
    f = lambda y, t: y / (1 + t**2)
    y = 1

    _, axes = plt.subplots(1, 2, sharey="row")

    # Define two subplots with different scales
    axis1 = axes[0]
    axis2 = axes[1].twinx()

    # Prepare to store the time results with the method's
    # name as the key
    times = {}

    t = np.linspace(0, 50, 200)
    axis1.plot(t, odeint(f, y, t), '--', label='Solution exacte')

    for k, v in {"euler": euler, "heun": heun, "rk4": rk4, "points_milieu": point_milieu}.items():
        # Measure the time needed to compute the solution
        start = time.time()

        # Compute the solution
        axis1.plot(t, v(f, y, t), label="Méthode de {}".format(k))

        # Store the time result
        times[k] = time.time() - start

    # Set the title and position the legend in the upper left corner
    axis1.title.set_text("Comparaison des méthodes")
    axis1.legend(loc='upper left')

    # Create a bar graph with the time results
    axis2.bar(list(times.keys()), list(times.values()), align='center')
    axis2.title.set_text("Temps de calcul (s)")

    # Show the graph
    plt.show()
    plt.close()

    # Implement the second differential equation of dimension 2
    # Then plot the solution
    # y(t) = [y_1(t) y_2(t)]
    # case y(0) = [1 0]
    # case y'(t) = [-y_2(t) y_1(t)]
    y = np.array([0, 0])
    t = 0
    f = lambda y, _: np.array([-y[1], y[0]])

    t = np.linspace(0, 50, 200)

    # Plot the exact solution
    plt.plot(t, odeint(f, y, t), '--', label='Solution exacte')

    # Plot the solutions of the different methods
    for k, v in {"euler": euler, "heun": heun, "rk4": rk4, "point milieu": point_milieu}.items():
        plt.plot(t, v(f, y, t), label="Méthode de {}".format(k))
    
    # Set the title and position the legend in the upper left corner
    plt.title("Comparaison des méthodes".format(k))
    plt.legend(loc='upper left')

    # Show the graph
    plt.show()


if __name__ == "__main__":
    main()
