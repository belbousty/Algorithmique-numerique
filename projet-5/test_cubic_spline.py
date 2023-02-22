from asyncio import new_event_loop
from airfoil_refinement import *
import numpy as np
import matplotlib.pyplot as plt
from database_to_array import *
plt.rcParams.update({'font.size': 15})


def plot_csi_nb_points():
    sin = np.sin
    X_N = [5, 10, 20]
    fig, ax = plt.subplots()
    X_plot = np.linspace(0, 2*np.pi, 1000)

    ax.plot(X_plot, sin(X_plot), label="sin(x)")

    for i in range(len(X_N)):
        X = np.linspace(0, 2*np.pi, X_N[i])
        F_test_e = f(X, [sin(x) for x in X])
        ax.plot(X_plot, [F_test_e(x) for x in X_plot], label="N = " + str(X_N[i]))

    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.title("Fonction sinus selon le nombre de points N en entrée de l'interpolation en spline cubique")
    plt.show()


def take_n_points_unformely(tab1, tab2, n):
    if n == 0:
        return ([], [])
    elif n == len(tab1):
        return (tab2, tab2)
    else:
        cpt = 0
        new_tab1 = []
        new_tab2 = []
        idx = 0
        while cpt < n or idx < len(tab1):
            new_tab1.append(tab1[cpt])
            new_tab2.append(tab2[cpt])
            cpt += 1
            idx += len(tab1)//n
        return (new_tab1, new_tab2)


def plot_error_cubic_spline():
    (dim,ex,ey,ix,iy) = load_foil("./b737a.dat")
    X_test = np.linspace(ex[0], ex[-1], 100)
    X = [k for k in range(2, len(ex))]
    E_real = f(ex, ey, "scipy")
    error = []
    for k in X:
        (x,y) = take_n_points_unformely(ex, ey, k)
        E = f(x, y)
        error.append(np.mean([abs(E_real(x_test) - E(x_test)) for x_test in X_test]))

    plt.xlabel("Pourcentage des valeurs prises")
    plt.ylabel("Erreur directe moyenne")
    plt.title("Erreur directe moyenne en fonction du pourcentage du nombre de points N")
    plt.plot([x/len(ex) * 100 for x in X], error)
    # plt.xscale("log")
    plt.yscale("log")
    plt.show()








if __name__ == "__main__":
    plot_csi_nb_points()
    plot_error_cubic_spline()