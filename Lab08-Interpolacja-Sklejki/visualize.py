import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

SAVE = True

for i in [5, 6, 7, 10, 20]:
    df = pd.read_csv('data/interpolation_results_{}.txt'.format(i), sep=' ')
    nodes = pd.read_csv('data/nodes_xi_{}.txt'.format(i), sep=' ')

    x = df["x"]
    y = df["y"]

    plt.rcParams.update({'font.size': 14})
    plt.plot(x, 1/(1+x**2), label='$f(x) = \\frac{1}{1+x^2} $ - Oryginalna')
    plt.plot(x, y, label='$f(x) = \\frac{1}{1+x^2} $ - Interpolowana')
    plt.scatter(nodes["x"], nodes["y"], label='Punkty Węzłowe $x_i$', color='red', marker='o')
    plt.scatter(nodes["x"], 1/(1+nodes["x"]**2), label='Punkty Węzłowe $x_i$', color='purple', marker='o')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Funkcja Interpolowana dla N = {}'.format(i))
    plt.legend()
    plt.grid(True)

    if SAVE:
        plt.savefig('charts/interpolation_plot_n_{}.png'.format(i))

    # Show plot
    plt.show()


for i in [ 6, 7, 14 ]:
    df = pd.read_csv('data/interpolation_results_{}f_2.txt'.format(i), sep=' ')
    nodes = pd.read_csv('data/nodes_xi_{}f_2.txt'.format(i), sep=' ')

    x = df["x"]
    y = df["y"]

    plt.rcParams.update({'font.size': 14})
    plt.plot(x, np.cos(2*x), label='$f(x) = \\cos(2x) $ - Oryginalna')
    plt.plot(x, y, label='$f(x) = \\cos(2x) $ - Interpolowana')
    plt.scatter(nodes["x"], nodes["y"], label='Punkty Węzłowe $x_i$', color='red', marker='o')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Funkcja Interpolowana dla N = {}'.format(i))
    plt.legend()
    plt.grid(True)

    if SAVE:
        plt.savefig('charts/interpolation_plot_n_{}_f_2.png'.format(i))

    # Show plot
    plt.show()