import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

SAVE = True

for i in [5, 10, 15]:
    nodes = pd.read_csv('data/nodes_xi_{}_czebyszew.txt'.format(i), sep=' ')

    df = pd.read_csv('data/interpolation_results_{}_czebyszew.txt'.format(i), sep=' ')

    # Extracting columns
    x = df["x"]
    y = df["y"]
    interpolated = df["interpolated"]

    plt.rcParams.update({'font.size': 14})

    # Plotting
    plt.plot(x, y, label='$f(x) = \\frac{x}{1+x^2} $ - Oryginalna')
    plt.plot(x, interpolated, label='F(x) - Funkcja Interpolowana')

    plt.scatter(nodes["x"], nodes["y"], label='Punkty Węzłowe $x_i$', color='red', marker='o')

    # Adding labels and legend
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Funkcja Interpolowana (Czebyszew) dla N = {}'.format(i))
    plt.legend()

    plt.grid(True)

    # Save plot to file
    if SAVE:
        plt.savefig('charts/interpolation_plot_n_{}_czebyszew.png'.format(i))

    # Show plot
    plt.show()
