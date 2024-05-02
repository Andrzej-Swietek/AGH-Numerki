import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

SAVE = True

# for i in [4, 2, 6]:
for i in [4]:
    df = pd.read_csv('data/approximation_results_{}_{}.csv'.format(i,i), sep=' ')

    x = df["x"]
    y = df["y"]

    plt.rcParams.update({'font.size': 14})
    plt.plot(x, np.cos(x), label='$f(x) = \\cos(x)$ - Oryginalna')
    plt.plot(x, y, label='$R_{' + str(i) + "," + str(i) + '}(x) $ - Aproksymowana')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Funkcja Aproksymowana dla N = {}, M = {}'.format(i,i))
    plt.legend()
    plt.grid(True)
    plt.gca().tick_params(width=1.5, length=5)

    # Dodanie siatki poprzez pogrubienie osi
    plt.axhline(0, color='black', linewidth=1.5)
    plt.axvline(0, color='black', linewidth=1.5)

    if SAVE:
        plt.savefig('charts/approximation_plot_{}_{}.png'.format(i,i))

    # Show plot
    plt.show()