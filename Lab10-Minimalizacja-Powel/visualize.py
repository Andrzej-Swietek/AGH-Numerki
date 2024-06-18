import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import os


parser = argparse.ArgumentParser(description='Results Visualization Lab 10',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--x_min', type=float, default=-1.5, help='Starting plot argument (x)')
parser.add_argument('--x_max', type=float, default=1.0, help='Ending plot argument (x)')
parser.add_argument("--save", default=True, help="Save plots to files")
parser.add_argument("--func", default='all', help="Which function is to be plotted")


args = parser.parse_args()

def plot_function_1():
    f1_v1_df = pd.read_csv('data/f1_results.csv', sep='\t')
    f1_v2_df = pd.read_csv('data/f1_results_2.csv', sep='\t')

    x_min_f1 = args.x_min
    x_max_f1 = args.x_max
    x = np.linspace(x_min_f1, x_max_f1)

    def f_x(x):
        return np.log(x**5 + 3*x**2 + x + 9)

    plt.rcParams.update({'font.size': 14})
    plt.plot(x, f_x(x), label='$f_1 (x) = \ln(x^5 + 3x^2 + x + 9)$')

    x = f1_v1_df["xm"]
    y = f_x(f1_v1_df["xm"])

    col_start = np.array([255, 0, 0]) / 255.0  # Convert to range [0, 1]
    col_end = np.array([0, 0, 255]) / 255.0    # Convert to range [0, 1]
    colors_v1 = [col_start + (col_end - col_start) * i / len(x) for i in range(len(x))]


    sc = plt.scatter(x, y, label='Kolejne przybliżenia minimum $x_m$ dla $x_1 = -0.5$', c=y, cmap='plasma', edgecolors='black', marker='o')
    plt.colorbar(sc, label='Kolejne przybliżenia minimum $x_m$ dla $x_1 = -0.5$')


    x = f1_v2_df["xm"]
    y = f_x(f1_v2_df["xm"])

    col_start = np.array([0, 125, 0]) / 255.0  # Convert to range [0, 1]
    col_end = np.array([0, 255, 0]) / 255.0    # Convert to range [0, 1]
    colors_v2 = [col_start + (col_end - col_start) * i / len(x) for i in range(len(x))]

    sc = plt.scatter(x, y, label='Kolejne przybliżenia minimum $x_m$ dla $x_1 = -0.9$', c=y, cmap='Greens', edgecolors='black' ,marker='o')
    plt.colorbar(sc, label='Kolejne przybliżenia minimum $x_m$ dla $x_1 = -0.9$')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('Szukanie minimum funkcji $f_1 (x) = \ln(x^5 + 3x^2 + x + 9)$')

    plt.show()



def plot_function_2():
    f2_df = pd.read_csv('data/f2_results.csv', sep='\t')

    x_min_f1 = -1.5
    x_max_f1 = 1
    x = np.linspace(x_min_f1, x_max_f1)

    def f_x(x):
        return x**6

    plt.rcParams.update({'font.size': 14})
    plt.plot(x, f_x(x), label='$f_2 (x) = x^6$')

    x = f2_df["xm"]
    y = f_x(f2_df["xm"])
    # colors = np.linspace(0, 1, len(x))
    col_start = np.array([0, 125, 0]) / 255.0  # Convert to range [0, 1]
    col_end = np.array([0, 255, 0]) / 255.0    # Convert to range [0, 1]
    colors = [col_start + (col_end - col_start) * i / len(x) for i in range(len(x))]

    sc = plt.scatter(x, y, label='Kolejne przybliżenia minimum xm', c=colors, marker='o')
    # plt.colorbar(sc, label='Różnica od minimum')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('Szukanie minimum funkcji $f_2 (x) = x^6 $')

    if args.save:
        plt.savefig('charts/function_2.png')

    plt.show()


def plot_iter_f1():
    f1_v1_df = pd.read_csv('data/f1_results.csv', sep='\t')

    x =[]
    for i in range(1,10+1):
        x.append(i)
    plt.plot(x, f1_v1_df["xm"], label="$x_m$ w zależności od iteracji")
    plt.scatter(x, f1_v1_df["xm"])

    plt.xlabel('x - Iteracja')
    plt.ylabel('y - Wartość $x_m$')
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('Szukanie minimum funkcji $f_1 (x)$')

    if args.save:
        plt.savefig('charts/function_1_iterations_v1.png')

    plt.show()

def plot_iter_f1_v2():
    f1_v1_df = pd.read_csv('data/f1_results_2.csv', sep='\t')

    x =[]
    for i in range(1,10+1):
        x.append(i)
    plt.plot(x, f1_v1_df["xm"], label="$x_m$ w zależności od iteracji")
    plt.scatter(x, f1_v1_df["xm"])

    plt.xlabel('x - Iteracja')
    plt.ylabel('y - Wartość $x_m$')
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('Szukanie minimum funkcji $f_1 (x) $')

    if args.save:
        plt.savefig('charts/function_1_iterations_v2.png')

    plt.show()

def plot_iter_f2():
    f1_v1_df = pd.read_csv('data/f2_results.csv', sep='\t')

    x =[]
    for i in range(1,100+1):
        x.append(i)
    plt.plot(x, f1_v1_df["xm"], label="$x_m$ w zależności od iteracji")
    plt.scatter(x, f1_v1_df["xm"])

    plt.xlabel('x - Iteracja')
    plt.ylabel('y - Wartość $x_m$')
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('Szukanie minimum funkcji $f_2 (x) $')

    if args.save:
        plt.savefig('charts/function_2_iterations.png')

    plt.show()


def difference_quotients_f1_v1():
    f1_v1_df = pd.read_csv('data/f1_results.csv', sep='\t')

    x =[]
    for i in range(1,10+1):
        x.append(i)
    plt.plot(x, f1_v1_df["F[x1,x2]"], label="Wartości ilorazów $F[x_1,x_2]$ w zależności od iteracji")
    plt.scatter(x, f1_v1_df["F[x1,x2]"])
    plt.plot(x, f1_v1_df["F[x1,x2,x3]"], label="Wartości ilorazów $F[x_1,x_2, x_3]$ w zależności od iteracji")
    plt.scatter(x, f1_v1_df["F[x1,x2,x3]"])

    plt.xlabel('x - Iteracja')
    plt.ylabel('y - Wartość Ilorazu Rożnicowego')
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('Wartości ilorazów $F[x_1,x_2]$ i $F[x_1,x_2,x_3]$ w zależności od iteracji dla funkcji $f_1(x)$')

    if args.save:
        plt.savefig('charts/function_1_quotients_v1.png')

    plt.show()

def difference_quotients_f1_v2():
    f1_v1_df = pd.read_csv('data/f1_results_2.csv', sep='\t')

    x =[]
    for i in range(1,10+1):
        x.append(i)
    plt.plot(x, f1_v1_df["F[x1,x2]"], label="Wartości ilorazów $F[x_1,x_2]$ w zależności od iteracji")
    plt.scatter(x, f1_v1_df["F[x1,x2]"])
    plt.plot(x, f1_v1_df["F[x1,x2,x3]"], label="Wartości ilorazów $F[x_1,x_2, x_3]$ w zależności od iteracji")
    plt.scatter(x, f1_v1_df["F[x1,x2,x3]"])

    plt.xlabel('x - Iteracja')
    plt.ylabel('y - Wartość Ilorazu Rożnicowego')
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('Wartości ilorazów $F[x_1,x_2]$ i $F[x_1,x_2,x_3]$ w zależności od iteracji dla funkcji $f_1(x)$')

    if args.save:
        plt.savefig('charts/function_1_quotients_v2.png')

    plt.show()

def difference_quotients_f2():
    f2_df = pd.read_csv('data/f2_results.csv', sep='\t')

    x =[]
    for i in range(1,100+1):
        x.append(i)

    plt.plot(x, f2_df["F[x1,x2]"], label="Wartości ilorazów $F[x_1,x_2]$ w zależności od iteracji")
    plt.scatter(x, f2_df["F[x1,x2]"])
    plt.plot(x, f2_df["F[x1,x2,x3]"], label="Wartości ilorazów $F[x_1,x_2, x_3]$ w zależności od iteracji")
    plt.scatter(x, f2_df["F[x1,x2,x3]"])

    plt.xlabel('x - Iteracja')
    plt.ylabel('y - Wartość Ilorazu Rożnicowego')
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('Wartości ilorazów $F[x_1,x_2]$ i $F[x_1,x_2,x_3]$ w zależności od iteracji dla $f_2(x)$')

    if args.save:
        plt.savefig('charts/function_2_quotients.png')

    plt.show()

if __name__ == '__main__':
    if args.func == 'all' or args.func == '1':
        plot_function_1()
    if args.func == 'all' or args.func == '2':
        plot_function_2()
    if args.func == 'all' or args.func == '3':
        plot_iter_f1()
    if args.func == 'all' or args.func == '4':
        plot_iter_f1_v2()
    if args.func == 'all' or args.func == '5':
        plot_iter_f2()
    if args.func == 'all' or args.func == '6':
        difference_quotients_f1_v1()
    if args.func == 'all' or args.func == '7':
        difference_quotients_f1_v2()
    if args.func == 'all' or args.func == '8':
        difference_quotients_f2()