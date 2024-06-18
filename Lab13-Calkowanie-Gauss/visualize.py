import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import os

# Parser to handle command line arguments
parser = argparse.ArgumentParser(description='Results Visualization Lab 13',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--save", default=True, help="Save plots to files")
parser.add_argument("--func", default='all', help="Which function is to be plotted")
args = parser.parse_args()

PI = np.pi
TWO_PI = 2.0 * np.pi

def plot_integral_1() -> None:
    # Read data from the CSV file
    df_result = pd.read_csv("data/integral_1.csv", sep="\t")

    plt.rcParams.update({'font.size': 14})

    # Plot Approx Integral vs. n
    plt.plot(df_result["n"], df_result["Approx Integral"], label="Approx Integral")
    plt.scatter(df_result["n"], df_result["Approx Integral"])

    # Label the axes
    plt.xlabel('Kolejna iteracja')
    plt.ylabel('Wartość przybliżenia całki')

    # Adding grid, legend, and title
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('Całkowanie przy użyciu kwadratury Gaussa-Legendre’a')

    # Save the plot if specified
    if args.save:
        if not os.path.exists('charts'):
            os.makedirs('charts')
        plt.savefig('charts/integral_1.png')

    # Show the plot
    plt.show()


def plot_integral_1_error() -> None:
    # Read data from the CSV file
    df_result = pd.read_csv("data/integral_1.csv", sep="\t")

    plt.rcParams.update({'font.size': 14})

    # Plot Approx Integral vs. n
    plt.plot(df_result["n"], df_result["Error"], label="$|c_1 - c_{1,a}|$")
    plt.scatter(df_result["n"], df_result["Error"])

    # Label the axes
    plt.xlabel('Kolejna iteracja')
    plt.ylabel('Wartość przybliżenia całki')

    # Adding grid, legend, and title
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('Błąd Całkowania przy użyciu kwadratury Gaussa-Legendre’a')

    # Save the plot if specified
    if args.save:
        if not os.path.exists('charts'):
            os.makedirs('charts')
        plt.savefig('charts/integral_1_error.png')

    # Show the plot
    plt.show()

def plot_integral_2a() -> None:
    # Read data from the CSV file
    df_result = pd.read_csv("data/integral_2_a.csv", sep="\t")

    plt.rcParams.update({'font.size': 14})

    # Plot Approx Integral vs. n
    plt.plot(df_result["n"], df_result["Exact Integral"], label="Exact Integral", c='blue')
    plt.plot(df_result["n"], df_result["Approx Integral"], label="Approx Integral", c='orange')
    plt.scatter(df_result["n"], df_result["Approx Integral"], c='orange')

    # Label the axes
    plt.xlabel('Kolejna iteracja')
    plt.ylabel('Wartość przybliżenia całki')

    # Adding grid, legend, and title
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('Całkowanie przy użyciu kwadratury Gaussa-Laguerre’a')

    # Save the plot if specified
    if args.save:
        if not os.path.exists('charts'):
            os.makedirs('charts')
        plt.savefig('charts/integral_2a.png')

    # Show the plot
    plt.show()


def plot_integral_2a_error() -> None:
    # Read data from the CSV file
    df_result = pd.read_csv("data/integral_2_a.csv", sep="\t")

    plt.rcParams.update({'font.size': 14})

    # Plot Approx Integral vs. n
    plt.plot(df_result["n"], df_result["Error"], label="$|c_2 - c_{2,a}|$")
    plt.scatter(df_result["n"], df_result["Error"])

    # Label the axes
    plt.xlabel('Kolejna iteracja')
    plt.ylabel('Błąd przybliżenia całki')

    # Adding grid, legend, and title
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('Błąd Całkowania przy użyciu kwadratury Gaussa-Legendre’a')

    # Save the plot if specified
    if args.save:
        if not os.path.exists('charts'):
            os.makedirs('charts')
        plt.savefig('charts/integral_2a_error.png')

    # Show the plot
    plt.show()

def plot_integral_2b_error() -> None:
    # Read data from the CSV file
    df_result = pd.read_csv("data/integral_2_b.csv", sep="\t")

    plt.rcParams.update({'font.size': 14})

    # Plot Approx Integral vs. n
    plt.plot(df_result["n"], df_result["Error"], label="$|c_2 - c_{2,a}|$")
    plt.scatter(df_result["n"], df_result["Error"])

    # Label the axes
    plt.xlabel('Kolejna iteracja')
    plt.ylabel('Błąd przybliżenia całki')

    # Adding grid, legend, and title
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('Błąd Całkowania przy użyciu kwadratury Gaussa-Legendre’a')

    # Save the plot if specified
    if args.save:
        if not os.path.exists('charts'):
            os.makedirs('charts')
        plt.savefig('charts/integral_2b_error.png')

    # Show the plot
    plt.show()

def plot_integral_2b() -> None:
    # Read data from the CSV file
    df_result = pd.read_csv("data/integral_2_b.csv", sep="\t")

    plt.rcParams.update({'font.size': 14})

    # Plot Approx Integral vs. n
    plt.plot(df_result["n"], df_result["Exact Integral"], label="Exact Integral", c='blue')
    plt.plot(df_result["n"], df_result["Approx Integral"], label="Approx Integral", c='orange')
    plt.scatter(df_result["n"], df_result["Approx Integral"], c='orange')

    # Label the axes
    plt.xlabel('Kolejna iteracja')
    plt.ylabel('Wartość przybliżenia całki')

    # Adding grid, legend, and title
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('Całkowanie przy użyciu kwadratury Gaussa-Laguerre’a')

    # Save the plot if specified
    if args.save:
        if not os.path.exists('charts'):
            os.makedirs('charts')
        plt.savefig('charts/integral_2b.png')

    # Show the plot
    plt.show()


def plot_integral_3() -> None:
    # Read data from the CSV file
    df_result = pd.read_csv("data/integral_3.csv", sep="\t")

    print(df_result.head())
    plt.rcParams.update({'font.size': 14})

    # Plot Approx Integral vs. n
    plt.plot(df_result["n"], df_result["Exact Integral"], label="Exact Integral", c='blue')
    plt.plot(df_result["n"], df_result["Approx Integral"], label="Approx Integral", c='orange')
    plt.scatter(df_result["n"], df_result["Approx Integral"], c='orange')

    # Label the axes
    plt.xlabel('Kolejna iteracja')
    plt.ylabel('Wartość przybliżenia całki')

    # Adding grid, legend, and title
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('Całkowanie przy użyciu kwadratury Gaussa-Hermite’a')

    # Save the plot if specified
    if args.save:
        if not os.path.exists('charts'):
            os.makedirs('charts')
        plt.savefig('charts/integral_3.png')

    # Show the plot
    plt.show()

def plot_integral_3_error() -> None:
    # Read data from the CSV file
    df_result = pd.read_csv("data/integral_3.csv", sep="\t")

    plt.rcParams.update({'font.size': 14})

    # Plot Approx Integral vs. n
    plt.plot(df_result["n"], df_result["Error"], label="$|c_3 - c_{3,a}|$")
    plt.scatter(df_result["n"], df_result["Error"])

    # Label the axes
    plt.xlabel('Kolejna iteracja')
    plt.ylabel('Błąd przybliżenia całki')

    # Adding grid, legend, and title
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('Błąd Całkowania przy użyciu kwadratury Gaussa-Hermite’a')

    # Save the plot if specified
    if args.save:
        if not os.path.exists('charts'):
            os.makedirs('charts')
        plt.savefig('charts/integral_3_error.png')

    # Show the plot
    plt.show()

if __name__ == '__main__':
    if args.func == '1':
        plot_integral_1()
    if args.func == '1-error':
        plot_integral_1_error()
    if args.func == '2a':
        plot_integral_2a()
    if args.func == '2b':
        plot_integral_2b()
    if args.func == '2-error':
        plot_integral_2a_error()
        plot_integral_2b_error()
    if args.func == '3':
        plot_integral_3()
    if args.func == '3-error':
        plot_integral_3_error()
    pass
