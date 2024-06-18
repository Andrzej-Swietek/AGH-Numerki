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

def gen_plot(gen_no: int, r: int) -> None:
    colors = ['', 'blue', 'orange', 'purple']
    gen1_df = pd.read_csv(f"data/wykres_gen_{gen_no}r_{r}.csv", sep=" ")
    plt.rcParams.update({'font.size': 14})

    # Label the axes
    plt.xlabel('$x_i$')
    plt.ylabel('$x_{i+' f'{r}' '}$')
    plt.scatter(gen1_df["x_1"], gen1_df["x_2"], label=f"generator {gen_no} r={r}", c=colors[gen_no])

    # Adding grid, legend, and title
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title("Wykres zależności $x_{i+" f'{r}' "}$ od $x_i$ dla generatora" + f'{gen_no}')

    # Save the plot if specified
    if args.save:
        if not os.path.exists('charts'):
            os.makedirs('charts')
        plt.savefig(f'charts/gen_{gen_no}_{r}.png')

    # Show the plot
    plt.show()


if __name__ == '__main__':
    if args.func == '1':
        gen_plot(1,1)
        gen_plot(1,2)
        gen_plot(1,3)
    if args.func == '2':
        gen_plot(2,1)
        gen_plot(2,2)
        gen_plot(2,3)
    if args.func == '3':
        gen_plot(3,1)
        gen_plot(3,2)
        gen_plot(3,3)
    pass