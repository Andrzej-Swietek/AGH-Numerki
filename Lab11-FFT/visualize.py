import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import os
import argparse
import os


parser = argparse.ArgumentParser(description='Results Visualization Lab 10',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--save", default=True, help="Save plots to files")
parser.add_argument("--func", default='all', help="Which function is to be plotted")
parser.add_argument("--k", default='all', help="Which probe pram is to be plotted")
args = parser.parse_args()


PI = np.pi
TWO_PI = 2.0 * np.pi

def plot_signals(k: int) -> None:
    df_clean = pd.read_csv(f"data/fft_clean_signal_{k}.csv", sep=" ")
    df_noiced = pd.read_csv(f"data/fft_noiced_signal_{k}.csv", sep=" ")

    plt.rcParams.update({'font.size': 14})
    plt.plot(df_clean["i"], df_clean["signal"], label='Clean Signal', c="red")
    plt.scatter(df_noiced["i"], df_noiced["signal"], label='Noiced Signal')
    plt.xlabel('i-ty element dyskretnego sygnalu')
    plt.ylabel('Wartość sygnału')
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('FFT - Sygnały czysty i zaszumiony')

    if args.save:
        plt.savefig(f'charts/fft-signals-{k}.png')

    plt.show()

def plot_transformed_Re_Im(k: int) -> None:
    df_result = pd.read_csv(f"data/fft_result_transformed_{k}.csv", sep=" ")

    plt.rcParams.update({'font.size': 14})
    plt.plot(df_result["i"], df_result["Im"], label='Im Signal', c="red")
    plt.plot(df_result["i"], df_result["Re"], label='Re Signal')

    plt.xlabel('i-ty element dyskretnego sygnalu')
    plt.ylabel('Wartość po transformacie')
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('FFT - Część rzeczywista i urojona sygnału po transformacie')

    if args.save:
        plt.savefig(f'charts/fft-transformed-re-im-{k}.png')

    plt.show()

def plot_threshold(k: int) -> None:
    df_result = pd.read_csv(f"data/fft_result_transformed_{k}.csv", sep=" ")

    plt.rcParams.update({'font.size': 14})
    # plt.plot(df_result["i"], df_result["Im"], label='Im Signal', c="red")
    # plt.plot(df_result["i"], df_result["Re"], label='Re Signal')

    df_result['Magnitude'] = np.sqrt(df_result['Re']**2 + df_result['Im']**2)
    plt.plot(df_result["i"], df_result["Magnitude"], label="|FFT|")
    threshold = df_result['Magnitude'].max() / 2

    # Rysowanie poziomej linii threshold
    plt.axhline(y=threshold, color='green', linestyle='--', label=f'Threshold = {threshold:.2f}')
    plt.xlabel('i')
    plt.ylabel('y')
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('FFT - Próg dykryminacji')

    if args.save:
        plt.savefig(f'charts/fft-transformed-treshold-{k}.png')

    plt.show()


# Definiowanie funkcji sygnału
def signal_fn(i, omega):
    return np.sin(omega * i) + np.sin(2 * omega * i) + np.sin(3 * omega * i)

# Generowanie oryginalnego sygnału
def generate_periodic_signal(N: int):
    omega = TWO_PI / N
    signal = [signal_fn(i, omega) for i in range(N)]
    return signal

def plot_cleaned(k: int) -> None:
    df_result = pd.read_csv(f"data/fft_result_wynik{k}.csv", sep=" ")
    N = 1 << k  # N = 2^k
    original_signal = generate_periodic_signal(N)

    df_result["Oryginal"] = original_signal
    plt.rcParams.update({'font.size': 14})
    plt.plot(df_result["i"], df_result["Re"], label='Sygnał oczyszczony')
    plt.plot(df_result["i"], df_result["Oryginal"], label="Sygnał oryginalny")

    plt.xlabel('i-ty element dyskretnego sygnalu')
    plt.ylabel('Wartość sygnału')
    plt.grid(True)
    plt.legend()
    plt.gca().tick_params(width=1.5, length=5)
    plt.title('FFT - Sygnały Oczyszczony i oryginalny')

    if args.save:
        plt.savefig(f'charts/fft-cleaned-{k}.png')
    plt.show()


if __name__ == '__main__':
    probing = []
    if args.k == 'all' or args.k == '8': probing.append(8)
    if args.k == 'all' or args.k == '10': probing.append(10)
    if args.k == 'all' or args.k == '12': probing.append(12)

    for k in probing:
        if args.func == 'all' or args.func == '1': plot_signals(k)
        if args.func == 'all' or args.func == '2': plot_transformed_Re_Im(k)
        if args.func == 'all' or args.func == '3': plot_threshold(k)
        if args.func == 'all' or args.func == '4': plot_cleaned(k)
