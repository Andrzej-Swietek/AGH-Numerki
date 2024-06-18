import os
import pandas as pd
import matplotlib.pyplot as plt

folder_path = './data' 

file_list = [file for file in os.listdir(folder_path) if file.endswith('.csv')]


for file in file_list:
    file_path = os.path.join(folder_path, file)
 
    df = pd.read_csv(file_path, sep='\t')

    is_modified = '_m' in file

    title_suffix = "pierwszego" if "_1" in file else "drugiego" if "_2" in file else "trzeciego"


    title_suffix += " (f. modyfikowana)" if is_modified else ""


    fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    # Subplot dla zmiany x wraz z iteracją
    axs[0].plot(df['iteracja'], df['x'], marker='o', color='b', label='Wartość x')
    axs[0].set_title(f'Zmiana wartości x wraz z iteracją dla {title_suffix}')
    axs[0].set_xlabel('Iteracja')
    axs[0].set_ylabel('Wartość x')
    axs[0].grid(True)
    axs[0].legend()

    # Subplot dla zmiany epsilon wraz z iteracją
    axs[1].plot(df['iteracja'], df['epsilon'], marker='s', color='r', label='Wartość epsilon')
    axs[1].set_title(f'Zmiana wartości epsilon wraz z iteracją dla {title_suffix}')
    axs[1].set_xlabel('Iteracja')
    axs[1].set_ylabel('Wartość epsilon')
    axs[1].grid(True)
    axs[1].legend()

    # Ustawienie większych czcionek dla osi i legendy
    for ax in axs:
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.tick_params(axis='both', which='minor', labelsize=10)
        ax.legend(fontsize=12)

    # Wyświetlenie wykresu
    plt.tight_layout()
    plt.savefig(file_path.replace('.csv', '.png'))  # Zapisanie wykresu do pliku PNG
    plt.close()  # Zamknięcie aktualnego wykresu

print("Wygenerowano wykresy dla wszystkich plików CSV.")

