import numpy as np
import matplotlib.pyplot as plt

# Wczytanie danych z plików
golden_data = np.loadtxt('data/golden_ratio_output.txt')
trisection_data = np.loadtxt('data/trisection_output.txt')

# Rozpakowanie danych
golden_iter, golden_xmin, golden_diff = golden_data[:, 0], golden_data[:, 1], golden_data[:, 2]
trisection_iter, trisection_xmin, trisection_diff = trisection_data[:, 0], trisection_data[:, 1], trisection_data[:, 2]

# Wykres
plt.figure(figsize=(10, 6))

# Wykres dla metody złotego podziału
plt.plot(golden_iter, golden_diff, label='Golden Ratio', marker='o')

# Wykres dla metody trisekcji
plt.plot(trisection_iter, trisection_diff, label='Trisection', marker='o')

# Skalowanie osi y na logarytmiczną
plt.yscale('log')

# Etykiety osi i tytuł
plt.xlabel('Numer Iteracji')
plt.ylabel('Moduł Różnicy')
plt.title('Moduł Różnicy w Funkcji Numeru Iteracji')

# Legenda
plt.legend()

# Wyświetlenie wykresu
plt.grid(True)
plt.show()
