import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt

# Wczytanie danych z plików tekstowych
data_alpha = np.genfromtxt("data/alpha_values.txt", skip_header=1)
data_eigenvalues = np.genfromtxt("data/eigenvalues.txt", skip_header=1)

# Wyodrębnienie wartości α i wartości pierwiastków własnych
alpha_values = data_alpha[:, np.newaxis]  # Dodanie nowej osi do tablicy
eigenvalues = data_eigenvalues[:, 1:]

# Tworzenie pierwszego wykresu: zmiana wartości własnych w funkcji parametru α
plt.figure(figsize=(8, 6))
for i in range(6):
    plt.plot(alpha_values, eigenvalues[:, i], label=f"eigenvalue_{i+1}")
plt.xlabel("Alpha")
plt.ylabel("Eigenvalues")
plt.title("Eigenvalues vs Alpha")
plt.legend()
plt.grid(True)
plt.show()

