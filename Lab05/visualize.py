import numpy as np
import matplotlib.pyplot as plt

# Funkcja do odczytu wektorów własnych z pliku
def read_eigenvectors(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        eigenvectors = []
        for line in lines:
            if line.startswith("Eigenvalue"):
                eigenvectors.append([])
            else:
                eigenvectors[-1].append(float(line.strip()))
    return eigenvectors


# Odczytanie wektorów własnych dla alpha = 0 i alpha = 100
eigenvectors_alpha_0 = read_eigenvectors("data/eigenvectors_alpha_0.txt")
eigenvectors_alpha_100 = read_eigenvectors("data/eigenvectors_alpha_100.txt")


# Wykreślenie wektorów własnych dla alpha = 0
plt.figure(figsize=(10, 6))
for i in range(6):
    plt.plot(eigenvectors_alpha_0[i], label=f"Eigenvalue {i+1}")
plt.xlabel('Index')
plt.ylabel('Eigenvalue')
plt.title('Eigenvalues for alpha = 0')
plt.legend()
plt.grid(True)
plt.savefig('charts/eigenvectors_alpha_0.png')
plt.show()



# Wykreślenie wektorów własnych dla alpha = 100
plt.figure(figsize=(10, 6))
for i in range(6):
    plt.plot(eigenvectors_alpha_100[i], label=f"Eigenvalue {i+1}")
plt.xlabel('Index')
plt.ylabel('Eigenvalue')
plt.title('Eigenvalues for alpha = 100')
plt.legend()
plt.grid(True)
plt.savefig('charts/eigenvectors_alpha_100.png')
plt.show()
