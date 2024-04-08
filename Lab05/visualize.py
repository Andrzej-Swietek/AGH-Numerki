import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


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
eigenvectors_alpha_0 = read_eigenvectors("/content/data/eigenvectors_alpha_0.txt")
eigenvectors_alpha_100 = read_eigenvectors("/content/data/eigenvectors_alpha_100.txt")



# Wykreślenie wektorów własnych dla alpha = 0
plt.figure(figsize=(10, 6))
for i in range(6):
    plt.plot(eigenvectors_alpha_0[i], label=f"Eigenvalue {i+1}")
plt.xlabel('Index')
plt.ylabel('Eigenvalue')
plt.title('Eigenvalues for alpha = 0')
plt.legend()
plt.grid(True)
plt.savefig('eigenvectors_alpha_0.png')
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
plt.savefig('eigenvectors_alpha_100.png')
plt.show()





# Wykreślenie wektorów własnych dla alpha = 0
plt.figure(figsize=(10, 6))
for i in range(6):
    plt.plot(eigenvectors_alpha_0[i], label=f"Eigenvalue {i+1}")
plt.xlabel('Index', fontsize=14)
plt.ylabel('Eigenvalue', fontsize=14)
plt.title('Eigenvalues for alpha = 0', fontsize=16)
plt.legend(fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.grid(True)
plt.savefig('eigenvectors_alpha_0.png')
plt.show()




# Wykreślenie wektorów własnych dla alpha = 100
plt.figure(figsize=(10, 6))
for i in range(6):
    plt.plot(eigenvectors_alpha_100[i], label=f"Eigenvalue {i+1}")
plt.xlabel('Index', fontsize=14)
plt.ylabel('Eigenvalue', fontsize=14)
plt.title('Eigenvalues for alpha = 100', fontsize=16)
plt.legend(fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.grid(True)
plt.savefig('eigenvectors_alpha_100.png')
plt.show()





eigenvalues_data = pd.read_csv('data/eigenvalues.txt', delimiter='\t')
alpha_values = eigenvalues_data['alpha']
eigenvalues = eigenvalues_data.drop(columns=['alpha'])

# Calculate square root of eigenvalues -> omega
eigenvalues_sqrt = np.sqrt(eigenvalues)





plt.figure(figsize=(10, 6))
for i in range(len(eigenvalues.columns)):
    plt.plot(alpha_values, eigenvalues_sqrt.iloc[:, i], label=f"Eigenvalue {i+1}")

plt.xlabel('Alpha', fontsize=14)
plt.ylabel('Square Root of Eigenvalue (Omega)', fontsize=14)
plt.title('Square Root of Eigenvalues', fontsize=16)
plt.legend(fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.grid(True)
plt.savefig('eigenvectors_alpha_vs_omega.png')
plt.show()

