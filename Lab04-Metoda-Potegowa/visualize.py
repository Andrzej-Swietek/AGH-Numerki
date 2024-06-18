import matplotlib.pyplot as plt

# Data
eigenvalues = [
    [3.58916, 3.59582, 3.59586, 3.59586, 3.59586, 3.59586, 3.59586, 3.59586, 3.59586, 3.59586, 3.59586, 3.59586],
    [0.000548461, 0.282512, 0.284531, 0.284903, 0.284972, 0.284985, 0.284987, 0.284988, 0.284988, 0.284988, 0.284988, 0.284988],
    [1.91485e-05, 0.121709, 0.122242, 0.122514, 0.12265, 0.122719, 0.122753, 0.12277, 0.122778, 0.122782, 0.122784, 0.122785],
    [7.64792e-07, 0.0865954, 0.0865954, 0.0865955, 0.0865995, 0.0867844, 0.0952313, 0.312156, 0.577366, 0.590103, 0.590384, 0.59039],
    [7.64792e-07, 0.0865955, 0.0865955, 0.0865955, 0.0865955, 0.0865955, 0.0865955, 0.0865955, 0.0865955, 0.0865955, 0.0865955, 0.0865955],
    [4.55772e-11, 3.06401e-07, 0.170944, 0.170974, 0.170974, 0.170974, 0.170974, 0.170974, 0.170974, 0.170974, 0.170974, 0.170974],
    [4.55772e-11, 3.03351e-07, 1.50615e-06, 0.0981543, 0.0981544, 0.0981544, 0.0981544, 0.0981544, 0.0981544, 0.0981544, 0.0981544, 0.0981544]
]

# Predefined list of colors
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

# # Plotting and Saving
for i in range(len(eigenvalues)):
    plt.figure(figsize=(8, 6))
    plt.plot(range(1, len(eigenvalues[i]) + 1), eigenvalues[i], marker='o', color=colors[i % len(colors)])  # Cycle through colors list
    plt.title(f"Change of Eigenvalue $\lambda_{i+1}$ with Iterations ", fontsize=16)
    plt.xlabel("Iterations", fontsize=14)
    plt.ylabel("Eigenvalues", fontsize=14)
    plt.grid(True)
    plt.xticks(fontsize=12)  # Adjust font size of x-axis numbers
    plt.yticks(fontsize=12)  # Adjust font size of y-axis numbers
    plt.tight_layout()

    # Save plot as PNG
    plt.savefig(f"eigenvalue_chart_{i+1}.png")
    plt.show()

plt.figure(figsize=(10, 8))
for i, eigenvalue in enumerate(eigenvalues):
    plt.plot(range(1, len(eigenvalue) + 1), eigenvalue, marker='o', color=colors[i % len(colors)], label=f"$\lambda_{i+1}$ ")
plt.title("Eigenvalues vs Iterations (All Charts)", fontsize=16)
plt.xlabel("Iterations", fontsize=14)
plt.ylabel("Eigenvalues", fontsize=14)
plt.grid(True)
plt.legend(fontsize=12)
plt.xticks(fontsize=12)  # Adjust font size of x-axis numbers
plt.yticks(fontsize=12)  # Adjust font size of y-axis numbers
plt.tight_layout()

# Save plot as PNG
plt.savefig("eigenvalues_all_charts.png")
plt.show()