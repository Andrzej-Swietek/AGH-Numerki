import numpy as np

# Define the matrix
A = np.array([[0.707107, 0.57735, 0.5, 0.447214, 0.408248, 0.377964, 0.353553],
              [0.57735, 0.707107, 0.57735, 0.5, 0.447214, 0.408248, 0.377964],
              [0.5, 0.57735, 0.707107, 0.57735, 0.5, 0.447214, 0.408248],
              [0.447214, 0.5, 0.57735, 0.707107, 0.57735, 0.5, 0.447214],
              [0.408248, 0.447214, 0.5, 0.57735, 0.707107, 0.57735, 0.5],
              [0.377964, 0.408248, 0.447214, 0.5, 0.57735, 0.707107, 0.57735],
              [0.353553, 0.377964, 0.408248, 0.447214, 0.5, 0.57735, 0.707107]])

# Calculate eigenvalues
eigenvalues = np.linalg.eigvals(A)
print("Eigenvalues:", eigenvalues)