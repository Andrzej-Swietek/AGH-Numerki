import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters
N = 2000
d = 3
K = 10

# Step 1: Generate 3D vectors on the surface of a sphere using the Box-Muller method
u1 = np.random.uniform(0, 1, N)
u2 = np.random.uniform(0, 1, N)
u3 = np.random.uniform(0, 1, N)
u4 = np.random.uniform(0, 1, N)

x = np.sqrt(-2 * np.log(1 - u1)) * np.cos(2 * np.pi * u2)
y = np.sqrt(-2 * np.log(1 - u1)) * np.sin(2 * np.pi * u2)
z = np.sqrt(-2 * np.log(1 - u3)) * np.cos(2 * np.pi * u4)

# Normalize the vectors to lie on the surface of a sphere
r = np.sqrt(x**2 + y**2 + z**2)
x /= r
y /= r
z /= r

# Save vectors to a file
sphere_vectors = pd.DataFrame({'x': x, 'y': y, 'z': z})
sphere_vectors.to_csv('sphere_vectors_v2.csv', index=False)

# Plotting 3D points on a sphere
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, s=1)
ax.set_title('Points on the surface of a sphere')
plt.show()

# Step 2: Transform the distribution to be uniform within a sphere
ui = np.random.uniform(0, 1, N)
si = ui**(1/d)
ri_x = si * x
ri_y = si * y
ri_z = si * z

# Save the transformed vectors to a file
uniform_sphere_vectors = pd.DataFrame({'x': ri_x, 'y': ri_y, 'z': ri_z})
uniform_sphere_vectors.to_csv('uniform_sphere_vectors.csv', index=False)

# Plotting 3D points within a sphere
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(ri_x, ri_y, ri_z, s=1)
ax.set_title('Points uniformly distributed within a sphere')
plt.show()

# Step 3: Check the uniformity of the distribution within the sphere
radii = np.sqrt(ri_x**2 + ri_y**2 + ri_z**2)
delta = 1 / K
bins = np.linspace(0, 1, K + 1)
density, _ = np.histogram(radii, bins)

# Normalize density by the volume of the spherical shells
volumes = 4/3 * np.pi * (bins[1:]**3 - bins[:-1]**3)
normalized_density = density / volumes

# Plotting the density
plt.bar(range(1, K + 1), normalized_density, width=0.8)
plt.xlabel('Radial bin')
plt.ylabel('Density')
plt.title('Density of points in spherical shells')
plt.show()
