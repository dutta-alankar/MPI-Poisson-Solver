import numpy as np
import matplotlib.pyplot as plt

# Define grid parameters (must match C code)
GRID_SIZE = 64
BOX_SIZE = 2.0

# Load potential data
# The C code outputs a 2D slice (z = GRID_SIZE / 2) of x, y, potential
data = np.loadtxt('build/potential.dat')

# Reshape the potential data into a 2D grid
potential_2d = data[:, 2].reshape((GRID_SIZE, GRID_SIZE))

# Load density data
# The C code outputs a 2D slice (z = GRID_SIZE / 2) of x, y, density
density_data = np.loadtxt('build/density.dat')

# Reshape the density data into a 2D grid
density_2d = density_data[:, 2].reshape((GRID_SIZE, GRID_SIZE))

# Create a meshgrid for plotting
x = np.linspace(0, BOX_SIZE, GRID_SIZE)
y = np.linspace(0, BOX_SIZE, GRID_SIZE)
X, Y = np.meshgrid(x, y)

# Plotting
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

# Plot Density
im1 = axes[0].imshow(density_2d, origin='lower', extent=[0, BOX_SIZE, 0, BOX_SIZE], cmap='viridis')
axes[0].set_title('Density Field (2D Slice)')
axes[0].set_xlabel('X')
axes[0].set_ylabel('Y')
fig.colorbar(im1, ax=axes[0], label='Density')

# Plot Potential
im2 = axes[1].imshow(potential_2d, origin='lower', extent=[0, BOX_SIZE, 0, BOX_SIZE], cmap='plasma')
axes[1].set_title('Potential Field (2D Slice)')
axes[1].set_xlabel('X')
axes[1].set_ylabel('Y')
fig.colorbar(im2, ax=axes[1], label='Potential')

plt.tight_layout()
plt.savefig('density_potential_slice.png')
plt.show()
