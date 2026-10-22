# Import libraries:
import numpy as np

def IsingEnergy(grid, J):
    """
    Compute the energy density of a spin configuration.
    """
    # Calculate the sum of nearest neighbors for each site
    neighbors = np.roll(grid, 1, axis=1) + np.roll(grid, -1, axis=1) + \
                np.roll(grid, 1, axis=0) + np.roll(grid, -1, axis=0)

    # Calculate total energy (each bond appears twice in the sum, hence the 1/2)
    energy = -J * np.sum(grid * neighbors) / (2 * grid.size)

    return energy
