# Import libraries:
import numpy as np

def IsingEnergy(grid, J):
    """
    Compute the energy per spin of a spin configuration,
        E/N^2 = -(J/N^2) * (sum over nearest-neighbor pairs <ij> of s_i*s_j),
    on an N x N lattice with periodic boundary conditions.
    """
    # Pair each spin with its right and down neighbors only, so that each bond is counted once
    # (np.roll wraps around the edges, which gives periodic boundary conditions)
    right = np.roll(grid, -1, axis=1)
    down = np.roll(grid, -1, axis=0)
    bondSum = np.sum(grid * right) + np.sum(grid * down)

    # Energy per spin:
    energy = -J * bondSum / grid.size

    return energy
