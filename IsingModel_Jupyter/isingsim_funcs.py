# Import this file to use the functions defined here

# Import necessary libraries
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

def IsingEnergy(grid, J):
    """
    Compute the energy density of a spin configuration.
    """
    # Calculate the sum of nearest neighbors for each site
    neighbors = np.roll(grid, 1, axis=1) + np.roll(grid, -1, axis=1) + \
                np.roll(grid, 1, axis=0) + np.roll(grid, -1, axis=0)

    # Calculate total energy
    energy = -J * np.sum(grid * neighbors) / grid.size

    return energy

def myNeighbors(s, N):
    """
    Take a list of linear indices s and return the linear indices of the
    neighbors of s on an N by N grid with periodic boundary conditions.
    """
    s = np.array(s) - 1  # index by zero
    adj = np.zeros((len(s), 4), dtype=int)

    # s = r*N + c
    r = s // N
    c = s % N

    # Calculate neighbor indices with periodic boundary conditions
    adj[:, 0] = np.mod(r + 1, N) * N + c     # down
    adj[:, 1] = np.mod(r - 1, N) * N + c     # up
    adj[:, 2] = r * N + np.mod(c + 1, N)     # right
    adj[:, 3] = r * N + np.mod(c - 1, N)     # left

    adj = adj + 1  # index by one again

    return adj

def SampleGrid(grid, kT, J, numTimePoints, everyT, sampleHow="Metropolis", timeLag=0, saveVideo=False):
    """
    Sampling algorithms for the 2D Ising model
    Returns an animation-ready function and updates data arrays
    """
    N = grid.shape[0]
    grid_history = []  # Store grid states for animation

    # Precompute the indices adjacent to each spin index
    adj = myNeighbors(range(1, N**2+1), N)

    # Initialize based on sampling method
    if sampleHow in ["HeatBath", "Metropolis"]:
        # Precompute a sequence of random spins (with a linear index)
        spin = np.random.randint(1, N**2+1, numTimePoints)
    elif sampleHow == "Wolff":
        p = 1 - np.exp(-2*J/kT)

    # Store for observables
    num_samples = numTimePoints // everyT + 1  # +1 for initial state
    M_store = np.zeros(num_samples)
    energyStore = np.zeros(num_samples)

    # Calculate initial values
    M_store[0] = np.sum(grid) / grid.size
    energyStore[0] = IsingEnergy(grid, J)
    grid_history.append(grid.copy())

    # Define update function for a single step
    def update_step(grid, t, spin=None, sampleHow=sampleHow):
        if sampleHow == "HeatBath":
            # Index, s, of the spin to consider flipping:
            s = spin[t-1]
            # Calculate the difference in energy between s up/down
            pUp = J * np.sum(grid.flat[adj[s-1]-1])
            pDown = -pUp
            z = np.exp(-pUp/kT) + np.exp(-pDown/kT)
            p = np.exp(-pUp/kT) / z
            # Decide whether to set this spin up or down:
            if np.random.random() <= p:
                grid.flat[s-1] = -1
            else:
                grid.flat[s-1] = 1

        elif sampleHow == "Metropolis":
            # Index, s, of the spin to consider flipping:
            s = spin[t-1]
            # Compute the change in energy from flipping this spin:
            deltaE = 2 * J * grid.flat[s-1] * np.sum(grid.flat[adj[s-1]-1])
            if deltaE < 0:
                # Always flip to lower energy
                grid.flat[s-1] = -grid.flat[s-1]
            else:
                # Calculate the transition probability
                p = np.exp(-deltaE/kT)
                # A transition to higher energy occurs with probability p:
                if np.random.random() <= p:
                    grid.flat[s-1] = -grid.flat[s-1]

        elif sampleHow == "Wolff":
            # Identify a cluster to flip using the Wolff algorithm
            p_wolff = 1 - np.exp(-2*J/kT)
            C, _ = WolffIteration(N, p_wolff, grid, adj)
            # Convert to 0-indexed for numpy array indexing
            C_idx = [c-1 for c in C]
            # Flip the cluster
            grid.flat[C_idx] = -grid.flat[C_idx]

        return grid

    # Run the simulation
    for t in range(1, numTimePoints+1):
        # Update the grid
        grid = update_step(grid, t, spin)

        # Store grid and observables at specified intervals
        if t % everyT == 0:
            idx = t // everyT
            # Calculate observables
            M = np.sum(grid) / grid.size
            E = IsingEnergy(grid, J)

            # Store values
            M_store[idx] = M
            energyStore[idx] = E
            grid_history.append(grid.copy())

    # Define animation update function
    def animate_func(i):
        return grid_history[i]

    return grid, energyStore, M_store, grid_history, animate_func

def CorrelationFun(A, doNorm=True):
    """
    Computes the correlation function of A
    INPUTS: - A: a square (spin) grid.
            - doNorm: (boolean), whether to subtract the lattice-average M squared
    Output is centered at size(A)/2 because of periodic boundary conditions
    """
    # Calculate correlation using convolution
    # Extend A by tiling a 180-degree rotated version to handle periodicity
    A_rot = np.rot90(A, 2)
    A_extended = np.tile(A_rot, (2, 2))

    # Use scipy's signal.convolve2d for the convolution
    cor = signal.convolve2d(A, A_extended, mode='same') / A.size

    if doNorm:
        cor = cor - np.mean(A)**2

    # Center the correlation function
    cor = np.roll(cor, (A.shape[0]//2, A.shape[1]//2), axis=(0, 1))

    return cor

def RadialAverage(cor, N):
    """
    Compute the radial average of the NxN correlation function, cor
    (average out the angular dependence from a 2D connected correlation function)
    """
    L = N // 2 + (N % 2)  # ceiling of N/2
    c = np.arange(1, N+1) - L
    X, Y = np.meshgrid(c, c)
    rho = np.sqrt(X**2 + Y**2)  # equivalent to cart2pol in MATLAB

    rbins = np.arange(-0.5, L+0.5)  # only go to N/2, because of periodic boundaries
    R = np.zeros(L)

    for j in range(L):
        r = rbins[j]
        ring = (r <= rho) & (rho < r+1)
        if np.sum(ring) > 0:  # Avoid division by zero
            R[j] = np.sum(cor[ring]) / np.sum(ring)

    # Replace NaN values with 0
    R[np.isnan(R)] = 0

    return R
