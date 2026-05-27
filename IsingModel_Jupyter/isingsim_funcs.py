# Import this file to use the functions defined here

# Import necessary libraries
import numpy as np
from scipy import ndimage, signal
import matplotlib.pyplot as plt
from tqdm import tqdm

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
    s = np.array(s)
    adj = np.zeros((len(s), 4), dtype=int)

    # s = r*N + c
    r = s // N
    c = s % N

    # Calculate neighbor indices with periodic boundary conditions
    adj[:, 0] = np.mod(r + 1, N) * N + c     # down
    adj[:, 1] = np.mod(r - 1, N) * N + c     # up
    adj[:, 2] = r * N + np.mod(c + 1, N)     # right
    adj[:, 3] = r * N + np.mod(c - 1, N)     # left

    return adj

def WolffIteration(N, p, grid, adj):
    """
    Find a cluster according to the Wolff sampling rule - MATLAB-like implementation
    """
    import numpy as np

    # Random seed spin (0-indexed for Python)
    i = np.random.randint(0, N**2)

    # The cluster (initialized with seed)
    C = [i]

    # The frontier of spins
    F = [i]

    # Seed spin direction
    s = grid.flat[i]

    # Indicator function for cluster elements
    Ci = np.zeros(N**2, dtype=int)

    while len(F) > 0:
        # Get all neighbors of all frontier spins (vectorized)
        neighbors = adj[F].reshape(-1)  # Flattened array of all neighbors

        # Only choose neighbors parallel to the seed spin (vectorized)
        parallel_mask = [grid.flat[n] == s for n in neighbors]
        F = neighbors[parallel_mask]

        # Find elements that aren't in the cluster (using indicator arrays)
        Ci[C] = 1
        Fi = np.zeros(N**2, dtype=int)
        Fi[F] = 1

        # Elements in frontier but not in cluster
        F = np.where(Fi - Ci > 0)[0]

        # Apply probability filter (vectorized)
        r = np.random.random(len(F))
        F = F[r < p]

        # Add to cluster (vectorized)
        C.extend(F)

    return C, i

def _checkerboard_update(grid, kT, J, sampleHow, update_masks):
    for update in update_masks:
        neighbor_sum = (
            np.roll(grid, 1, axis=0)
            + np.roll(grid, -1, axis=0)
            + np.roll(grid, 1, axis=1)
            + np.roll(grid, -1, axis=1)
        )

        if sampleHow == "HeatBath":
            p_up = 1 / (1 + np.exp(-2 * J * neighbor_sum[update] / kT))
            grid[update] = np.where(np.random.random(np.sum(update)) < p_up, 1, -1)
        elif sampleHow == "Metropolis":
            spins = grid[update]
            deltaE = 2 * J * spins * neighbor_sum[update]
            flips = (deltaE < 0) | (np.random.random(len(deltaE)) < np.exp(-deltaE / kT))
            spins[flips] = -spins[flips]
            grid[update] = spins

    return grid

def SampleGrid(grid, kT, J, numTimePoints, everyT, sampleHow="Metropolis", timeLag=0, saveVideo=False):
    """
    Sampling algorithms for the 2D Ising model
    Returns an animation-ready function and updates data arrays
    """
    N = grid.shape[0]
    grid_history = []  # Store grid states for animation

    def store_sample(idx):
        M_store[idx] = np.sum(grid) / grid.size
        energyStore[idx] = IsingEnergy(grid, J)
        grid_history.append(grid.copy())

    sweep_size = N**2
    can_use_sweeps = (
        sampleHow in ["HeatBath", "Metropolis"]
        and numTimePoints % sweep_size == 0
        and everyT % sweep_size == 0
    )

    if can_use_sweeps:
        num_sweeps = numTimePoints // sweep_size
        sample_every_sweeps = everyT // sweep_size
        num_samples = num_sweeps // sample_every_sweeps + 1
        M_store = np.zeros(num_samples)
        energyStore = np.zeros(num_samples)
        store_sample(0)

        checkerboard = np.indices(grid.shape).sum(axis=0) % 2
        update_masks = [checkerboard == 0, checkerboard == 1]

        sample_idx = 1
        for sweep in tqdm(range(1, num_sweeps + 1)):
            grid = _checkerboard_update(grid, kT, J, sampleHow, update_masks)
            if sweep % sample_every_sweeps == 0:
                store_sample(sample_idx)
                sample_idx += 1

        def animate_func(i):
            return grid_history[i]

        return grid, energyStore, M_store, grid_history, animate_func

    # Precompute the indices adjacent to each spin index
    adj = myNeighbors(range(0, N**2), N)

    # Initialize based on sampling method
    if sampleHow in ["HeatBath", "Metropolis"]:
        # Precompute a sequence of random spins (with a linear index)
        spin = np.random.randint(0, N**2, numTimePoints)
    elif sampleHow == "Wolff":
        p = 1 - np.exp(-2*J/kT)
        spin = None  # We don't need spin for Wolff algorithm

    # Store for observables
    num_samples = numTimePoints // everyT + 1
    M_store = np.zeros(num_samples)
    energyStore = np.zeros(num_samples)

    # Calculate initial values
    store_sample(0)

    # Define update function for a single step
    def update_step(grid, t, spin_idx=None):
        if sampleHow == "HeatBath":
            # Index, s, of the spin to consider flipping:
            s = spin_idx
            # Calculate the difference in energy between s up/down
            pUp = J * np.sum(grid.flat[adj[s]])
            pDown = -pUp
            z = np.exp(-pUp/kT) + np.exp(-pDown/kT)
            p = np.exp(-pUp/kT) / z
            # Decide whether to set this spin up or down:
            if np.random.random() <= p:
                grid.flat[s] = -1
            else:
                grid.flat[s] = 1

        elif sampleHow == "Metropolis":
            # Index, s, of the spin to consider flipping:
            s = spin_idx
            # Compute the change in energy from flipping this spin:
            deltaE = 2 * J * grid.flat[s] * np.sum(grid.flat[adj[s]])
            if deltaE < 0:
                # Always flip to lower energy
                grid.flat[s] = -grid.flat[s]
            else:
                # Calculate the transition probability
                p = np.exp(-deltaE/kT)
                # A transition to higher energy occurs with probability p:
                if np.random.random() <= p:
                    grid.flat[s] = -grid.flat[s]

        elif sampleHow == "Wolff":
            # Identify a cluster to flip using the Wolff algorithm
            p_wolff = 1 - np.exp(-2*J/kT)
            C, _ = WolffIteration(N, p_wolff, grid, adj)
            grid.flat[C] = -grid.flat[C]


        return grid

    # Run the simulation
    sample_idx = 1
    for t in tqdm(range(0, numTimePoints)):
        # Update the grid based on sampling method
        if sampleHow in ["HeatBath", "Metropolis"]:
            grid = update_step(grid, t, spin[t])
        elif sampleHow == "Wolff":
            grid = update_step(grid, t, None)

        # Store grid and observables at specified intervals
        if (t + 1) % everyT == 0:
            store_sample(sample_idx)
            sample_idx += 1

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


def _sample_ising_grid(size=200, sweeps=250, kT=None, J=1, seed=0):
    """
    Generate a representative Ising grid for running this module as a script.
    """
    if kT is None:
        kT = 2 / np.log(1 + np.sqrt(2))

    rng = np.random.default_rng(seed)
    grid = rng.choice([-1, 1], size=(size, size))
    checkerboard = np.indices(grid.shape).sum(axis=0) % 2

    for _ in range(sweeps):
        for parity in (0, 1):
            neighbor_sum = (
                np.roll(grid, 1, axis=0)
                + np.roll(grid, -1, axis=0)
                + np.roll(grid, 1, axis=1)
                + np.roll(grid, -1, axis=1)
            )
            p_up = 1 / (1 + np.exp(-2 * J * neighbor_sum / kT))
            update = checkerboard == parity
            grid[update] = np.where(rng.random(np.sum(update)) < p_up[update], 1, -1)

    return grid


def ClusterSizeStats(grid):
    """
    Get statistics on cluster sizes
    -------------------------------------------------------------------------------
    """
    # Get connected components for positive values
    CC1, num_features_pos = ndimage.label(grid > 0)
    # Get area of each labeled region
    unique_pos, SS1 = np.unique(CC1[CC1 > 0], return_counts=True)

    # Get connected components for negative values
    CC2, num_features_neg = ndimage.label(grid < 0)
    # Get area of each labeled region
    unique_neg, SS2 = np.unique(CC2[CC2 > 0], return_counts=True)

    # Combine cluster areas
    clusterSizes = np.concatenate([SS1, SS2])
    if len(clusterSizes) == 0:
        raise ValueError("ClusterSizeStats requires a grid with positive or negative spins.")

    # -------------------------------------------------------------------------------
    # Logarithmic binning
    # -------------------------------------------------------------------------------
    numBins = min(len(np.unique(clusterSizes))//5, 12)
    numBins = max(numBins, 2)  # np.histogram needs at least two bin edges

    # log10-spaced bin edges:
    binEdges = np.logspace(np.log10(min(clusterSizes)*0.9999),
                          np.log10(max(clusterSizes)*1.0001),
                          numBins)

    # Bin the data using custom bin edges:
    N, binEdges = np.histogram(clusterSizes, bins=binEdges)

    # Bin centers as middle points between bin edges:
    binCenters = (binEdges[1:] + binEdges[:-1]) / 2

    # Convert counts to probabilities:
    Nnorm = N / np.sum(N)

    # -------------------------------------------------------------------------------
    # PLOTTING
    # -------------------------------------------------------------------------------
    plt.figure(figsize=(15, 5))

    plt.subplot(1, 3, 1)
    plt.imshow(grid, cmap='binary')
    plt.axis('image')

    plt.subplot(1, 3, 2)
    # Create colored labels for visualization
    L1 = CC1
    L2 = CC2

    # Color the positive clusters using 'jet' colormap
    RGB1 = np.zeros((*L1.shape, 3))
    for i in range(1, num_features_pos + 1):
        RGB1[L1 == i] = plt.cm.jet(i / num_features_pos)[:3]

    # Color the negative clusters using a different colormap
    lines = plt.get_cmap('tab20')
    RGB2 = np.zeros((*L2.shape, 3))
    for i in range(1, num_features_neg + 1):
        RGB2[L2 == i] = lines(i / num_features_neg)[:3]

    # Combine the two images
    plt.imshow(RGB1 + RGB2)
    plt.axis('image')

    plt.subplot(1, 3, 3)
    plt.loglog(binCenters, Nnorm, 'o-k')
    plt.xlabel('Cluster size')
    plt.ylabel('Probability')

    plt.tight_layout()
    plt.show()

    return clusterSizes


# def main():
#     grid = _sample_ising_grid()
#     ClusterSizeStats(grid)


# if __name__ == "__main__":
#     main()
