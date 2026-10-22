# Helper functions for the IsingSim and WolffDemo notebooks.
# These are the same functions as in ../python/ (one file per function there);
# keep the two copies in sync (tests/test_exact.py checks both).

# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from tqdm import tqdm


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
    Grow one cluster of spins according to the Wolff rule.

    Starting from a random seed spin, the cluster grows outwards one shell at a time.
    Each bond between a cluster spin and a parallel neighbor outside the cluster is
    activated with probability p = 1 - exp(-2J/kT), and a neighbor joins the cluster
    if at least one of its bonds to the cluster is activated.

    INPUTS: - N: linear size of the lattice (N x N spins)
            - p: probability of activating a bond
            - grid: the N x N spin configuration (not modified here)
            - adj: N^2 x 4 array of the neighbors of each spin (from myNeighbors)
    OUTPUTS: - cluster: linear indices of the spins in the cluster
             - seed: linear index of the seed spin
    """
    # Pick a random seed spin (as a linear index into grid.flat):
    seed = np.random.randint(0, N**2)
    seedSpin = grid.flat[seed]  # the cluster contains only spins pointing this way

    # Keep track of which spins are in the cluster:
    inCluster = np.zeros(N**2, dtype=bool)
    inCluster[seed] = True

    # The frontier: spins that joined in the last step, whose bonds have not been tried yet
    frontier = np.array([seed])

    while len(frontier) > 0:
        # All bonds from the frontier to its neighbors. There is one entry per bond, so a
        # neighbor touching two frontier spins appears twice and gets two chances to join:
        neighbors = adj[frontier].reshape(-1)

        # Only bonds to parallel spins that are not already in the cluster can be activated:
        candidates = neighbors[(grid.flat[neighbors] == seedSpin) & ~inCluster[neighbors]]

        # Activate each bond with probability p:
        activated = candidates[np.random.random(len(candidates)) < p]

        # The newly joined spins (each counted once) form the next frontier:
        frontier = np.unique(activated)
        inCluster[frontier] = True

    cluster = np.flatnonzero(inCluster)

    return cluster, seed


def SampleGrid(grid, kT, J, numTimePoints, everyT, sampleHow="Metropolis"):
    """
    Evolve a 2D Ising spin configuration using a Markov chain Monte Carlo sampler.

    INPUTS: - grid: N x N array of spins (+1 or -1); it is updated in place
            - kT: temperature (in units where Boltzmann's constant is 1)
            - J: coupling strength
            - numTimePoints: number of update steps. For 'HeatBath' and 'Metropolis',
                  each step updates one spin (so N^2 steps make one sweep of the lattice).
                  For 'Wolff', each step flips one cluster.
            - everyT: store the energy, magnetization, and grid every everyT steps
            - sampleHow: 'HeatBath', 'Metropolis', or 'Wolff'
    OUTPUTS: - grid: the final spin configuration
             - energies: energy per spin at each stored time
             - magnetizations: magnetization per spin at each stored time
             - gridHistory: list of the stored spin configurations (for animation)
    """
    if sampleHow not in ["HeatBath", "Metropolis", "Wolff"]:
        raise ValueError(f"Unknown sampling method '{sampleHow}'")

    N = grid.shape[0]

    # Precompute the linear indices of the four neighbors of each spin:
    adj = myNeighbors(range(N**2), N)

    # Probability of activating a bond in the Wolff algorithm:
    pBond = 1 - np.exp(-2 * J / kT)

    # Store the energy, magnetization, and grid, starting with the initial configuration:
    energies = [IsingEnergy(grid, J)]
    magnetizations = [np.mean(grid)]
    gridHistory = [grid.copy()]

    # Evolve the Markov chain for numTimePoints steps:
    for t in tqdm(range(numTimePoints)):
        if sampleHow == "HeatBath":
            # Pick a random spin, s (as a linear index into grid.flat):
            s = np.random.randint(N**2)
            # Sum of its four neighbors:
            h = np.sum(grid.flat[adj[s]])
            # Energy of spin s pointing up or down, given its neighbors:
            E_up = -J * h
            E_down = J * h
            # Boltzmann probability that spin s points up:
            P_up = np.exp(-E_up / kT) / (np.exp(-E_up / kT) + np.exp(-E_down / kT))
            # Set spin s to up with probability P_up (regardless of its current value):
            if np.random.random() < P_up:
                grid.flat[s] = 1
            else:
                grid.flat[s] = -1

        elif sampleHow == "Metropolis":
            # Pick a random spin, s (as a linear index into grid.flat):
            s = np.random.randint(N**2)
            # Change in energy if spin s were flipped:
            deltaE = 2 * J * grid.flat[s] * np.sum(grid.flat[adj[s]])
            # Accept the flip with probability min(1, exp(-deltaE/kT)),
            # i.e., always flip if the energy does not increase:
            if deltaE <= 0 or np.random.random() < np.exp(-deltaE / kT):
                grid.flat[s] = -grid.flat[s]

        elif sampleHow == "Wolff":
            # Grow a cluster using the Wolff rule, then flip all of its spins:
            cluster, _ = WolffIteration(N, pBond, grid, adj)
            grid.flat[cluster] = -grid.flat[cluster]

        # Store the energy, magnetization, and grid every everyT steps:
        if (t + 1) % everyT == 0:
            energies.append(IsingEnergy(grid, J))
            magnetizations.append(np.mean(grid))
            gridHistory.append(grid.copy())

    return grid, np.array(energies), np.array(magnetizations), gridHistory


def CorrelationFun(A, doNorm=True):
    """
    Computes the spin-spin correlation function of a square spin grid, A, with periodic boundaries:
        C(dx, dy) = mean over all sites (x, y) of A[x, y] * A[x + dx, y + dy]
    INPUTS: - A: a square (spin) grid.
            - doNorm: (boolean), whether to subtract the squared lattice-average magnetization,
                      giving the connected correlation function
    The output is centered: zero displacement (dx = dy = 0) is at index [N//2, N//2].
    """
    # The straightforward way loops over every displacement, shifting the grid each time:
    #     cor[dx, dy] = np.mean(A * np.roll(A, (-dx, -dy), axis=(0, 1)))
    # The Fourier transform gives the same result for all displacements at once, because the
    # autocorrelation is the inverse transform of the power spectrum, |FFT(A)|^2
    # (the Wiener-Khinchin theorem; the FFT assumes periodic boundaries, as our lattice has):
    cor = np.real(np.fft.ifft2(np.abs(np.fft.fft2(A))**2)) / A.size

    if doNorm:
        cor = cor - np.mean(A)**2

    # Move zero displacement from the corner [0, 0] to the center [N//2, N//2]:
    cor = np.fft.fftshift(cor)

    return cor


def RadialAverage(cor, N):
    """
    Compute the radial average of the NxN correlation function, cor
    (average out the angular dependence from a 2D connected correlation function).
    Zero displacement must be at the center, cor[N//2, N//2], as returned by CorrelationFun.
    R[r] is the average of cor over all displacements at distance r (to the nearest integer).
    """
    # Distance of each displacement (dx, dy) from zero displacement at the center:
    d = np.arange(N) - N // 2
    dx, dy = np.meshgrid(d, d)
    distance = np.sqrt(dx**2 + dy**2)

    # Only go out to about N/2: beyond that, displacements wrap around the periodic boundaries
    numR = (N + 1) // 2
    R = np.zeros(numR)

    for r in range(numR):
        # Average over the ring of displacements with r - 0.5 <= distance < r + 0.5:
        ring = (r - 0.5 <= distance) & (distance < r + 0.5)
        R[r] = np.mean(cor[ring])

    return R


def ClusterSizeStats(grid):
    """
    Get statistics on cluster sizes, where a cluster is a connected region of parallel spins.
    Plots the grid, its clusters, and the distribution of cluster sizes; returns the cluster sizes.
    -------------------------------------------------------------------------------
    """
    # Label the connected regions of up spins (labelsUp is 0 for down spins, 1, 2, ... for each cluster):
    labelsUp, numClustersUp = ndimage.label(grid > 0)
    # Count the spins in each cluster:
    _, sizesUp = np.unique(labelsUp[labelsUp > 0], return_counts=True)

    # Same for the down spins:
    labelsDown, numClustersDown = ndimage.label(grid < 0)
    _, sizesDown = np.unique(labelsDown[labelsDown > 0], return_counts=True)

    # Combine cluster sizes
    clusterSizes = np.concatenate([sizesUp, sizesDown])

    # -------------------------------------------------------------------------------
    # Logarithmic binning
    # -------------------------------------------------------------------------------
    numBins = min(len(np.unique(clusterSizes))//5, 12)
    numBins = max(numBins, 2)  # np.histogram needs at least two bin edges

    # log10-spaced bin edges:
    binEdges = np.logspace(np.log10(min(clusterSizes)*0.9999),
                           np.log10(max(clusterSizes)*1.0001),
                           numBins)

    # Count the clusters in each bin:
    counts, binEdges = np.histogram(clusterSizes, bins=binEdges)

    # Bin centers as middle points between bin edges:
    binCenters = (binEdges[1:] + binEdges[:-1]) / 2

    # Convert counts to probabilities:
    probabilities = counts / np.sum(counts)

    # -------------------------------------------------------------------------------
    # PLOTTING
    # -------------------------------------------------------------------------------
    plt.figure(figsize=(15, 5))

    # The spin configuration:
    plt.subplot(1, 3, 1)
    plt.imshow(grid, cmap='binary')
    plt.axis('image')

    # Each cluster in its own color (up clusters from 'jet', down clusters from 'tab20'):
    plt.subplot(1, 3, 2)
    colorsUp = np.zeros((*labelsUp.shape, 3))
    for i in range(1, numClustersUp + 1):
        colorsUp[labelsUp == i] = plt.cm.jet(i / numClustersUp)[:3]
    colorsDown = np.zeros((*labelsDown.shape, 3))
    for i in range(1, numClustersDown + 1):
        colorsDown[labelsDown == i] = plt.cm.tab20(i / numClustersDown)[:3]
    plt.imshow(colorsUp + colorsDown)
    plt.axis('image')

    # The distribution of cluster sizes:
    plt.subplot(1, 3, 3)
    plt.loglog(binCenters, probabilities, 'o-k')
    plt.xlabel('Cluster size')
    plt.ylabel('Probability')

    plt.tight_layout()
    plt.show()

    return clusterSizes
