import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage


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


def main():
    grid = _sample_ising_grid()
    ClusterSizeStats(grid)


if __name__ == "__main__":
    main()
