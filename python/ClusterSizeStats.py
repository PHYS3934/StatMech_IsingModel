# Import libraries:
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage


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


if __name__ == "__main__":
    # Example: cluster statistics of a grid sampled at the critical temperature
    from SampleGrid import SampleGrid
    N = 100
    kT = 2 / np.log(1 + np.sqrt(2))  # critical temperature
    grid = np.sign(0.5 - np.random.random((N, N)))
    grid, _, _, _ = SampleGrid(grid, kT, 1, 200 * N**2, N**2, "HeatBath")
    ClusterSizeStats(grid)
