# Import libraries:
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy import ndimage


def ClusterSizes(grid):
    """
    Find the clusters of a spin grid: connected regions of parallel spins.
    (Clusters are not joined across the periodic boundaries.)
    OUTPUTS: - clusterSizes: the number of spins in each cluster
             - sizeMap: an array like grid, giving the size of the cluster that each spin belongs to
    """
    clusterSizes = []
    sizeMap = np.zeros(grid.shape, dtype=int)

    for spin in [1, -1]:
        # Label the clusters of this spin direction: 0 for the other direction, then 1, 2, ...
        labels, _ = ndimage.label(grid == spin)
        # Count the spins with each label:
        sizes = np.bincount(labels.ravel())
        sizes[0] = 0  # label 0 is the other spin direction, not a cluster
        # Give each spin the size of its cluster:
        sizeMap = sizeMap + sizes[labels]
        clusterSizes.append(sizes[1:])

    return np.concatenate(clusterSizes), sizeMap


def ClusterSizeDistribution(clusterSizes, numSpins, numBins=20):
    """
    Bin the cluster sizes into logarithmically spaced bins from 1 to numSpins (the largest possible
    cluster). Using the same bins for every grid of a given size lets us compare grids directly.
    OUTPUTS: - binCenters: (geometric) center of each bin
             - probabilities: the fraction of clusters in each bin
    """
    # log10-spaced bin edges:
    binEdges = np.logspace(0, np.log10(numSpins + 1), numBins + 1)

    # Count the clusters in each bin, and convert counts to probabilities:
    counts, _ = np.histogram(clusterSizes, bins=binEdges)
    probabilities = counts / np.sum(counts)

    # Bin centers, as the geometric mean of the bin edges (the midpoint on a log scale):
    binCenters = np.sqrt(binEdges[1:] * binEdges[:-1])

    return binCenters, probabilities


def ClusterColors(grid, sizeMap):
    """
    Color each spin by the size of its cluster: up clusters from red (small) to yellow (large),
    down clusters from blue (small) to green (large). Neighboring clusters always point in
    opposite directions, so they always get different colors. Returns an RGB image.
    """
    # Position of each cluster size on a log scale, from 0 (one spin) to 1 (the whole grid):
    logSize = np.log(sizeMap) / np.log(grid.size)

    colorsUp = plt.cm.autumn(logSize)[:, :, :3]
    colorsDown = plt.cm.winter(logSize)[:, :, :3]
    return np.where((grid > 0)[:, :, np.newaxis], colorsUp, colorsDown)


def _plot_clusters(axes, grid, sizeMap, binCenters, probabilities):
    """
    Draw the three panels: the spin grid, its clusters colored by size, and the
    cluster-size distribution. Returns the plotted objects so an animation can update them.
    """
    ax_grid, ax_clusters, ax_dist = axes

    # The spin configuration (black = up, white = down):
    im_grid = ax_grid.imshow(grid, cmap='binary', vmin=-1, vmax=1)
    ax_grid.set_title('Spins')
    ax_grid.axis('off')

    # Each cluster colored by its direction and size:
    im_clusters = ax_clusters.imshow(ClusterColors(grid, sizeMap))
    ax_clusters.set_title('Clusters, colored by size')
    ax_clusters.set_xticks([])
    ax_clusters.set_yticks([])
    ax_clusters.set_xlabel('up: red (small) to yellow (large)\ndown: blue (small) to green (large)')

    # The distribution of cluster sizes (empty bins are left out of the log-log plot):
    show = probabilities > 0
    line, = ax_dist.loglog(binCenters[show], probabilities[show], 'o-k')
    ax_dist.set_xlim(0.7, grid.size * 1.5)
    ax_dist.set_xlabel('Cluster size')
    ax_dist.set_ylabel('Probability')
    ax_dist.set_title('Cluster-size distribution')
    ax_dist.set_box_aspect(1)

    return im_grid, im_clusters, line


def ClusterSizeStats(grid):
    """
    Get statistics on cluster sizes, where a cluster is a connected region of parallel spins.
    Plots the grid, its clusters, and the distribution of cluster sizes; returns the cluster sizes.
    """
    clusterSizes, sizeMap = ClusterSizes(grid)
    binCenters, probabilities = ClusterSizeDistribution(clusterSizes, grid.size)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5), layout='constrained')
    _plot_clusters(axes, grid, sizeMap, binCenters, probabilities)
    plt.show()

    return clusterSizes


def AnimateClusters(gridHistory, interval=200):
    """
    Animate the clusters and the cluster-size distribution of a sequence of grids
    (e.g., the gridHistory returned by SampleGrid). The distributions of earlier frames
    stay on the plot in light gray, so you can see the distribution settle down.
    Returns the animation (keep a reference to it, or it stops).
    """
    numSpins = gridHistory[0].size

    # Compute the clusters and their size distribution for every frame up front:
    frames = []
    for grid in gridHistory:
        clusterSizes, sizeMap = ClusterSizes(grid)
        binCenters, probabilities = ClusterSizeDistribution(clusterSizes, numSpins)
        frames.append((sizeMap, probabilities))

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5), layout='constrained')
    im_grid, im_clusters, line = _plot_clusters(axes, gridHistory[0], frames[0][0], binCenters, frames[0][1])

    # Fix the vertical axis so that it covers every frame:
    smallest = min(np.min(p[p > 0]) for _, p in frames)
    axes[2].set_ylim(smallest / 2, 1.5)

    def update(i):
        sizeMap, probabilities = frames[i]
        im_grid.set_array(gridHistory[i])
        im_clusters.set_array(ClusterColors(gridHistory[i], sizeMap))
        # Leave the previous distribution behind in light gray, then draw the current one:
        x, y = line.get_data()
        axes[2].loglog(x, y, '-', color='0.85', zorder=1)
        show = probabilities > 0
        line.set_data(binCenters[show], probabilities[show])
        line.set_zorder(2)
        axes[0].set_title(f'Spins (frame {i + 1} of {len(gridHistory)})')
        return im_grid, im_clusters, line

    ani = FuncAnimation(fig, update, frames=len(gridHistory), interval=interval, repeat=False)

    return ani


if __name__ == "__main__":
    # Example: watch clusters grow as a random grid evolves at the critical temperature
    from SampleGrid import SampleGrid
    N = 150
    kT_c = 2 / np.log(1 + np.sqrt(2))  # critical temperature
    grid = np.sign(0.5 - np.random.random((N, N)))
    _, _, _, gridHistory = SampleGrid(grid, kT_c, 1, 300 * N**2, 5 * N**2, "HeatBath")
    ani = AnimateClusters(gridHistory)
    plt.show()
