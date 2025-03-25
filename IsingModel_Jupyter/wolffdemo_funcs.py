# Import necessary libraries
import numpy as np
from scipy import signal
from scipy import ndimage
import matplotlib.pyplot as plt 

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
    
    # -------------------------------------------------------------------------------
    # Logarithmic binning
    # -------------------------------------------------------------------------------
    numBins = min(len(np.unique(clusterSizes))//5, 12)
    numBins = max(numBins, 1)  # Ensure at least one bin
    
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
    plt.axis('square')
    
    plt.subplot(1, 3, 2)
    # Create colored labels for visualization
    L1 = CC1
    L2 = CC2
    
    # Color the positive clusters using 'jet' colormap
    RGB1 = np.zeros((*L1.shape, 3))
    for i in range(1, num_features_pos + 1):
        RGB1[L1 == i] = plt.cm.jet(i / num_features_pos)[:3]
    
    # Color the negative clusters using a different colormap
    RGB2 = np.zeros((*L2.shape, 3))
    for i in range(1, num_features_neg + 1):
        RGB2[L2 == i] = plt.cm.lines(i / num_features_neg)[:3]
        
    # Combine the two images
    plt.imshow(RGB1 + RGB2)
    plt.axis('square')
    
    plt.subplot(1, 3, 3)
    plt.loglog(binCenters, Nnorm, 'o-k')
    plt.xlabel('Cluster size')
    plt.ylabel('Probability')
    
#     plt.tight_layout()
    plt.show()
    
    return clusterSizes

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
