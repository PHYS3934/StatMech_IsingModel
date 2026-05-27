import numpy as np

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