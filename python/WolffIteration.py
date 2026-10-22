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
        # Get all neighbors of all frontier spins (vectorized), keeping
        # duplicates: each entry is one bond, and each bond gets its own trial
        neighbors = adj[F].reshape(-1)

        # Only choose neighbors parallel to the seed spin and not already in the cluster
        Ci[C] = 1
        neighbors = neighbors[(grid.flat[neighbors] == s) & (Ci[neighbors] == 0)]

        # Activate each bond with probability p, then remove duplicate sites
        # (a site touching k cluster spins joins with probability 1 - (1-p)^k)
        r = np.random.random(len(neighbors))
        F = np.unique(neighbors[r < p])

        # Add to cluster (vectorized)
        C.extend(F)
    
    return C, i