import numpy as np

def WolffIteration(N, p, grid, adj):
    """
    Find a cluster, C, according the the Wolff sampling rule
    """
    # Random seed spin
    i = np.random.randint(0, N**2)
    
    # The cluster and frontier initialization
    C = [i]
    F = [i]
    s = grid.flat[i]  # seed spin direction
    
    # Indicator arrays to track cluster membership
    Ci = np.zeros(N**2, dtype=bool)
    Ci[i] = True
    
    while F:
        # Get all neighbors of frontier spins
        neighbors = adj[F].flatten()
        
        # Only choose spins parallel to the seed spin
        neighbors = [n for n in neighbors if grid.flat[n-1] == s]
        
        # Find elements not already in the cluster
        Fi = np.zeros(N**2, dtype=bool)
        Fi[np.array(neighbors)-1] = True  # Adjust for 0-indexing
        new_spins = np.where(Fi & ~Ci)[0]
        
        # Keep spins only with probability p
        F = []
        for spin in new_spins:
            if np.random.random() < p:
                F.append(spin)
                C.append(spin)
                Ci[spin] = True
    
    # Convert from 0-indexed to 1-indexed for consistency with MATLAB
    C = [c+1 for c in C]
    return C, i+1