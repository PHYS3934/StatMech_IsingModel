# Import libraries:
import numpy as np

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
