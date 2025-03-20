# Import libraries:
import numpy as np
from tqdm import tqdm

# Import necessary functions:
from myNeighbors import myNeighbors
from IsingEnergy import IsingEnergy
from WolffIteration import WolffIteration


def SampleGrid(grid, kT, J, numTimePoints, everyT, sampleHow="Metropolis", timeLag=0, saveVideo=False):
    """
    Sampling algorithms for the 2D Ising model
    Returns an animation-ready function and updates data arrays
    """
    N = grid.shape[0]
    grid_history = []  # Store grid states for animation
    
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
    num_samples = numTimePoints // everyT + 1  # +1 for initial state
    M_store = np.zeros(num_samples)
    energyStore = np.zeros(num_samples)
    
    # Calculate initial values
    M_store[0] = np.sum(grid) / grid.size
    energyStore[0] = IsingEnergy(grid, J)
    grid_history.append(grid.copy())
    
    # Define update function for a single step
    def update_step(grid, t, spin_idx=None):
        if sampleHow == "HeatBath":
            # Index, s, of the spin to consider flipping:
            s = spin_idx
            # Calculate the difference in energy between s up/down
            pUp = J * np.sum(grid.flat[adj[s-1]])
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
            deltaE = 2 * J * grid.flat[s-1] * np.sum(grid.flat[adj[s-1]])
            if deltaE < 0:
                # Always flip to lower energy
                grid.flat[s-1] = -grid.flat[s-1]
            else:
                # Calculate the transition probability
                p = np.exp(-deltaE/kT)
                # A transition to higher energy occurs with probability p:
                if np.random.random() <= p:
                    grid.flat[s-1] = -grid.flat[s-1]
        
        elif sampleHow == "Wolff":
            # Identify a cluster to flip using the Wolff algorithm
            p_wolff = 1 - np.exp(-2*J/kT)
            C, _ = WolffIteration(N, p_wolff, grid, adj)
            grid.flat[C] = -grid.flat[C]

            
        return grid
    
    # Run the simulation
    for t in tqdm(range(1, numTimePoints+1)):
        # Update the grid based on sampling method
        if sampleHow in ["HeatBath", "Metropolis"]:
            grid = update_step(grid, t, spin[t-1])
        elif sampleHow == "Wolff":
            grid = update_step(grid, t, None)
        
        # Store grid and observables at specified intervals
        if t % everyT == 0:
            idx = t // everyT
            # Calculate observables
            M = np.sum(grid) / grid.size
            E = IsingEnergy(grid, J)
            
            # Store values
            M_store[idx] = M
            energyStore[idx] = E
            grid_history.append(grid.copy())
    
    # Define animation update function
    def animate_func(i):
        return grid_history[i]
    
    return grid, energyStore, M_store, grid_history, animate_func