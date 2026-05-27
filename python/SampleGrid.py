# Import libraries:
import numpy as np
from tqdm import tqdm

# Import necessary functions:
from myNeighbors import myNeighbors
from IsingEnergy import IsingEnergy
from WolffIteration import WolffIteration


def _checkerboard_update(grid, kT, J, sampleHow, update_masks):
    for update in update_masks:
        neighbor_sum = (
            np.roll(grid, 1, axis=0)
            + np.roll(grid, -1, axis=0)
            + np.roll(grid, 1, axis=1)
            + np.roll(grid, -1, axis=1)
        )

        if sampleHow == "HeatBath":
            p_up = 1 / (1 + np.exp(-2 * J * neighbor_sum[update] / kT))
            grid[update] = np.where(np.random.random(np.sum(update)) < p_up, 1, -1)
        elif sampleHow == "Metropolis":
            spins = grid[update]
            deltaE = 2 * J * spins * neighbor_sum[update]
            flips = (deltaE < 0) | (np.random.random(len(deltaE)) < np.exp(-deltaE / kT))
            spins[flips] = -spins[flips]
            grid[update] = spins

    return grid


def SampleGrid(grid, kT, J, numTimePoints, everyT, sampleHow="Metropolis", timeLag=0, saveVideo=False):
    """
    Sampling algorithms for the 2D Ising model
    Returns an animation-ready function and updates data arrays
    """
    N = grid.shape[0]
    grid_history = []  # Store grid states for animation

    def store_sample(idx):
        M_store[idx] = np.sum(grid) / grid.size
        energyStore[idx] = IsingEnergy(grid, J)
        grid_history.append(grid.copy())

    sweep_size = N**2
    can_use_sweeps = (
        sampleHow in ["HeatBath", "Metropolis"]
        and numTimePoints % sweep_size == 0
        and everyT % sweep_size == 0
    )

    if can_use_sweeps:
        num_sweeps = numTimePoints // sweep_size
        sample_every_sweeps = everyT // sweep_size
        num_samples = num_sweeps // sample_every_sweeps + 1
        M_store = np.zeros(num_samples)
        energyStore = np.zeros(num_samples)
        store_sample(0)

        checkerboard = np.indices(grid.shape).sum(axis=0) % 2
        update_masks = [checkerboard == 0, checkerboard == 1]

        sample_idx = 1
        for sweep in tqdm(range(1, num_sweeps + 1)):
            grid = _checkerboard_update(grid, kT, J, sampleHow, update_masks)
            if sweep % sample_every_sweeps == 0:
                store_sample(sample_idx)
                sample_idx += 1

        def animate_func(i):
            return grid_history[i]

        return grid, energyStore, M_store, grid_history, animate_func
    
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
    num_samples = numTimePoints // everyT + 1
    M_store = np.zeros(num_samples)
    energyStore = np.zeros(num_samples)
    
    # Calculate initial values
    store_sample(0)
    
    # Define update function for a single step
    def update_step(grid, t, spin_idx=None):
        if sampleHow == "HeatBath":
            # Index, s, of the spin to consider flipping:
            s = spin_idx
            # Calculate the difference in energy between s up/down
            pUp = J * np.sum(grid.flat[adj[s]])
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
            deltaE = 2 * J * grid.flat[s] * np.sum(grid.flat[adj[s]])
            if deltaE < 0:
                # Always flip to lower energy
                grid.flat[s] = -grid.flat[s]
            else:
                # Calculate the transition probability
                p = np.exp(-deltaE/kT)
                # A transition to higher energy occurs with probability p:
                if np.random.random() <= p:
                    grid.flat[s] = -grid.flat[s]
        
        elif sampleHow == "Wolff":
            # Identify a cluster to flip using the Wolff algorithm
            p_wolff = 1 - np.exp(-2*J/kT)
            C, _ = WolffIteration(N, p_wolff, grid, adj)
            grid.flat[C] = -grid.flat[C]

            
        return grid
    
    # Run the simulation
    sample_idx = 1
    for t in tqdm(range(0, numTimePoints)):
        # Update the grid based on sampling method
        if sampleHow in ["HeatBath", "Metropolis"]:
            grid = update_step(grid, t, spin[t])
        elif sampleHow == "Wolff":
            grid = update_step(grid, t, None)
        
        # Store grid and observables at specified intervals
        if (t + 1) % everyT == 0:
            store_sample(sample_idx)
            sample_idx += 1
    
    # Define animation update function
    def animate_func(i):
        return grid_history[i]
    
    return grid, energyStore, M_store, grid_history, animate_func
