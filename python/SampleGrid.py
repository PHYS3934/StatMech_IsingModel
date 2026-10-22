# Import libraries:
import numpy as np
from tqdm import tqdm

# Import necessary functions:
from myNeighbors import myNeighbors
from IsingEnergy import IsingEnergy
from WolffIteration import WolffIteration


def SampleGrid(grid, kT, J, numTimePoints, everyT, sampleHow="Metropolis"):
    """
    Evolve a 2D Ising spin configuration using a Markov chain Monte Carlo sampler.

    INPUTS: - grid: N x N array of spins (+1 or -1); it is updated in place
            - kT: temperature (in units where Boltzmann's constant is 1)
            - J: coupling strength
            - numTimePoints: number of update steps. For 'HeatBath' and 'Metropolis',
                  each step updates one spin (so N^2 steps make one sweep of the lattice).
                  For 'Wolff', each step flips one cluster.
            - everyT: store the energy, magnetization, and grid every everyT steps
            - sampleHow: 'HeatBath', 'Metropolis', or 'Wolff'
    OUTPUTS: - grid: the final spin configuration
             - energies: energy per spin at each stored time
             - magnetizations: magnetization per spin at each stored time
             - gridHistory: list of the stored spin configurations (for animation)
    """
    if sampleHow not in ["HeatBath", "Metropolis", "Wolff"]:
        raise ValueError(f"Unknown sampling method '{sampleHow}'")

    N = grid.shape[0]

    # Precompute the linear indices of the four neighbors of each spin:
    adj = myNeighbors(range(N**2), N)

    # Probability of activating a bond in the Wolff algorithm:
    pBond = 1 - np.exp(-2 * J / kT)

    # Store the energy, magnetization, and grid, starting with the initial configuration:
    energies = [IsingEnergy(grid, J)]
    magnetizations = [np.mean(grid)]
    gridHistory = [grid.copy()]

    # Evolve the Markov chain for numTimePoints steps:
    for t in tqdm(range(numTimePoints)):
        if sampleHow == "HeatBath":
            # Pick a random spin, s (as a linear index into grid.flat):
            s = np.random.randint(N**2)
            # Sum of its four neighbors:
            h = np.sum(grid.flat[adj[s]])
            # Energy of spin s pointing up or down, given its neighbors:
            E_up = -J * h
            E_down = J * h
            # Boltzmann probability that spin s points up:
            P_up = np.exp(-E_up / kT) / (np.exp(-E_up / kT) + np.exp(-E_down / kT))
            # Set spin s to up with probability P_up (regardless of its current value):
            if np.random.random() < P_up:
                grid.flat[s] = 1
            else:
                grid.flat[s] = -1

        elif sampleHow == "Metropolis":
            # Pick a random spin, s (as a linear index into grid.flat):
            s = np.random.randint(N**2)
            # Change in energy if spin s were flipped:
            deltaE = 2 * J * grid.flat[s] * np.sum(grid.flat[adj[s]])
            # Accept the flip with probability min(1, exp(-deltaE/kT)),
            # i.e., always flip if the energy does not increase:
            if deltaE <= 0 or np.random.random() < np.exp(-deltaE / kT):
                grid.flat[s] = -grid.flat[s]

        elif sampleHow == "Wolff":
            # Grow a cluster using the Wolff rule, then flip all of its spins:
            cluster, _ = WolffIteration(N, pBond, grid, adj)
            grid.flat[cluster] = -grid.flat[cluster]

        # Store the energy, magnetization, and grid every everyT steps:
        if (t + 1) % everyT == 0:
            energies.append(IsingEnergy(grid, J))
            magnetizations.append(np.mean(grid))
            gridHistory.append(grid.copy())

    return grid, np.array(energies), np.array(magnetizations), gridHistory
