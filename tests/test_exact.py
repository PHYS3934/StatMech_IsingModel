"""
Check the samplers against exact enumeration of a 4x4 periodic lattice.
Run from the repo root with: python tests/test_exact.py
"""
# Import libraries:
import itertools
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'python'))
from IsingEnergy import IsingEnergy
from myNeighbors import myNeighbors
from WolffIteration import WolffIteration
from SampleGrid import SampleGrid

N = 4
J = 1
kT = 2.5


def exact_averages():
    """
    Exact <|m|> and energy per spin, summing over all 2^(N^2) states.
    """
    states = np.array(list(itertools.product([-1, 1], repeat=N**2))).reshape(-1, N, N)
    # Count each bond once (right and down neighbors only):
    E = -J * np.sum(states * np.roll(states, 1, axis=1) + states * np.roll(states, 1, axis=2), axis=(1, 2))
    w = np.exp(-(E - E.min()) / kT)
    w = w / np.sum(w)
    absM = np.abs(np.sum(states, axis=(1, 2))) / N**2
    return np.sum(w * absM), np.sum(w * E) / N**2


def check(name, x, x_exact, tol):
    print(f'{name}: {x:.4f} (exact {x_exact:.4f})')
    assert abs(x - x_exact) < tol, f'{name} differs from the exact value by more than {tol}'


def test_energy_ground_state():
    # All spins up: 2 bonds per spin, so the energy per spin is -2J
    assert IsingEnergy(np.ones((N, N)), J) == -2 * J


def test_wolff():
    np.random.seed(0)
    absM_exact, E_exact = exact_averages()
    adj = myNeighbors(range(0, N**2), N)
    p = 1 - np.exp(-2 * J / kT)
    grid = np.ones((N, N))
    absM = []
    E = []
    for t in range(60000):
        C, _ = WolffIteration(N, p, grid, adj)
        grid.flat[C] = -grid.flat[C]
        if t >= 1000:
            absM.append(np.abs(np.mean(grid)))
            E.append(IsingEnergy(grid, J))
    # Run-to-run standard deviations are about 0.002 (|m|) and 0.004 (energy);
    # the old per-site Wolff rule gave <|m|> = 0.59
    check('Wolff <|m|>', np.mean(absM), absM_exact, 0.01)
    check('Wolff energy', np.mean(E), E_exact, 0.025)


def test_single_spin():
    absM_exact, E_exact = exact_averages()
    for sampleHow in ['Metropolis', 'HeatBath']:
        np.random.seed(0)
        grid = np.ones((N, N))
        _, energies, M, _, _ = SampleGrid(grid, kT, J, 40000 * N**2, N**2, sampleHow)
        # Run-to-run standard deviations are about 0.003 (|m|) and 0.007 (energy)
        check(f'{sampleHow} <|m|>', np.mean(np.abs(M[100:])), absM_exact, 0.015)
        check(f'{sampleHow} energy', np.mean(energies[100:]), E_exact, 0.035)


if __name__ == '__main__':
    test_energy_ground_state()
    test_wolff()
    test_single_spin()
    print('All tests passed')
