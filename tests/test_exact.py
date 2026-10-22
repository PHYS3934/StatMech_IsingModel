"""
Check the samplers against exact enumeration of a 4x4 periodic lattice, and the
correlation functions against a direct calculation. Each check runs on both copies
of the functions: the modules in python/ and jupyter/isingsim_funcs.py.
Run from the repo root with: python tests/test_exact.py
"""
# Import libraries:
import importlib.util
import itertools
import os
import sys
from types import SimpleNamespace
import numpy as np

repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')

# The functions in python/ (one module per function):
sys.path.insert(0, os.path.join(repo, 'python'))
from IsingEnergy import IsingEnergy
from myNeighbors import myNeighbors
from WolffIteration import WolffIteration
from SampleGrid import SampleGrid
from CorrelationFun import CorrelationFun
from RadialAverage import RadialAverage
python_funcs = SimpleNamespace(IsingEnergy=IsingEnergy, myNeighbors=myNeighbors,
                               WolffIteration=WolffIteration, SampleGrid=SampleGrid,
                               CorrelationFun=CorrelationFun, RadialAverage=RadialAverage)

# The same functions in jupyter/isingsim_funcs.py:
spec = importlib.util.spec_from_file_location('isingsim_funcs', os.path.join(repo, 'jupyter', 'isingsim_funcs.py'))
jupyter_funcs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(jupyter_funcs)

ALL_FUNCS = {'python': python_funcs, 'jupyter': jupyter_funcs}

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


def test_energy_ground_state(f):
    # All spins up: 2 bonds per spin, so the energy per spin is -2J
    assert f.IsingEnergy(np.ones((N, N)), J) == -2 * J


def test_wolff(f):
    np.random.seed(0)
    absM_exact, E_exact = exact_averages()
    adj = f.myNeighbors(range(0, N**2), N)
    p = 1 - np.exp(-2 * J / kT)
    grid = np.ones((N, N))
    absM = []
    E = []
    for t in range(60000):
        cluster, _ = f.WolffIteration(N, p, grid, adj)
        grid.flat[cluster] = -grid.flat[cluster]
        if t >= 1000:
            absM.append(np.abs(np.mean(grid)))
            E.append(f.IsingEnergy(grid, J))
    # Run-to-run standard deviations are about 0.002 (|m|) and 0.004 (energy);
    # the old per-site Wolff rule gave <|m|> = 0.59
    check('Wolff <|m|>', np.mean(absM), absM_exact, 0.01)
    check('Wolff energy', np.mean(E), E_exact, 0.025)


def test_single_spin(f):
    absM_exact, E_exact = exact_averages()
    for sampleHow in ['Metropolis', 'HeatBath']:
        np.random.seed(0)
        grid = np.ones((N, N))
        _, energies, M, _ = f.SampleGrid(grid, kT, J, 40000 * N**2, N**2, sampleHow)
        # Run-to-run standard deviations are about 0.003 (|m|) and 0.007 (energy)
        check(f'{sampleHow} <|m|>', np.mean(np.abs(M[100:])), absM_exact, 0.015)
        check(f'{sampleHow} energy', np.mean(energies[100:]), E_exact, 0.035)


def test_correlation(f):
    for n in [8, 9]:  # even and odd lattice sizes
        np.random.seed(0)
        A = np.sign(0.5 - np.random.random((n, n)))
        cor = f.CorrelationFun(A, doNorm=False)
        # Direct calculation, with zero displacement at [n//2, n//2]:
        for dx in range(-(n // 2), (n + 1) // 2):
            for dy in range(-(n // 2), (n + 1) // 2):
                direct = np.mean(A * np.roll(A, (-dx, -dy), axis=(0, 1)))
                assert np.isclose(cor[n // 2 + dx, n // 2 + dy], direct)
        # The radial average at distance 0 is <s^2> = 1:
        assert np.isclose(f.RadialAverage(cor, n)[0], 1)
    print('CorrelationFun and RadialAverage match the direct calculation')


if __name__ == '__main__':
    for name, f in ALL_FUNCS.items():
        print(f'--- {name}')
        test_energy_ground_state(f)
        test_correlation(f)
        test_wolff(f)
        test_single_spin(f)
    print('All tests passed')
