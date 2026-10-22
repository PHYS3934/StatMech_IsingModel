# Simulate the 2D Ising model, H = -J * (sum over nearest-neighbor pairs <ij> of s_i*s_j),
# on an N x N lattice with periodic boundary conditions, then plot an animation,
# spin-spin correlations, and the energy and magnetization over time.

# Import necessary libraries:
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import mplcursors

# Import necessary functions:
from SampleGrid import SampleGrid
from RadialAverage import RadialAverage
from CorrelationFun import CorrelationFun

# -------------------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------------------
N = 100  # linear lattice size (N x N spins)
J = 1  # coupling strength
kT_c = 2 / np.log(1 + np.sqrt(2))  # critical temperature, about 2.269 (for J = 1)
kT = 0.85 * kT_c  # temperature (here below the critical temperature)

fractionUp = 0.5  # average fraction of +1 spins in the initial configuration
samplingMethod = 'Metropolis'  # 'HeatBath', 'Metropolis', or 'Wolff'

# Number of steps to run (numTimePoints), and how often to store the grid (everyT):
if samplingMethod in ['Metropolis', 'HeatBath']:
    # Each step updates one spin, so N^2 steps make one sweep of the lattice
    numTimePoints = 100 * N**2  # 100 sweeps
    everyT = N**2  # store once per sweep
    stepsPerTimeUnit, timeUnit = N**2, 'sweeps'
elif samplingMethod == 'Wolff':
    # Each step flips one cluster
    # (starting from a random grid, clusters are small at first and grow as the grid orders)
    numTimePoints = 20 * N
    everyT = N // 5
    stepsPerTimeUnit, timeUnit = 1, 'cluster flips'

# -------------------------------------------------------------------------------
# Run the simulation
# -------------------------------------------------------------------------------
# Random initial configuration: each spin is +1 with probability fractionUp
grid = np.sign(fractionUp - np.random.random((N, N)))

finalGrid, energies, magnetizations, gridHistory = SampleGrid(
    grid, kT, J, numTimePoints, everyT, samplingMethod)

# Time of each stored sample:
times = np.arange(len(energies)) * everyT / stepsPerTimeUnit

# -------------------------------------------------------------------------------
# Animation of the stored spin configurations
# -------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 8))

im = ax.imshow(gridHistory[0], cmap='binary', vmin=-1, vmax=1)
ax.set_title(f'2D Ising Model with {N} by {N} lattice', size=15)
plt.axis('off')

# Text annotation for energy and magnetization
time_text = ax.text(0.5, 0, '', transform=ax.transAxes, ha='center', fontsize=15,
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))


def update(i):
    im.set_array(gridHistory[i])
    time_text.set_text(f'J = {J}, kT = {kT:.2f}, M = {magnetizations[i]:.3f}, E = {energies[i]:.3f}')
    return im, time_text


ani = FuncAnimation(fig, update, frames=len(gridHistory),
                    interval=200, blit=True, repeat=False)  # Change speed by varying "interval"

# -------------------------------------------------------------------------------
# Spin-spin correlations of the final grid
# -------------------------------------------------------------------------------
fig = plt.figure(figsize=(10, 7))
plt.subplots_adjust(wspace=0.4)

ax1 = fig.add_subplot(121, projection="3d")
ax2 = fig.add_subplot(122)

# Connected correlation function, with zero displacement at the center:
corrMatrix = CorrelationFun(finalGrid)

# Displacements (dx, dy) matching the entries of corrMatrix:
d = np.arange(N) - N // 2
X, Y = np.meshgrid(d, d)

ax1.plot_surface(X, Y, corrMatrix, cmap='coolwarm')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('Corr')
ax1.set_title('Correlation', fontweight='bold')

# Average over directions, giving the correlation as a function of distance:
R = RadialAverage(corrMatrix, N)
ax2.plot(R)
ax2.set_xlabel('Distance')
ax2.set_ylabel('Correlation')
ax2.set_title('Radial average', fontweight='bold')

# Add cursor to click over line plot and get coordinate info
cursor = mplcursors.cursor(ax2.get_lines(), hover=False)

# -------------------------------------------------------------------------------
# Energy and magnetization over time
# -------------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
plt.subplots_adjust(wspace=0.3)

ax_en, ax_en_hist = axes[0]
ax_mag, ax_mag_hist = axes[1]

custom_blue = (0, 0.4470, 0.7410)

# Plot energy
ax_en.plot(times, energies, color=custom_blue)
ax_en.set_xlabel(f'time ({timeUnit})')
ax_en.set_ylabel('energy per spin')

# Plot energy histogram
ax_en_hist.hist(energies, bins='sturges', color=custom_blue, alpha=0.7, edgecolor='black')
ax_en_hist.set_xlabel('energy per spin')
ax_en_hist.set_ylabel('frequency')

# Plot magnetization
ax_mag.plot(times, magnetizations, color=custom_blue)
ax_mag.set_xlabel(f'time ({timeUnit})')
ax_mag.set_ylabel('magnetization per spin')

# Plot magnetization histogram
ax_mag_hist.hist(magnetizations, bins='sturges', color=custom_blue, alpha=0.7, edgecolor='black')
ax_mag_hist.set_xlabel('magnetization per spin')
ax_mag_hist.set_ylabel('frequency')

# Enable cursors upon clicking on a point:
cursor = mplcursors.cursor(ax_en.get_lines(), hover=False)
cursor = mplcursors.cursor(ax_mag.get_lines(), hover=False)

plt.show()
