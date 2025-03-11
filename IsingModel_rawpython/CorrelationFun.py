# Import libraries:
import numpy as np
from scipy import signal

def CorrelationFun(A, doNorm=True):
    """
    Computes the correlation function of A
    INPUTS: - A: a square (spin) grid.
            - doNorm: (boolean), whether to subtract the lattice-average M squared
    Output is centered at size(A)/2 because of periodic boundary conditions
    """
    # Calculate correlation using convolution
    # Extend A by tiling a 180-degree rotated version to handle periodicity
    A_rot = np.rot90(A, 2)
    A_extended = np.tile(A_rot, (2, 2))
    
    # Use scipy's signal.convolve2d for the convolution
    cor = signal.convolve2d(A, A_extended, mode='same') / A.size
    
    if doNorm:
        cor = cor - np.mean(A)**2
    
    # Center the correlation function
    cor = np.roll(cor, (A.shape[0]//2, A.shape[1]//2), axis=(0, 1))
    
    return cor