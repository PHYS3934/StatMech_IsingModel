# Import libraries:
import numpy as np

def CorrelationFun(A, doNorm=True):
    """
    Computes the spin-spin correlation function of a square spin grid, A, with periodic boundaries:
        C(dx, dy) = mean over all sites (x, y) of A[x, y] * A[x + dx, y + dy]
    INPUTS: - A: a square (spin) grid.
            - doNorm: (boolean), whether to subtract the squared lattice-average magnetization,
                      giving the connected correlation function
    The output is centered: zero displacement (dx = dy = 0) is at index [N//2, N//2].
    """
    # The straightforward way loops over every displacement, shifting the grid each time:
    #     cor[dx, dy] = np.mean(A * np.roll(A, (-dx, -dy), axis=(0, 1)))
    # The Fourier transform gives the same result for all displacements at once, because the
    # autocorrelation is the inverse transform of the power spectrum, |FFT(A)|^2
    # (the Wiener-Khinchin theorem; the FFT assumes periodic boundaries, as our lattice has):
    cor = np.real(np.fft.ifft2(np.abs(np.fft.fft2(A))**2)) / A.size

    if doNorm:
        cor = cor - np.mean(A)**2

    # Move zero displacement from the corner [0, 0] to the center [N//2, N//2]:
    cor = np.fft.fftshift(cor)

    return cor
