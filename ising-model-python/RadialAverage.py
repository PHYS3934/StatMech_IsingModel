import numpy as np

def RadialAverage(cor, N):
    """
    Compute the radial average of the NxN correlation function, cor
    (average out the angular dependence from a 2D connected correlation function)
    """
    L = N // 2 + (N % 2)  # ceiling of N/2
    c = np.arange(1, N+1) - L
    X, Y = np.meshgrid(c, c)
    rho = np.sqrt(X**2 + Y**2)  # equivalent to cart2pol in MATLAB
    
    rbins = np.arange(-0.5, L+0.5)  # only go to N/2, because of periodic boundaries
    R = np.zeros(L)
    
    for j in range(L):
        r = rbins[j]
        ring = (r <= rho) & (rho < r+1)
        if np.sum(ring) > 0:  # Avoid division by zero
            R[j] = np.sum(cor[ring]) / np.sum(ring)
    
    # Replace NaN values with 0
    R[np.isnan(R)] = 0
    
    return R
