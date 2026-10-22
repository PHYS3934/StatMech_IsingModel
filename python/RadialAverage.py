# Import libraries:
import numpy as np

def RadialAverage(cor, N):
    """
    Compute the radial average of the NxN correlation function, cor
    (average out the angular dependence from a 2D connected correlation function).
    Zero displacement must be at the center, cor[N//2, N//2], as returned by CorrelationFun.
    R[r] is the average of cor over all displacements at distance r (to the nearest integer).
    """
    # Distance of each displacement (dx, dy) from zero displacement at the center:
    d = np.arange(N) - N // 2
    dx, dy = np.meshgrid(d, d)
    distance = np.sqrt(dx**2 + dy**2)

    # Only go out to about N/2: beyond that, displacements wrap around the periodic boundaries
    numR = (N + 1) // 2
    R = np.zeros(numR)

    for r in range(numR):
        # Average over the ring of displacements with r - 0.5 <= distance < r + 0.5:
        ring = (r - 0.5 <= distance) & (distance < r + 0.5)
        R[r] = np.mean(cor[ring])

    return R
