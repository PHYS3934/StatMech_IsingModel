# Import libraries:
import numpy as np

def myNeighbors(s, N):
    """
    Take a list of linear indices s and return the linear indices of the
    neighbors of s on an N by N grid with periodic boundary conditions.
    """
    s = np.array(s) - 1  # index by zero
    adj = np.zeros((len(s), 4), dtype=int)
    
    # s = r*N + c
    r = s // N
    c = s % N
    
    # Calculate neighbor indices with periodic boundary conditions
    adj[:, 0] = np.mod(r + 1, N) * N + c     # down
    adj[:, 1] = np.mod(r - 1, N) * N + c     # up
    adj[:, 2] = r * N + np.mod(c + 1, N)     # right
    adj[:, 3] = r * N + np.mod(c - 1, N)     # left
    
    adj = adj + 1  # index by one again
    
    return adj
