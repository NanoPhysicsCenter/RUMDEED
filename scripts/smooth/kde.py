import matplotlib.pyplot as plt
import numpy as np
import numba
import math

@numba.jit
def kernel(a, b, c, x):
    # a: height
    # b: center
    # c: width
    # x: x-axis grid
    return a*np.exp(-(x-b)**2/(2*c**2))

@numba.jit
def kde_kt(x, y):

    FWHM = 0.25 # Full Width at half maximum
    c = FWHM/(2.0*np.sqrt(2.0*np.log(2.0)))

    N = 1
    grid = np.linspace(np.min(x), np.max(x), x.size*N)
    kde = np.zeros(len(grid))
    M = len(x)
    for i in range(M):
        a = y[i]
        b = x[i]
        kde += kernel(a, b, c, grid)

    return kde, grid

filename = '../../data/emitted.dt'
x, y = np.loadtxt(filename, usecols=(0, 2), unpack=True)

#x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
#y = np.array([0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0])

kde, grid = kde_kt(x, y)
plt.plot(grid, kde, x, y)
plt.show()
