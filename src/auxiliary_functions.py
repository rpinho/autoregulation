import header
reload(header)
from header import *

########################
## Auxiliaty functions #
########################

def nondiag(x):
    return triu(x, k=1) + tril(x, k=-1)

def flatnondiag(x):
    return x.flatten()[delete(xrange(x.size), xrange(0, x.size, len(x)+1))]

def difference_diag_nondiag(x, function=array, bits=N):
    if not any(x):
        return NaN
    if x.ndim == 1:
        x = x.reshape((bits, bits))
    return function(diag(x)) - function(nondiag(x))

def get_int_datapoints(start, stop, num, filter_i=None):
    x = unique(logspace(start, stop, num).astype(int))
    if any(filter_i):
        return delete(x, filter_i)
    return x

def print_min_max(x1, x2):
    return '%d_%d' %(x1[0], x2[-1])

def gradient_norm(x):
    Gx, Gy = gradient(x)
    return (Gx**2+Gy**2)**.5

def find_nearest_i(x, value):
    return (abs(x-value)).argmin()

def find_nearest(x, value):
    return x[find_nearest_i(x, value)]

def add_random_noise(x, noise=0.2):
    return x + noise*rand() - .1
