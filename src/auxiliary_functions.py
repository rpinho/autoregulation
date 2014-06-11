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

# fixed for all elements of array x
def fixed_jitter(x, jitter=0.2):
    return x + jitter*rand() - .1

# random for each element of array x
def random_jitter(x, jitter=0.2):
    return x + jitter*rand(len(x)) - .1

# from https://stackoverflow.com/questions/8671808/matplotlib-preventing-overlaying-datapoints
def rand_jitter(arr):
    stdev = arr.max()/100.
    return arr + randn(len(arr)) * stdev

def jitter(x, y, s=20, c='b', marker='o', cmap=None, norm=None, vmin=None, vmax=None, alpha=None, linewidths=None, verts=None, hold=None, **kwargs):
    return scatter(rand_jitter(x), rand_jitter(y), s=20, c='b', marker='o', cmap=None, norm=None, vmin=None, vmax=None, alpha=None, linewidths=None, verts=None, hold=None, **kwargs)
