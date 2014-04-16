"""
What follows is a simplified version of some of the main functions in my model.
If you're just looking for the nifty algorithm, scroll down to brent()
However, I didn't invent brent(), I just implemented it.
My original insight was to translate my problem to the form where brent() could be used.
To do that I had to restrict some parameters and check that results didn't change,
and to develop a mapping function.
The code is simple, but allowed for 5x speed-up,
and the complete enumeration of the state space (from just 7% for previous authors).
It also allowed for longer simulations, that revealed new phases of evolution.

At the core of my model is a boolean network with threshold update rules.
In essence, its a matrix multiplication between a NxN matrix and a 1xN vector,
followed by a step or sigmoid function that normalizes the output vector to [-1, 1].
This process is iterated until a steady state equilibrium (attractor) is found.

If the state vectors are binary {-1, 1} or {0, 1},
the state space of the model is finite and the dynamics deterministic.
The system will eventually reach an attractor given an initial vector.
The attractor may be a fixed point or a limit cycle.

Previous authors have resorted to a convergence rule,
analogous to calculating the mean variance of the 10 last iterations (vectors),
and requiring it to be less than a small residual error (1e-4).
This approach is slow, and only identifies fixed points, leaving out limit cycles
(that turn out to constitute the majority of the sate space).
It requires keeping a list of 10 or more vectors, calculating mean variance
and the tuning of multiple parameters.

I noticed that the update rule is a function that maps a finite set to itself.
Given any initial value, the sequence of iterated function values
must eventually use the same value twice.
Once this happens, the sequence must continue by repeating a cycle of values.
This is a cycle detection problem, and I decided to use Brent's algorithm.

But first, I had to map the binary vectors of size N to either
binary/gray code or decimal representation (integers)
"""

# dependencies and parameters for this sample code
from numpy import *
N = 10
def_base     = array([-1.,  1.])
def_basename = array2string(def_base, precision = 2)
def_gpmap    = sign


# create and save maps: vector_space <-> decimal_base
def build_base_mapping(base = def_base, bits = N, extension = ''):

    # full enumeration of vector space: size = base.size^bits
    # its the full Cartesian Product of its base with itself (permutations with replacement)
    # from itertools: equivalent to nested for-loops in a generator expression
    vector_space = list(itertools.product(base, repeat = bits)) # list because we iterate twice

    # build maps
    decimal_to_vector = dict(((decimal, array(vector)) for decimal, vector in enumerate(vector_space)))
    vector_to_decimal = dict(((tuple(vector), decimal) for decimal, vector in enumerate(vector_space)))

    # save maps
    basename = array2string(base, precision = 2)
    filename = maps_dir + basename + extension + '.dict'
    try:
        file = open(filename)
    # new file
    except IOError:
        basemap = {'decimal_to_vector': {}, 'vector_to_decimal': {}}
    else:
        basemap = cPickle.load(file)
        file.close()
    basemap['decimal_to_vector'].update({bits: decimal_to_vector})
    basemap['vector_to_decimal'].update({bits: vector_to_decimal})
    save_file(basemap, filename) # cPickle.dump()

    return basemap


# You can then do the mapping with:
# after loading basemaps with load_base(basename) - not shown
def get_vector(decimal, basename, bits = N):
    return basemaps[basename]['decimal_to_vector'][bits][decimal].copy()

def get_decimal(vector, basename, bits = N):
    return basemaps[basename]['vector_to_decimal'][vector.size][tuple(vector)]


# Now that we handle mapping, we can code the basic model (one iteration).
# Notation: genotype = NxN matrix, phe_decimal = decimal representation of 1xN vector
#           gpmap = sign, step or sigmoid function (filter/normalization)
def devo(genotype, phe_decimal, basename = def_basename, gpmap = def_gpmap,
         bits = N, *args):

    phenotype = get_vector(phe_decimal, basename, bits)
    # the development model: matrix multiplication (additive/linear)
    phenotype = dot(genotype, phenotype)
    # normalization/filter/step or sigmoid function: introduces nonlinearity
    phenotype = gpmap(phenotype)
    return get_decimal(phenotype, basename, bits)


# cycle detection algorithm: Brent's (the "tortoise and the hare")
# Notation: f = devo, x0 = initial state, max_power = maximum cycle size
def brent(f, genotype, x0, basename = def_basename, gpmap = def_gpmap,
          bits = N, max_power = inf):

    # lam is the period (length of the cycle)
    power = lam = 1
    # f(x0) is the element/node next to x0.
    tortoise, hare = x0, f(genotype, x0, basename, gpmap, bits)

    # Main phase: search successive powers of two
    while tortoise != hare:
        if power == lam:   # time to start a new power of two?
            tortoise = hare
            power *= 2 
            lam = 0
        hare = f(genotype, hare, basename, gpmap, bits)
        lam += 1

    # Find the position of the first repetition of length lambda (where the cycle starts: mu)
    mu = 0
    tortoise = hare = x0

    # lets advance hare one period ahead of the tortoise, which is at x0
    # we have to do this step by step (our model)
    for i in range(lam):
        hare = f(genotype, hare, basename, gpmap, bits)

    while tortoise != hare:
        # now that hare and tortoise are one period away from each other,
        # lets advance them together, and find mu where the tortoise stops
        tortoise = f(genotype, tortoise, basename, gpmap, bits)
        hare = f(genotype, hare, basename, gpmap, bits)
        mu += 1

    return lam, mu, tortoise # period, length, phe_decimal


"""
Brent increased the speed of the simulation by 5x, and allows to identify all
attractors, including cycles, which are the majority (93%).

When matrices (genotype) are binary and N is even,
and when vectors (phenotype) are binary {-1, 1},
dot(genotype, phenotype) in devo() can return 0.
Python's native function returns sign(0) = 0,
but our mapping doesn't allow for 0, only {-1, 1}.
To handle this exception I coded the following in C, using scipy.weave:
it compiles once, and runs faster than python code.
"""

from scipy import weave

# 0 -> 1
def csign(input):
    size   = input.size
    output = input.copy()
    code = """
           for (int i=0; i<size; ++i)
               output(i) = output(i) < 0 ? -1 : 1;
           return_val = py_output;
           """
    return weave.inline(code, ['output','size'], type_converters = weave.converters.blitz)
