## Generalization to more states or generalization of the step functions

import cfunctions
reload(cfunctions)
from cfunctions import *
import file_functions
reload(file_functions)
from file_functions import *

# Load already built default maps
basemaps = {}
for basename in bases:
    basemaps[basename] = load_file(maps_dir + basename + '.dict')
del basename

# Quantiles
quantiles = load_file(maps_dir + 'quantiles.dict')

# State space sizes
sizes = load_file(maps_dir + 'state_space_size.dict')

# maximum devo times (both for brent and develop)
devo_times = load_file(maps_dir + 'devo_times.dict')

# tf_gene_map from RegulonDB's E. coli TFSet
tf_gene_map = load_file(maps_dir + 'tf_gene_map.dict')


###########################
## STATE SPACE FUNCTIONS ##
###########################

## returns equally spaced states (steps yy) given min_y, max_y and
#  number of steps (states)
def get_equalprobable_states(
        n=n_states, min=min_state, max=max_state, dtype=Dtype):
    return linspace(min, max, n).astype(dtype)

## uniform jumps are always between [-1, 1], so they only depend on n_states
def get_jumps(n=n_states):
    step = 2./n
    return arange(-1 + step, 1, step)

## full enumeration of phenotype (vector) space: size = base.size^bits
#  its the full Cartesian Product of its base with itself
#  (permutations with replacement)
#  returns GENERATOR
def get_vector_space(base, bits=N):
 # from itertools: equivalent to nested for-loops in a generator expression
    return iproduct(base, bits)

## create and save maps: vector_space <-> decimal_base
def build_base_mapping(base, bits=N, extension=''):

    # list because we iterate twice
    vector_space = list(get_vector_space(base, bits))

    # build maps
    decimal_to_vector = dict(((decimal, array(vector))
                              for decimal, vector in enumerate(vector_space)))
    vector_to_decimal = dict(((tuple(vector), decimal)
                              for decimal, vector in enumerate(vector_space)))

    # save maps
    basename = get_basename(base)
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
    save_file(filename, basemap)

    return basemap

## Load or build state (non-default) maps
def load_base(basename, base, bits=N, extension=''):
    filename = maps_dir + basename + extension + '.dict'
    try:
        file = open(filename, 'rb')
    except IOError:
        basemaps[basename] = build_base_mapping(base, bits, extension)
    else:
        basemaps[basename] = cPickle.load(file)
        file.close()

    if (bits not in basemaps[basename]['decimal_to_vector'].keys() or
        bits not in basemaps[basename]['vector_to_decimal'].keys()):
        basemaps[basename] = build_base_mapping(base, bits)

##
def get_base(
        k=n_states, bits=N, min=min_state, steps='uniform', max=max_state):
    base = get_equalprobable_states(k, min, max)
    basename = get_basename(base)
    jumps = get_jumps(k) if steps == 'uniform' else quantiles[basename]
    ## ATTENTION: bias corrected
    extension = '_big' if bits > N else ''
    if k <= k_dict(bits) and bits <= N_dict:
        if basename not in basemaps:
            load_base(basename, base, bits, extension)
        elif (bits not in basemaps[basename]['decimal_to_vector'].keys() or
              bits not in basemaps[basename]['vector_to_decimal'].keys()):
            load_base(basename, base, bits, extension)
            #basemaps[basename] = build_base_mapping(base, bits, extension)
    return base, basename, jumps


'''
basemaps[basename]['vector_to_decimal'][10].items()[0], basemaps['gray']['vector_to_decimal'][10].items()[0]
basemaps[basename]['decimal_to_vector'][10].items()[0], basemaps['gray']['decimal_to_vector'][10].items()[0]
'''


################################
## STEP and SIGMOID FUNCTIONS ##
################################

# returns the step function conditional on the number of states (steps) and
# implementation (C or python)
def get_step_function(n=n_states, step_function_name='cstep'):
    return globals()[step_function_name + '%s' %(n-1)]

# input is 0-dim, jump is 1-dim, state is 2-dim
def step1(input, jump, state):
    if input < jump[0]:
        return state[0]
    else:
        return state[1]

# input is 0-dim, jump is 2-dim, state is 3-dim
def step2(input, jump, state):
    if input < jump[0]:
        return state[0]
    elif input < jump[1]:
        return state[1]
    else:
        return state[2]

# Mike's sigmoid function [0,1]
def sigmoid(x):
    return 1 / ( 1 + exp(-a*x) )

# Siegal's sigmoid function [-1,1]
def sigmoid2(x):
    return 2 / ( 1 + exp(-a*x) ) - 1

# too slow!!!
def sign0(x):
    v = sign(x)
    put(v, where(v == 0), -1.)
    return v


###############################################
## GRAY, BINARY and other bases CODE MAPPING ##
###############################################

def get_vector(decimal, basename, base, bits=N):
    # TODO: consider using binary_repr(decimal, width = bits) for base = (0,1)
    if base.size <= k_dict(bits) and bits <= N_dict:
        return basemaps[basename]['decimal_to_vector'][bits][decimal].copy()
    else:
        return decimal_to_vector(decimal, base, bits)

def get_decimal(vector, basename, base, bits=N):
    if base.size <= k_dict(bits) and bits <= N_dict:
        return basemaps[basename]['vector_to_decimal'][vector.size][
            tuple(vector)]
    else:
        return vector_to_decimal(vector, base)

## ATTENTION: MAX BASE SIZE = 10 !!!
def vector_to_decimal(vector, base, epsilon=1e-15):

    index = list(base).index # index funtion

    # I'm suffering from a floating-point precision problem
    # Defining tolerance as a work-around, probably very slow
    def index(x): return nonzero(abs(base - x) < epsilon)[0][0]

    decimal_str = "".join((str(index(x)) for x in vector))

    # Note: plain integer literals that are above the largest representable
    # plain integer (e.g., 2147483647 when using 32-bit arithmetic) are
    # accepted as if they were long integers instead 'L'.
    # But .__str__() works the same, i.e., it strips the 'L'.

    # Will use decimal base and not the one corresponding to the number of
    # states. It doesn't matter, as long as it is in 1:1 correspondence
    return int(decimal_str)#, base.size)

def decimal_str_to_vector(decimal_str, base):
    return array([base[int(x)] for x in decimal_str])

def decimal_to_vector(decimal, base, bits = N):
    ## int(decimal_str) strips the leading '0's of decimal_str,
    #  so we have to make sure to put them back to get the same number of bits
    decimal_str = decimal.__str__().rjust(bits, '0')
    return decimal_str_to_vector(decimal_str, base)

def binary_to_decimal(phenotype):
    boolean = sigmoid(phenotype).astype(int)  # fix(sigmoid(phenotype))
    #binary  = ''.join([str(x) for x in boolean])
    binary  = ''.join(map(str,boolean))
    #decimal = string.atoi(binary, 2)
    decimal = int(binary, 2)
    return decimal

def decimal_to_binary(d, size=N):
    if d == 0:
        return -ones(size, dtype=int)

    string = np.binary_repr(d, width=size)
    binary = array(tuple(string)).astype(int)

    # map from Masel's and Palmer's {0, 1} back to Wagner's and Siegal's {-1,1}
    binary[where(binary == 0)] = -1
    return binary

'''
decimal_to_binaryN10 = dict(((i,decimal_to_binary(i)) for i in range(pow(2,N))))
decimal_to_binary_map = {10:decimal_to_binaryN10}
file = open(maps_dir + 'decimal_to_binary.dict', 'wb')
cPickle.dump(decimal_to_binary_map, file)

binary_to_decimalN10 = dict(((tuple(decimal_to_binary(i)),i) for i in range(pow(2,N))))
binary_to_decimal_map = {10:binary_to_decimalN10}
file = open(maps_dir + 'binary_to_decimal.dict', 'wb')
cPickle.dump(binary_to_decimal_map, file)
'''

def gray_to_decimal(phenotype):
    size = len(phenotype) # N
    weights = array([pow(2,i) - 1 for i in range(1,size+1)][::-1])

    # map from Wagner's and Siegal's {-1,1} to Masel's and Palmer's {0, 1}
    boolean = sigmoid(phenotype).astype(int) #;print boolean

    g = zeros(size, dtype=int)
    indices = nonzero(boolean)[0]
    g[indices[::2]] = boolean[indices[::2]]
    g[indices[1::2]] = -boolean[indices[1::2]]
    #print g, weights

    return dot(g, weights)

def decimal_to_gray(d, size=N):
    if d == 0:
        return -ones(size, dtype=int)

    weights = array([pow(2,i) - 1 for i in range(1,size+1)][::-1])
    g = zeros(size, dtype=int)

    while True:
        ln2 = int(log2(d))
        w = pow(2, 1+ln2) - 1
        g[where(weights == w)] = 1
        d = w - d
        if not d: break

    # map from Masel's and Palmer's {0, 1} back to Wagner's and Siegal's {-1,1}
    g[where(g == 0)] = -1

    return g

'''
gray_to_decimalN4 = dict(((tuple(decimal_to_gray(i)),i) for i in range(pow(2,4))))
gray_to_decimal_map = {10:gray_to_decimalN10, 4:gray_to_decimalN4}
gray_to_decimal_map[9] = gray_to_decimalN9
file = open(maps_dir + 'decimal_to_gray.dict', 'wb')
cPickle.dump(decimal_to_gray_map, file)
binary_to_grayN10 = dict(((binary_to_decimal_map[N][tuple(decimal_to_gray_map[N][i])],i) for i in range(pow(2,4))))
'''

def gray(i):
    """
    This function returns the i'th Gray Code.
    It is recursive and operates in O(log n) time.
    """
    if i == 0: return 0
    if i == 1: return 1
 
    ln2 = int(log2(i))
    # the grey code of index i is the same as the gray code of an index an 
    # equal distance on the other side of ln2-0.5, but with bit ln2 set
    pivot = 2**(ln2) - 0.5 # TODO: double everything so that we use no floats
    delta = i - pivot
    mirror = int(pivot - delta)
    x = gray(mirror)    # get the grey code of the 'mirror' value
    x = x + 2**(ln2)    # set the high bit
    return x

## Return log base 2 of x.
'''def log2(x):
    return math.log(x) / math.log(2)
'''


############################
## RegulonDB's E. coli TRN #
############################

def parse_list(list, parser):
    return map(parser, list)

def parse_network(network, parser):
    network[:,0] = parse_list(network[:,0], parser)
    network[:,1] = parse_list(network[:,1], parser)
    return network

def tf_to_gene(tf):
    try:
        # by convention I will only use the first gene
        return tf_gene_map['tf_to_gene'][tf][0]
        #return name[0].lower() + name[1:]
    except KeyError:
        # tf is already gene
        return tf

def gene_to_tf(gene):
    try:
        return tf_gene_map['gene_to_tf'][gene]
    except KeyError:
        return gene
    #return name[0].upper() + name[1:]

def build_tf_gene_map():

    # load file
    TFs = read_network('TFSet')

    # build map
    tf_gene_map = {'tf_to_gene':{}, 'gene_to_tf':{}}
    for tf, gene in TFs:
        genes = gene.split(',')
        # some TFs correspond to more than one gene
        tf_gene_map['tf_to_gene'][tf] = genes[:-1]
        for gene in genes:
            if gene:
                tf_gene_map['gene_to_tf'][gene] = tf

    # save map
    save_file(maps_dir + 'tf_gene_map.dict', tf_gene_map)

    return tf_gene_map
