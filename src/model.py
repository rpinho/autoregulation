import mapping
reload(mapping)
from mapping import *

# Defines default base map (ATTENTION: may change the one defined in header.py).
# Note: '[-1.  1.]' faster than 'gray' or binary'
def_base = get_equalprobable_states(n_states, min_state, max_state)
def_basename = array2string(def_base, precision=2)
if def_basename not in basemaps:
    load_base(def_basename, def_base)

# Defines default normalization function
# (ATTENTION: may change the one defined in header.py)
#step = get_step_function(n_states, 'cstep')
#def def_gpmap(x): return step(x, def_jumps, def_base)
# Note: cstep 1 slower than sign, but faster than cstep.
# (cstep2 is also faster than cstep)

# compile as much code as possible to speedup simulation
#psyco.full()
#psyco.log()
#psyco.profile()

# Set seed based on current time
seed()


##
def get_gpmap(
        base=def_base, basename=def_basename, density=c, bits=N, binary=False,
        jumps=array([ 0.]), alt_gpmap=csign, heaviside=cHeaviside):

    if base.size == 2:
        # if density < 1 or small N or binary matrices you might need to
        # correct for sign(0) = 0.
        # For [-1.  1.] csign is faster than cstep1 and still corrects
        # (sign is the fastest!)
        if basename == def_basename:
            # if matrix is binary:
            if binary:
                # if there are 0s in the matrix
                if density < 1:
                    return alt_gpmap
                # if N is even
                if not bits%2:
                    return alt_gpmap
            return def_gpmap

        # for [ 0.  1.] cHeaviside is faster than cstep1
        else:
            return heaviside

    else:
        if base.size <= k_dict(bits):
            step = get_step_function(base.size, 'cstep')
        else:
            step = cstep
        def gpmap(x):
            """cstep"""#%step.__name__
            return step(x, jumps, base)
        return gpmap

# return name of the development model used
def get_model_name(gpmap):
    if gpmap is sign:
        return '_wagners'
    elif gpmap is sigmoid2:
        return '_siegals%s' %(a)
    else:
        return '_mikes%s' %(a)

# return founder function
def get_founder_function(gpmap):
    if gpmap is sign:
        return wagners_founder
    else:
        return siegals_founder



#################################
## EUCLIDEAN DISTANCE FUNCTIONS #
#################################

## For real vectors - Siegal and Mike

# Euclidean vector distance function
## ATTENTION - squared!
def squared_vector_distance(v1, v2):
    d = v1 - v2
    return vdot(d,d) #square(d).sum() #sum(d*d) # in order of faster to slower

# Normalized between [0,1] - Siegal /4N, Mike /N
def normed_squared_vector_distance(v1, v2):
    return squared_vector_distance(v1, v2) / len(v2)

# a measure analogous to a variance E[(x - E(x))^2]
def avg_dist(phenotypes):
    avg_phenotype = mean(phenotypes, axis=0)
    d = 0.
    for phenotype in phenotypes:
        d += normed_squared_vector_distance(phenotype, avg_phenotype)
    return d/len(phenotypes)

# convergence test
def convergence(phenotypes, converge_threshold=converge):
    distance = mean(var(phenotypes, axis=0))
    # SLOW. custom slower alternative = avg_dist(phenotypes).
    # Var is the Bottleneck!
    return distance < converge_threshold

# Notice: genotypes are real number matrices, so can't use hamming distance
def genotype_distance(a, b):
    d = 0.
    for gene1, gene2 in zip(a, b):
        d += sqrt(squared_vector_dist(gene1, gene2))
    return d/len(a)

#
def distance_matrix(population, type=float64):
    d = array([], dtype=type)
    for i, ind1 in enumerate(population.individuals):
        for ind2 in population.individuals[i+1:]:
            d = append(d, genotype_distance(ind1.genotype, ind2.genotype))
    return d



###############################
## HAMMING DISTANCE FUNCTIONS #
###############################

## For binary vectors - Wagner

# for binary strings a and b the Hamming distance is equal to
# the number of ones in a XOR b. [Wikipedia]
def int_decimal_hamming_distance(decimal1, decimal2, bits=N):
    return binary_repr(bitwise_xor(decimal1, decimal2), width=bits).count('1')

# FASTER! int_decimal_hamming_distance == int_vector_hamming_distance
def int_vector_hamming_distance(v1, v2):
    return where(v1 != v2)[0].size #nonzero(v1 != v2)[0].size

# from [Wikipedia]: phenotypes are strings. STRINGS only, not ints!!
def int_str_hamming_distance(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

## Normalized between [0,1]

# decimals
def normed_decimal_hamming_distance(decimal1, decimal2, bits=N):
    return int_decimal_hamming_distance(decimal1, decimal2, bits) / float(bits)

# vectors
def normed_vector_hamming_distance(v1, v2):
    return int_vector_hamming_distance(v1, v2) / float(v1.size)



######################
## FITNESS FUNCTIONS #
######################

# Siegal and Mike
## ATTENTION - in truth, this is also gaussian because the distance is squared!
def exponent_fitness(
        a, b, distance_function=normed_squared_vector_distance,
        sel_strength = s):
    return exp(-distance_function(a,b) / sel_strength)

# Wagner
def gaussian_fitness(
        a, b, distance_function=normed_vector_hamming_distance,
        sel_strength=s):
    return exp(-distance_function(a,b)**2 / sel_strength)



######################
## COMPARE FUNCTIONS #
######################

#
def compare_vectors(v1, v2):
    if   list( v1) == list(v2):
        return 1, 0 # equal
    elif list(-v1) == list(v2):
        return 0, 1 # symmetrical
    else:
        return 0, 0 # different

# vector by vector: returns tuple (# equal vectors, # symmetrical vectors)
def compare_matrices(g1, g2):
    # control shape/size
    if shape(g1) != shape(g2): return 0, 0
    return array([
        compare_vectors(v1, v2) for v1, v2 in zip(g1, g2)]).sum(axis=0)

# entry by entry: returns fraction of equal
def compare_genotypes(g1, g2):
    # control shape/size
    if shape(g1) != shape(g2): return 0
    return array([1 for x1, x2 in zip(g1.flatten(), g2.flatten()) if
                  x1 == x2]).sum() / float(g1.size)



###########################################
## DENSITIES, DEGREES AND NUMBER OF ZEROS #
###########################################

def get_in_degrees(genotype):
    return array([nonzero(gene)[0].size for gene in genotype])

def get_out_degrees(genotype):
    return get_in_degrees(genotype.T)

def get_zeros(genotype):
    return len(genotype) - get_in_degrees(genotype)

def get_connectivity(genotype):
    return get_in_degrees(genotype) / float(len(genotype))

# return: density, number of zeros, degree
def get_density_degree_nzeros(density, n_zeros, degree, bits=N):

    if degree == 'N-1':
        degree = bits-1

    # there's a precision issue here:
    # using a workaround (1+1e-15) while i don't figure it out
    c = 1 + 1e-15

    if density is not None:
        n_zeros = rint(bits * (c - density)).astype(int)
        degree = rint(density * bits).astype(int)
        return density, n_zeros, degree

    if n_zeros is not None:
        density = c - n_zeros/float(bits)
        degree = bits - n_zeros
        return density, n_zeros, degree

    if degree  is not None:
        density = degree/float(bits)
        n_zeros = bits - degree
        return density, n_zeros, degree

    # else:
    return c, 0, bits

# get number non-zero elements in vector
def get_nonzeros(size, density, deterministic):
    return rint(density*size) if deterministic else binomial(size, density)

# get number of pos, neg and zero elements in vector
def get_pos_neg_zeros(size, nonzeros, p, deterministic, qp_scale=qp_scale):
    p = p if deterministic else normal(p, qp_scale)
    pos      = rint(p*nonzeros).astype(int)
    neg      = rint(nonzeros - pos).astype(int)
    n_zeros  = rint(size - nonzeros).astype(int)
    return pos, neg, n_zeros

# custom made binary vector
def generate_pos_neg_zeros_binary_vector(size, density, p, deterministic):
    nonzeros = get_nonzeros(size, density, deterministic)
    pos, neg, n_zeros = get_pos_neg_zeros(size, nonzeros, p, deterministic)
    vector = [1] * pos + [-1] * neg + [0] * n_zeros
    shuffle(vector)
    return vector

# ATTENTION: this only works for binary! and its slower
def count_ones(vec):
    return list(vec).count(1)

def get_activating_fraction_binary(vec):
    return count_ones(vec) / float(vec.size)

# general funtion
def get_number_of_positive(vec):
    return where(array(vec) > 0)[0].size

def get_number_of_negative(vec):
    return where(array(vec) < 0)[0].size

def get_activating_fraction(vec):
    p = get_number_of_positive(vec)
    q = get_number_of_negative(vec)
    if   p: return p / float(p + q)
    elif q: return 0
    else:   return NaN

def get_number_of_positive_diagonal(genotype):
    return get_number_of_positive(diag(genotype))

def get_activating_fraction_diagonal(genotype):
    return get_activating_fraction(diag(genotype))

# nondiag(x) defined in src.mapping.py
def get_number_of_positive_non_diagonal(genotype):
    return get_number_of_positive(nondiag(genotype))

def get_activating_fraction_non_diagonal(genotype):
    return get_activating_fraction(nondiag(genotype))

def get_in_number_of_positive(genotype):
    return array([get_number_of_positive(gene) for gene in genotype])

def get_out_number_of_positive(genotype):
    return get_in_number_of_positive(genotype.T)

def get_activating_fraction_length2_autoregulatory(genotype):
    x = genotype * genotype.T
    return get_activating_fraction(x[triu_indices_from(x,1)])



##############
## GENOTYPES #
##############

## randomize existing network with node degree conservation:
# consists in shuffling the edges between the nodes
# keeps the same number of nodes and edges as in the original graph.
# Moreover, each node keeps the same degree
def randomize_genotype(genotype):
    new = genotype.flatten()
    i = nonzero(new)[0]
    edges = new[i]
    shuffle(edges)
    new[i] = edges
    return new.reshape(shape(genotype))

# gene by gene - all genes CAN have the same connectivity = c if
# deterministic == True
def gene_by_gene_genotype(
        distribution, bits=N, density=c, deterministic=True, dtype=Dtype,
        **args):
    genotype = zeros((bits, bits), dtype=dtype)
    nonzeros = get_nonzeros(bits, density, deterministic)
    for gene in genotype:
        gene[permutation(bits)[:nonzeros]] = array(
            distribution(size=nonzeros), dtype=dtype)
    return genotype

# real-valued
# Standard Normal distribution (mean=0, stdev=1).
# Note: consider using normal(mean, std) and transform act_fraction -> mean
def generate_genotype(
        bits=N, density=c, deterministic=True, dtype=Dtype, **args):
    return gene_by_gene_genotype(
        standard_normal, bits, density, deterministic, dtype)

# the whole matrix at once - produces variation in each gene's connectivity
# (average = c) EVEN if deterministic == True
def generate_whole_genotype(
        bits=N, density=c, deterministic=False, dtype=Dtype, **args):
    size = bits*bits
    genotype = zeros(size, dtype=dtype)
    nonzeros = get_nonzeros(bits, density, deterministic)
    genotype[permutation(size)[:nonzeros]] = array(
        standard_normal(nonzeros), dtype=dtype)
    # torna-lo bi-dimensional NxN
    return genotype.reshape(bits,bits)

# integer-valued -1 or 1:
# all genes CAN have the same connectivity = c if deterministic == True
# activating_fraction(genotype) normal dist around p
def spin_genotype(
        bits=N, density=c, act_fraction=p, deterministic=True, dtype=Dtype,
        **args):

    if act_fraction == p:
        return sign(generate_genotype(bits, density, deterministic, dtype))
    # another solution would be the transformation:
    # 2*randint(0,2,100)-1 [McDonald 2008]

    distribution = rv_discrete(
        name='spin', values=((1.,-1.), (act_fraction, 1 - act_fraction)))

    return gene_by_gene_genotype(
        distribution.rvs, bits, density, deterministic, dtype)

# Elhanan - produces variation in each gene's connectivity (average = c)
## ATTENTION - ALLOWS FOR c TO DECREASE OVER EVOLUTIONARY TIME
# - effective rewiring
## ATTENTION: activating_fraction(genotype) = p ONLY IF density = 1,
# else normal dist around p
def spin_whole_genotype(
        bits=N, density=c, act_fraction=p, dtype=Dtype, **args):

    size = bits*bits
    genotype = -ones(size, dtype=dtype)
    genotype[permutation(size)[:rint(act_fraction*size)]] = 1
    genotype[permutation(size)[:rint((1-density) *size)]] = 0
    return genotype.reshape(bits,bits)

#
def change_diagonal(genotype, diagonal, bits):
    # flatten genotype to change it
    genotype = genotype.flatten()
    # slice gets diagonal
    genotype[0::bits+1] = diagonal
    return genotype.reshape(bits,bits)

#
def ldu_genotype(diagonal, triangle, bits=N):
    genotype = empty(bits**2)
    ind_all  = range(bits**2)
    ind_diag = range(0, bits**2, bits+1)
    ind_tri  = delete(ind_all, ind_diag)
    genotype[ind_diag] = diagonal
    genotype[ind_tri ] = triangle
    return genotype.reshape(bits, bits)

# doesnt work
'''def qp_normal_genotype(
        bits=N, density=c, p=p, q=p, deterministic=True, **args):
    diagonal = sign(normal(mut_biases_dict[p],    1,           bits*density))
    triangle = sign(normal(mut_biases_dict[q], bits, (bits**2-bits)*density))
    return ldu_genotype(diagonal, triangle, bits)
'''
#
def qp_genotype(
        bits=N, density=c, p=p, q=p, deterministic=True, **args):

    diagonal = generate_pos_neg_zeros_binary_vector(
        bits, density, p, deterministic)

    triangle = generate_pos_neg_zeros_binary_vector(
        bits**2-bits, density, q, deterministic)

    return ldu_genotype(diagonal, triangle, bits)

## ATTENTION: only binary! integer-valued -1 or 1
##            here act_fraction is for DIAGONAL ONLY!
## For the non-diagonal elements it's the default of the genotype_function
def diagonal_p_genotype(
        bits=N, density=c, p=p, act_fraction=p,
        binary=True, deterministic=True,
        nonzeros=None, dtype=Dtype, **args):

    # get genotype
    if binary:
        genotype = spin_genotype(bits, density, act_fraction)
    else:
        genotype = generate_genotype(bits, density)
        genotype = sign(genotype)

    diagonal = generate_pos_neg_zeros_binary_vector(
        bits, density, p, deterministic)

    return change_diagonal(genotype, diagonal, bits)

def diagonal_0_genotype(bits=N, density=c, binary=False, **args):
    # get genotype
    if binary:
        genotype = spin_genotype(bits, density)
    else:
        genotype = generate_genotype(bits, density)

    return change_diagonal(genotype, 0, bits)

def ldu_composition_genotype(lu_genotype, d_genotype):
    return nondiag(lu_genotype) + diag(diag(d_genotype))

# x is 1 or -1 -> 0 for label_clusters
def zerox_genotype(genotype, x):
    genotype[where(genotype == x)] = 0
    return genotype

# return labels, n
def label_clusters(genotype, x):
    return ndimage.label(zerox_genotype(genotype.copy(), x))

# return slices, labels
def find_objects(genotype, x):
    labels = label_clusters(genotype, x)[0]
    return ndimage.find_objects(labels), labels

# biological network: exponential in-degree and
# power-law out-degree distributions
# generated by networkx.directed_configuration_model()
def generate_graph_genotype(
        bits=N, degree=2, f_round=rint, loc=0,
        genotype_function=generate_genotype,
        p1=poly1d([-0.37373024,  0.26167554]),
        p2=poly1d([-0.04509132,  0.36909606]),
        d={4: 0.0415, 5: -0.053, 6: -0.0975, 7: 0.04271429, 8: 0.004125,
           9: -0.02555556}):

    # get degree correction as a function of bits:
    # empirical power-law + linear fit + delta
    if p1: loc = exp(p1(log(bits)))
    if p2 and bits < 7: loc += p2(loc)
    if bits in d:       loc += d[bits]

    # get exponential in_degree_sequence (array)
    z_in  = f_round(
        stats.expon.rvs(loc = loc, scale = degree, size = bits)).astype(int)
    # get power-law out_degree_sequence
    #z_out = nx.create_degree_sequence(
     #   size, nx.utils.powerlaw_sequence,exponent=2) # list
    z_out = f_round(
        stats.pareto.rvs(degree, loc = loc, size = bits)).astype(int) # array
    # generates a random directed pseudograph with the given degree sequences
    D = nx.directed_configuration_model(list(z_in), list(z_out))
    # to remove parallel edges (still keeps self-loops)
    D = nx.DiGraph(D)
    # return the graph adjacency matrix as a NumPy array
    genotype = array(nx.to_numpy_matrix(D))
    # multiply by random or binary weights
    genotype *= genotype_function(bits)
    return D, genotype

# generated by networkx.generators.directed.scale_free_graph
def generate_scale_free_genotype(
        bits=N, alpha=.45, gamma=.05, delta_in=0.2, delta_out=0,
        genotype_function=generate_genotype, correction=True,
        p002=poly1d([-0.05217295,  1.0168033]),
        p02=poly1d([-0.0523467,  0.9082144]),
        gammas=[ 0.97,  0.94,  0.93,  0.92,  0.9 ,  0.9 ,  0.89,  0.9 ,  0.87,
                 0.87,  0.89,  0.86,  0.85,  0.85,  0.87,  0.87,  0.86,  0.85],
        p01=array([[ 6.   ,   0.1 ,   0.01,   0.9 ,   0.2 ],
                   [ 7.   ,   0.1 ,   0.02,   0.9 ,   0.2 ],
                   [ 8.   ,   0.1 ,   0.03,   0.9 ,   0.2 ],
                   [ 9.   ,   0.1 ,   0.03,   0.7 ,   0.2 ],
                   [ 10.  ,   0.1 ,   0.05,   0.8 ,   0.2 ],
                   [ 11.  ,   0.1 ,   0.07,   1.  ,   0.2 ],
                   [ 12.  ,   0.1 ,   0.07,   0.9 ,   0.2 ],
                   [ 13.  ,   0.1 ,   0.05,   1.  ,   0.2 ],
                   [ 14.  ,   0.1 ,   0.07,   0.8 ,   0.2 ],
                   [ 15.  ,   0.1 ,   0.07,   0.7 ,   0.2 ],
                   [ 16.  ,   0.1 ,   0.08,   0.7 ,   0.2 ],
                   [ 17.  ,   0.1 ,   0.08,   0.8 ,   0.2 ],
                   [ 18.  ,   0.1 ,   0.09,   0.9 ,   0.2 ],
                   [ 19.  ,   0.1 ,   0.11,   1.  ,   0.2 ],
                   [ 20.  ,   0.1 ,   0.12,   1.  ,   0.2 ],
                   [ 21.  ,   .2  ,   0.02,   0.8 ,   0.2]])):

    # get graph parameters as a function of bits (empirical power-law)
    # and swap parameters (i screwed up)
    if correction:
        if p002:
            alpha = .02
            gamma = p002(log(bits))
        if gammas and bits < 6:
            gamma = gammas[bits-4]
        if p02 and bits > 21:
            alpha = .2
            gamma = p02(log(bits))
        delta_out = delta_in
        delta_in  = gamma
        gamma     = 1 - alpha - gamma
        if any(p01) and bits in p01.T[0]:
            alpha, gamma, delta_in, delta_out = p01[bits-6][1:]

    # generates a scale free directed pseudograph
    D = nx.generators.directed.scale_free_graph(
        bits, alpha, 1 - alpha - gamma, gamma, delta_in, delta_out)
    # to remove parallel edges (still keeps self-loops)
    D = nx.DiGraph(D)
    # return the graph adjacency matrix as a NumPy array
    genotype = array(nx.to_numpy_matrix(D))
    # multiply by random or binary weights
    genotype *= genotype_function(bits)
    return D, genotype

def test_period(period, stable=True, sel_period=1):
    return ((period == sel_period and stable) or
            (period != sel_period and not stable))

# returns matrix and initial state with desired period
# (fixpoint by default, i.e., stable - for at least one initial)
def generate_stable_genotype(
        bits=N, density=c, basename=def_basename, base=def_base,
        genotype_function=generate_genotype, p=p, q=p, act_fraction=p,
        binary=True, stable=True, sel_period=1, max_trials=F, dtype=Dtype):

    gpmap = get_gpmap(base, basename, density, bits, binary)

    #while True:
    for n in range(max_trials):
        initial = generate_decimal_phenotype(base, bits)
        genotype = genotype_function(
            bits, density, p=p, q=q, act_fraction=act_fraction, binary=binary,
            dtype=dtype)
        period, length, phe_decimal = brent(
            devo, genotype, initial, basename, gpmap, base, bits)
        if test_period(period, stable, sel_period):
            return initial, genotype, phe_decimal

    return 0, 0, 0

# full enumeration of (binary) genotype space
def generate_full_genotype_space(bits=4, n_zeros=0):
    # returns GENERATORS
    # n_zeros = number of zeros
    if not n_zeros:
        return (reshape(genotype, (bits, bits))
                for genotype in selections(def_base, bits*bits))
    else:
        # base is '[-1.  0.  1.]'
        base, basename = get_base(3, bits)[:2]
        # all possible collumn vectors in decimal form
        decimal_vector_map = basemaps[basename]['decimal_to_vector'][bits]
        decimals = [decimal
                    for decimal, vector in decimal_vector_map.iteritems()
                    if nonzero(vector)[0].size is bits - n_zeros]
        # all possible combinations of the decimal base collumns
        # in groups of bits
        genotypes = selections(decimals, bits)
        # get vectors
        return (array([get_vector(decimal, basename, base, bits).transpose()
                       for decimal in genotype]) for genotype in genotypes)

# symmetry constraints: can exchange collumns, order does not matter
def generate_all_symm_genotypes(bits=4, n_zeros=0):
    # returns GENERATORS

    # n_zeros = number of zeros
    if not n_zeros:
        # base is '[-1.  1.]'
        #base = def_base
        #basename = def_basename
        # all possible collumn vectors in decimal form
        #decimals = xrange(base.size**bits)
        return generate_all_symm_vector_genotypes(bits)

    else:
        # base is '[-1.  0.  1.]'
        base, basename = get_base(3, bits)[:2]
        # all possible collumn vectors in decimal form
        decimal_vector_map = basemaps[basename]['decimal_to_vector'][bits]
        decimals = [decimal
                    for decimal, vector in decimal_vector_map.iteritems()
                    if nonzero(vector)[0].size is bits - n_zeros]

    # all possible combinations of the decimal base collumns in groups of bits
    genotypes = selections(decimals, bits)

    # memory is OK
    if bits is 4:
        genotypes = list(genotypes)
        # order of collumns does not matter
        for genotype in genotypes:
            genotype.sort()
        # list objects are unhashable
        genotypes = unique(map(tuple, genotypes))

    # memory > 4GB
    else:
        symm_genotypes = []
        for genotype in genotypes:
            # order of collumns does not matter
            genotype.sort()
            if genotype not in symm_genotypes:
                symm_genotypes.append(genotype)
        genotypes = symm_genotypes

    # get vectors
    return (array([
        get_vector(decimal, basename, base, bits).transpose()
        for decimal in genotype]) for genotype in genotypes)

def generate_all_symm_decimal_genotypes(bits=4):
    return itertools.combinations_with_replacement(
        xrange(def_base.size**bits), bits)

def generate_all_symm_vector_genotypes(bits=4):
    decimal_genotypes = generate_all_symm_decimal_genotypes(bits)
    return (array([
        get_vector(decimal, def_basename, def_base, bits).T
        for decimal in genotype]) for genotype in decimal_genotypes)

# symmetry constraints: in anti-symmetric maps {-1,1}, dot(g, -p) = -dot(g, p)
# for N = 4 there's 8 symmetrical genotypes in a total of 3876, i.e.,
# we can neglect it
def filter_symm_genotypes(genotypes):
    size = shape(genotypes)[1] # bits
    b = []
    # all against all
    for i, g1 in enumerate(genotypes):
        for g2 in genotypes[i+1:]:
            b.append(compare_matrices(g1, g2))

    return b#genotypes


###############
## PHENOTYPES #
###############

## Vectors

# integer-valued (-1,1) or (0,1): A. Wagner and Siegal.
# Note: fraction of repressing (-1 or 0) given by parameter rep_fraction = 1-p
def spin_phenotype(bits=N, min_=-1., rep_fraction=1-p, dtype=Dtype):
    phenotype = ones(bits, dtype=dtype)
    phenotype[permutation(bits)[:rep_fraction*bits]] = min_
    # ATTENTION: has to be permutation because randint produces repeated numbers
    return phenotype

# multiple integer-value, i.e., more than 2 states.
# All states are equally-probable
def int_phenotype(base, bits=N, dtype=Dtype):
    return array([base[randint(base.size)] for i in range(bits)], dtype=dtype)

# real-valued [min, 1]
def real_phenotype(bits=N, min_=-1., **args):
    return array(uniform(min_, 1, bits))

# general caller function
def generate_vector_phenotype(
        base, bits=N, noise=0, rep_fraction=1-p, dtype=Dtype, **args):
    # discrete
    if not noise:#if any(base):
        # binary phenotypes: '[-1.  1.]' or '[ 0.  1.]'
        if base.size is 2:
            return spin_phenotype(bits, base[0], rep_fraction, dtype)
        # multiple states
        else:
            return int_phenotype(base, bits, dtype)
    # real-valued
    else:
        return real_phenotype(bits, base[0])

## Decimals

# decimal form of int_phenotype(). All integer states are equally-probable.
# ATTENTION: MAX SIZE IS BASE.SIZE = 10 !!!
def decimal_phenotype(base, bits=N):
    # phenotype vector ALREADY in the STATE BASE
    # (no need to call vector_to_decimal to convert from real to state base)
    vector = randint(base.size, size=bits)
    decimal_str = "".join( (str(x) for x in vector) )
    return int(decimal_str)

# general caller function
def generate_decimal_phenotype(base, bits=N, initials=None):
    if initials:
        return choice(initials)
    elif base.size <= k_dict(bits) and bits <= N_dict:
        return randint(base.size**bits)
    else:
        return decimal_phenotype(base, bits)

# full enumeration of phenotype space
# distinguishes between with (vectors) or without (decimals) noise
def generate_full_phenotype_space(noise, basename, base, bits=N):
    # returns GENERATORS.
    # ATTENTION: array(xrange) works fine but array((int...)) does not!
    if not noise:
        # decimals
        if base.size <= k_dict(bits) and bits <= N_dict:
            return xrange(base.size**bits)
        else:
            return (int("".join( (str(x) for x in vector) ))
                    for vector in get_vector_space(range(base.size), bits))
    # vectors
    else:
        return (array(vector) for vector in get_vector_space(base, bits))


## Symmetry constraints

# combinatorical symmetric, i.e., order does NOT matter
def get_symm(decimal, basename, base, bits):
    vector = get_vector(decimal, basename, base, bits)
    vector.sort()
    return get_decimal(vector, basename, base, bits)

# returns arrays
def generate_all_symm_phenotypes(noise, basename, base, bits=N):
    if (base.size <= k_dict(bits) and bits <= N_dict and
        basename not in basemaps):
        load_base(basename, base, bits)
    # list of phe_decimals because unique only works good with decimals
    # (unique is 1d)
    phenotypes = (get_symm(decimal, basename, base, bits)
                  for decimal
                  in generate_full_phenotype_space(0, basename, base, bits))
    # get unique: we have repeated decimals after get_symn()
    phenotypes = unique(phenotypes)
    # return either decimals or vectors
    if not noise:
        # decimals
        return phenotypes
    else:
        return array([
            get_vector(decimal, basename, base, bits)
            for decimal in phenotypes]) # vectors

def generate_all_symm_vector_phenotypes(base, bits=N):
    return itertools.combinations_with_replacement(base, bits)


#
def get_nearest_neighbours(
        ini_decimal, n=1, bits=N, basename=def_basename, base=def_base):
    initial = get_vector(ini_decimal, basename, base, bits)
    neighbours = []
    for i in range(initial.size):
        perturbation = initial.copy()
        perturbation[i] = -perturbation[i]
        perturb_decimal = get_decimal(perturbation, basename, base, bits)
        neighbours.append(perturb_decimal)
        if n-1:
            neighbours.extend(
                get_nearest_neighbours(perturb_decimal, n-1, bits, basename))
    return unique(neighbours) # list of ini_decimals


###################################
## GENOTYPIC AND SPACE GENERATORS #
###################################

# Note: if get_initials: have to untuple afterwards.
def get_initials_and_genotypes(
        full_enum, stable, samples, density=c, basename=def_basename,
        base=def_base, bits=N, genotype_function=generate_genotype,
        genotypic_space_function=generate_all_symm_genotypes,
        get_initials=False,
        phenotypic_space_function=generate_full_phenotype_space,
        noise=0, p=p, q=p, act_fraction=p,
        network=None, filters=[], knockouts=[],
        random=False, initials=None, binary=False,
        graph=False, scale_free=False,
        sel_period=1, max_trials=F, deterministic=True, shuffling=False):

    if not noise:
        phenotype_function = generate_decimal_phenotype
    else:
        phenotype_function = generate_vector_phenotype

    # all genotypes (binary): full enumeration of genotypic and
    # phenotypic spaces, with or without symmetry constraints
    if full_enum and bits <= enum_bits:
        n_zeros = int(bits * (1 - density)) # number of zeros
        genotypes = genotypic_space_function(bits, n_zeros)

        if not get_initials:
            # GENERATOR
            return genotypes

        else:
            # GENERATOR if generate_full_phenotype_space else array

            initials = phenotypic_space_function(noise, basename, base, bits)

            # if we don't cycle in order, it should increase the search speed of
            # unique(phenotypes)
            # but actually, NOT. And shuffle() is too slow, so use False
            if shuffling:
                genotypes = list(genotypes)
                initials  = list(initials)
                shuffle(genotypes)
                shuffle(initials)

            try:
                # equivalent to nested for-loops in a generator expression.
                # So even if you shuffle, one variable keeps constant while the
                # other loops
                # initials are fixed for every genotype.
                # NOTE: the discovery of unique(phenot) is faster if
                # you cycle first through all genotypes and only then
                # through initials
                return itertools.product(initials, genotypes) # GENERATOR

            except AttributeError:
                # if python < 2.6
                genotypes = list(genotypes)
                initials  = list(initials)
                return ((x,y) for x in initials for y in genotypes)
                # Note: random_product(initials, genotypes) repeats items,
                # so its basically random binary, not full_enum
                # (from http://docs.python.org/library/itertools.html#recipes)

    # graph with biological topology
    elif graph:
        if scale_free:
            graph_genotype = generate_scale_free_genotype
        else:
            graph_genotype = generate_graph_genotype
        if not get_initials:
            return (graph_genotype(bits, genotype_function=genotype_function)[1]
                    for i in xrange(int(samples)) )
        else:
            return ((phenotype_function(base, bits, initials=initials),
                     graph_genotype(
                         bits, genotype_function=genotype_function)[1])
                    for i in xrange(int(samples)))

    # RegulonDB's E. coli TRN
    elif network:
        if not get_initials:
            return (get_network(network, filters, knockouts, random)[0]
                    for i in xrange(int(samples)) )
        else:
            return ((phenotype_function(base, bits, initials=initials),
                     get_network(network, filters, knockouts, random)[0])
                    for i in xrange(int(samples)) )

    # conditional sample of (real) genotype space
    elif stable:
        if not get_initials:
            return (generate_stable_genotype(
                bits, density, basename, base, genotype_function, p, q,
                act_fraction, binary, stable, sel_period, max_trials)[1]
                    for i in xrange(int(samples)) )
            # [0] is initial, [1] genotype, [2] phenotype
        else:
            return (generate_stable_genotype(
                bits, density, basename, base, genotype_function, p, q,
                act_fraction, binary, stable, sel_period, max_trials)[:2]
                    for i in xrange(int(samples)) )
            # [0] is initial, [1] genotype, [2] phenotype
        # ATTENTION: is returning DECIMALS only!

    # random sample of (real) genotype space
    else:
        if not get_initials:
            return (genotype_function(
                bits, density, p=p, q=q, act_fraction=act_fraction,
                binary=binary, deterministic=deterministic)
                    for i in xrange(int(samples)) )
        else:
            return ((phenotype_function(base, bits, initials=initials),
                     genotype_function(
                         bits, density, p=p, q=q, act_fraction=act_fraction,
                         binary=binary, deterministic=deterministic))
                    for i in xrange(int(samples)) )


# General function
def get_individuals(
        full_enum, stable, samples, density=c,
        basename=def_basename, base=def_base, bits=N,
        genotype_function=generate_genotype,
        genotypic_space_function=generate_all_symm_genotypes,
        get_initials=True,
        phenotypic_space_function=generate_full_phenotype_space,
        noise=0, p=p, q=p, act_fraction=p,
        genotypes=None, network=None, filters=[], knockouts=[],
        random=False, initials=None, binary=False,
        graph=False, scale_free=False,
        sel_period=1, max_trials=F, deterministic=True, shuffling=False):

    # if genotypes are saved in a file
    if genotypes:
        individuals = cPickle.load(open(genotypes))
        get_initials = False
    else:
        individuals = get_initials_and_genotypes(
            full_enum, stable, samples, density, basename, base, bits,
            genotype_function, genotypic_space_function, get_initials,
            phenotypic_space_function,
            noise, p, q, act_fraction, network, filters, knockouts, random,
            initials, binary, graph, scale_free, sel_period, max_trials,
            deterministic, shuffling)

    # if individuals are only genotypes
    n_samples = samples
    if not get_initials and not initials:
        # GENERATOR if generate_full_phenotype_space else array
        initials = array(
            list(phenotypic_space_function(noise, basename, base, bits)))
    if any(initials):
        n_samples *= initials.size

    return individuals, initials, n_samples


#############
## MUTATION #
#############

def mutation(population):
# note: 1) if you want to keep old population call with a copy, e.g. pop[:]
    new_pop = []
    while population:
        new = mutate(population.pop())
        new_pop.append(new)
    new_pop.reverse() # reverse is optional, but easier to debug
    return new_pop

def get_mutating(genotype, n_mutations):
    # only mutate non-zero elements
    mutable = nonzero(genotype.flatten())[0]
    # random, no repetition
    mutating = permutation(mutable.size)[:n_mutations]
    return mutable[mutating]

# allows zero, single or multiple mutations.
# either Normal, change of sign or deletion
# mutation rate per non-zero entry = u/c*N^2
def mutate(
        genotype, n_mutations=None, rate=u,
        change_sign=False, deletion=False, bias=False):

    if not n_mutations:
        # ~ poisson(rate)
        n_mutations = binomial(c*genotype.size, rate/(c*genotype.size))

    # only mutate non-zero elements, randomly
    mutating = get_mutating(genotype, n_mutations)

    # house-of-cards assumption (large effects, independence)
    mutant = random_or_sign_or_deletion_mutation(
        genotype, mutating, change_sign, deletion, bias)

    # return it bi-dimensional NxN
    size = len(genotype)
    return mutant.reshape(size, size)

# NOTE: returns FLAT genotype
def random_or_sign_or_deletion_mutation(
        genotype, mutating, change_sign=False, deletion=False, bias=False):

    mutant = genotype.flatten() # makes a flat copy
    mutating = array(mutating)
    if deletion:
        mutant[mutating] = 0.
    elif change_sign:
        mutant[mutating]*= -1
    elif bias:
        mutant[mutating] = sign(normal(bias, 1, mutating.size))
    else:
        mutant[mutating] = array(
            standard_normal(mutating.size), dtype=mutant.dtype)
    return mutant


############
## MUTANTS #
############

# multiple mutant sample:
# each genotype differs from the PREVIOUS by one-point mutation
def get_multiple_mutants(genotype, samples, change_sign=False, deletion=False):
    mutants = []
    for i in range(int(samples)):
        genotype = mutate(genotype, 1, None, change_sign, deletion)
        mutants.append(genotype)
    return mutants

# single mutant sample:
# each genotype differs from the WILDTYPE by one-point mutation
def get_single_mutants(
        genotype, samples, change_sign=False, deletion=False):

    return (mutate(genotype, 1, None, change_sign, deletion)
            for i in range(int(samples)))

# NOTE: returns FLAT genotypes (generator)
def get_all_single_mutants(
        genotype, change_sign=False, deletion=False, *arg):

    wildtype = genotype.flatten()#reshape(genotype.size)
    mutable = nonzero(wildtype)[0]

    return (random_or_sign_or_deletion_mutation(
        wildtype.copy(), mutating, change_sign, deletion)
            for mutating in mutable)

#
def get_all_multiple_mutants(
        genotype, change_sign=False, deletion=False, n_multiple=2):

    wildtype = genotype.flatten()#reshape(genotype.size)
    mutable = nonzero(wildtype)[0]
    mutants = []

    for i in range(mutable.size):
        mutant = random_or_sign_or_deletion_mutation(
            wildtype.copy(), mutable[i], change_sign, deletion)
        mutants.append(mutant)
        if n_multiple - 1:
            mutants.extend(
                get_all_multiple_mutants(
                    mutant, change_sign, deletion, n_multiple - 1))
    return mutants

# General function for mutants
def get_mutants(
        initial, genotype, samples, all_mutants=False, change_sign=False,
        deletion=False, multiple=False):

    if all_mutants:
        # sampling of mutants
        if all_mutants > 1:
            generate_mutant_function = get_all_multiple_mutants
        else:
            generate_mutant_function = get_all_single_mutants

        # get mutants
        mutants = generate_mutant_function(
            genotype, change_sign, deletion, all_mutants)

        # get unique
        if all_mutants > 1:
            mutants = unique(map(tuple, mutants))

        # reshape
        mutants = reshape(list(mutants), (-1, len(genotype), len(genotype)))

    else:
        # sampling of mutants
        if multiple:
            generate_mutant_function = get_multiple_mutants
        else:
            generate_mutant_function = get_single_mutants

        # get mutants
        mutants = generate_mutant_function(
            genotype, samples, change_sign, deletion)

    return ((initial, mutant) for mutant in mutants)


##############
## EPISTASIS #
##############

# get single mutant genotype
def get_single_mutant(
        genotype, mutating, change_sign=False, deletion=False):
    return random_or_sign_or_deletion_mutation(
        genotype, mutating, change_sign, deletion).reshape(shape(genotype))

# get double mutant genotype
def get_double_mutant(mut1, mut2, mutating):
    size = len(mut1)
    mutant = mut1.flatten()
    mutant[mutating] = mut2.flatten()[mutating]
    return mutant.reshape(size, size)

# genotypes
def get_double_mutants(
        genotype, mutations, change_sign=False, deletion=False):
    mutants = []
    # first single mutant
    # gene_i is index of non-zero element of flattened genotype
    for i, gene_i in enumerate(mutations):
        mutant_i = get_single_mutant(genotype, gene_i, change_sign, deletion)
        # second single mutant
        for gene_j in mutations[i+1:]:
            mutant_j = get_single_mutant(
                genotype, gene_j, change_sign, deletion)
            # double mutant
            mutant_ij = get_double_mutant(mutant_i, mutant_j, gene_j)
            # store
            mutants.append((mutant_i, mutant_j, mutant_ij))
    return mutants

# ATTENTION: hamming distance is NOT normalized
def get_epistasis(hamming_i, hamming_j, hamming_ij, bits=N):
    bits = float(bits)
    return (1-hamming_ij/bits) - (1-hamming_i/bits)*(1-hamming_j/bits)


#############################
## INSERTIONS AND DELETIONS #
#############################

# the whole matrix at once
def insert_and_delete(genotype, ins_rate=uA, del_rate=uD):

    # ~ poisson(ins_rate)
    n_insertions = binomial(genotype.size, float(ins_rate)/genotype.size)
    # ~ poisson(del_rate)
    n_deletions  = binomial(genotype.size, float(del_rate)/genotype.size)

    # a one-dimensional copy makes it easier to handle
    network = genotype.reshape(genotype.size).copy()

    zeros = where(network == 0.)[0]
    edges = nonzero(network)[0]

    inserting = permutation(zeros)[:n_insertions]
    deleting  = permutation(edges)[:n_deletions]

    network[inserting] = standard_normal()
    network[deleting]  = 0.

    size = len(genotype)
    return network.reshape(size, size)

# this is wrong: creates bias always towards 0.5, like an elastic force of a
# spring. Because it samples in proportion to number of existing edges or zeros
'''def insert_and_delete(genotype, n_mutations = None, rate = uA):
    if not n_mutations: n_mutations = binomial(genotype.size, float(rate)/genotype.size) # ~ poisson(rate)
    # a one-dimensional copy makes it easier to handle
    network = genotype.reshape(genotype.size).copy()
    # mutation rate per entry = uA/N^2
    rewiring = permutation(genotype.size)[:n_mutations]
    # house-of-cards assumption (large effects, independence)
    for edge in rewiring:
        network[edge] = 0. if network[edge] else standard_normal()
    # return it bi-dimensional NxN
    size = len(genotype)
    return network.reshape(size, size)
'''
# gene by gene: this is still wrong. Gene.size = 10 is sufficiently small to create bias towards 0.5 again. Because edges = [] or edges = [] is suf. likely
'''insert_and_delete(genotype, ins_rate = uA, del_rate = uD):
    network = genotype.copy()
    for gene in network:
        r = uniform()
        # insert/add connection
        if   r < float(ins_rate) / gene.size:
            zeros = where(gene == 0.)[0]
            if zeros.any():
                inserting = zeros[randint(zeros.size)] #permutation[zeros][0]
                gene[inserting] = standard_normal()
        # delete connection
        elif r < float(ins_rate + del_rate) / gene.size:
            edges = nonzero(gene)[0]
            if edges.any():
                deleting = edges[randint(edges.size)] #permutation[edges][0]
                gene[deleting]  = 0.
    return network
'''


####################
## NOISE FUNCTIONS #
####################

def random_noise(phenotype, noise, noise_size, **args):
    # NOTE: np.random.normal(mu, sigma, size), sigma (scale) == std
    return phenotype + array(
        normal(0, noise, noise_size), dtype=phenotype.dtype)

def get_phenotype_mutating(bits, n_flips, noise_time=''):
    if noise_time == 'shmulevich':
        # random number of flips (from Shmulevich2002b)
        # NOTE: this is NOT uniform!!!
        #return randint(2, size=bits).nonzero()[0]# == 1
        # this one is uniform:
        n_flips = randint(1, bits+1)
    else:
        # deterministic number of flips
        pass
    return permutation(bits)[:n_flips]

# call with .copy() if you do not want to replace
# not garanteed to be always different:
# mean(int_vector_hamming_distance(phenotype, noise(phenotype))) = noise_size/2
# for def_base
def random_flip(
        phenotype, noise, noise_size=None, base=def_base, noise_time=''):

    # noise_size = noise
    mutating = get_phenotype_mutating(phenotype.size, noise, noise_time)

    phenotype[mutating] = base[randint(base.size, size=mutating.size)]

    return phenotype

# garanteed to be always different:
# int_vector_hamming_distance(phenotype, noise(phenotype)) = noise_size
def force_different_flip(
        phenotype, noise, noise_size=None, base=def_base, noise_time=''):

    # noise_size = noise
    mutating = get_phenotype_mutating(phenotype.size, noise, noise_time)

    for i in mutating:
        mutate_to = delete(base, where(phenotype[i] == base)[0])
        phenotype[i] = choice(mutate_to)

    return phenotype

# garanteed to be always different:
# int_vector_hamming_distance(phenotype, noise(phenotype)) = noise_size
def neighbour_flip(
        phenotype, noise, noise_size=None, base=def_base, noise_time=''):

    # noise_size = noise
    mutating = get_phenotype_mutating(phenotype.size, noise, noise_time)

    for i in mutating:
        i_wildtype = where(phenotype[i] == base)[0]
        # you can only move up (0 -> 1)
        if i_wildtype == 0:
            phenotype[i] = base[1]
        # you can only move down (-1 -> -2)
        elif i_wildtype == base.size - 1: # == -1
            phenotype[i] = base[-2]
        # up or down
        else:
            mutate_to = delete(base, where(phenotype[i] == base)[0])
            phenotype[i] = choice(mutate_to)

    return phenotype


################
## DEVELOPMENT #
################

# for populations
def development(population, initial, target, gpmap):
# note: 1) if you want to keep old population call with a copy, e.g. pop[:]
#       2) calculates fitness of stable individuals and returns it in a list of
#       tuples (genotype,fitness)
    new_pop = []
    path_length = 0.
    while population:
        #print(path_length)
        individual = population.pop()
        phenotype, length = develop(individual, initial, gpmap)
        if length:
            new = (individual, fitness(phenotype,target))
            new_pop.append(new)
            path_length += length
    path_length /= len(new_pop) # attention dividing by zero
    #print(path_length)
    new_pop.reverse() # reverse is optional, but easier to debug
    return new_pop, path_length


# for VECTORS only
# 4 times SLOWER than brent
#  (using sign or cSigmoid2; 5 times slower using sigmoid2)

def devo_vector(
        genotype, phenotype, gpmap=def_gpmap, base=def_base,
        noise=1e-100, noise_size=N, noise_function=random_noise,
        noise_time='before'):

    # noise at the beginning of each development step
    if noise_time == 'before':
        phenotype = noise_function(phenotype, noise, noise_size, base=base)

    # the development model: matrix multiplication (additive/linear)
    phenotype = dot(genotype, phenotype)

    # noise at the end of each development step
    if noise_time == 'after':
        phenotype = noise_function(phenotype, noise, noise_size, base=base)

    # normalization/filter/step or sigmoid function: introduces nonlinearity
    return gpmap(phenotype)

def develop(
        genotype, phenotype, gpmap=def_gpmap, noise=1e-100,
        window=T, converge_threshold=converge, devo_time=M,
        noise_size=N, noise_function=random_noise, noise_time='before',
        base=def_base, verbose=0):

    if noise_size > phenotype.size:
        noise_size = phenotype.size

    #if noise_function is not random_noise and isinstance(noise, int):
     #   noise_size = noise

    phenotypes = phenotype
    # alternative = array([phenotype]) that gives shape (1, N).
    # Using lists is also slower.
    # Fix size should be faster (no need to resize/realloc)

    for i in arange(1, devo_time + 1):

        # flip (no devo): from Shmulevich2002b (with noise)
        if (noise_time == 'shmulevich' and
            random() < 1 - (1-noise)**phenotype.size):
            phenotype = noise_function(
                phenotype, noise, noise_size, base, noise_time)

        # devo with or without noise
        else:
            phenotype = devo_vector(
                genotype, phenotype, gpmap, base,
                noise, noise_size, noise_function, noise_time)

        # store phenotypes for variance calculation
        phenotypes = append(phenotypes, phenotype) # SLOW
        # alternative: append(phenotypes, [phenotype], axis=0),
        # no need to reshape, but surprisingly slower.
        # Consider fix size arrays to speed up!

        # convergence test after 'window' steps
        if i >= window:
            phenotypes = phenotypes.reshape(i+1, -1) # shape ( , N)
            if convergence(phenotypes[-window:], converge_threshold):
                del phenotypes
                # stable: period, length, phenotype
                return 1, i+1-window, phenotype

    # no convergence
    del phenotypes
    # NOTE: Use inf because None < 2: True, and inf < 2: False.
    # BUT bool(inf) = True, and bool(None) = False - ATENTION!!
    # unstable: period, length, phenotype.
    return inf, inf, phenotype


## the following devo functions are for DECIMALS only, i.e.,
# only DISCRETE phenotypes. They are called by brent. Faster!
# note: you can't have noise in devo -> brent because it
# assumes DETERMINISTIC dynamics!!
def devo_verbose(
        genotype, phe_decimal,
        basename=def_basename, gpmap=def_gpmap, base=def_base, bits=N):

    phenotype = get_vector(phe_decimal, basename, base, bits)
    print phe_decimal, phenotype
    phenotype = dot(genotype, phenotype)
    print phenotype
    phenotype = gpmap(phenotype)
    phe_decimal = get_decimal(phenotype, basename, base, bits)
    print phe_decimal, phenotype
    return phe_decimal

def devo_fixed_state(
        genotype, phe_decimal,
        basename=def_basename, gpmap=def_gpmap, base=def_base, bits=N,
        fixed_states=None):

    phenotype = get_vector(phe_decimal, basename, base, bits)
    phenotype = dot(genotype, phenotype)
    phenotype = gpmap(phenotype)
    phenotype[fixed_states[0]] = [fixed_states[1]]
    return get_decimal(phenotype, basename, base, bits)

def devo(
        genotype, phe_decimal, basename=def_basename, gpmap=def_gpmap,
        base=def_base, bits=N, *args):
    phenotype = get_vector(phe_decimal, basename, base, bits)
    phenotype = dot(genotype, phenotype)
    phenotype = gpmap(phenotype)
    return get_decimal(phenotype, basename, base, bits)

def brent(
        f, genotype, x0, basename=def_basename, gpmap=def_gpmap, base=def_base,
        bits=N, max_power=inf, fixed_states=None, sampling_power=False):

    # main phase: search successive powers of two;
    # lam is the period (length of the cycle)
    power = lam = 1
    # f(x0) is the element/node next to x0.
    tortoise, hare = x0, f(
        genotype, x0, basename, gpmap, base, bits, fixed_states)
    while tortoise != hare:
        if power == lam:   # time to start a new power of two?
            tortoise = hare
            power *= 2
            # put a cap on development time
            if power > max_power:
                return inf, inf, power
            lam = 0
        hare = f(genotype, hare, basename, gpmap, base, bits, fixed_states)
        lam += 1

    # Find the position of the first repetition of length lambda
    # (where the cycle starts: mu)
    mu = 0
    tortoise = hare = x0
    for i in range(lam):
    # lets advance hare one period ahead of the tortoise, which is at x0;
    # we have to do this step by step
        hare = f(genotype, hare, basename, gpmap, base, bits, fixed_states)
    while tortoise != hare:
        # now that hare and tortoise are one period away from each other,
        # lets advance them together, and find mu where the tortoise stops
        tortoise = f(
            genotype, tortoise, basename, gpmap, base, bits, fixed_states)
        hare = f(genotype, hare, basename, gpmap, base, bits, fixed_states)
        mu += 1

    if sampling_power:
        return lam, mu, tortoise, power
    else:
        return lam, mu, tortoise # period, length, phe_decimal

# general function
def general_develop(
        genotype, initial, basename, gpmap, base, bits=N, noise=0,
        devo_function=devo, converge_time=T, converge_threshold=converge,
        devo_time=M, fixed_states=None,
        noise_function=random_noise, noise_time='before'):

    if not noise:
        if isinstance(initial, ndarray):
            initial = get_decimal(initial, basename, base, bits)

        period, length, phe_decimal = brent(
            devo_function, genotype, initial, basename, gpmap, base, bits,
            devo_time, fixed_states)

        phenotype = None

    else:
        if not isinstance(initial, ndarray):
            initial = get_vector(initial, basename, base, bits)

        period, length, phenotype = develop(
            genotype, initial, gpmap, noise, converge_time, converge_threshold,
            devo_time, bits, noise_function, noise_time, base)
        phe_decimal = None

    return period, length, phe_decimal, phenotype

# ORBITS
def get_orbit(
        period, genotype, phe_decimal,
        basename=def_basename, gpmap=def_gpmap, base=def_base,
        bits=N, devo_function=devo, fixed_states=None):

    # consider the option to use python sets: can't! sets are NOT indexable
    orbit = []
    while len(orbit) < period:
        phe_decimal = devo_function(
            genotype, phe_decimal, basename, gpmap, base, bits, fixed_states)
        orbit.append(phe_decimal)
    #orbit.sort() # This is wrong! Order matters!!
    # Will mantain order but start by the smallest state:
    # get the index of the smallest state
    # (smallest in terms of vector-to-decimal code)
    i = orbit.index(min(orbit))
    new = orbit[i:]
    new.extend(orbit[:i])
    return tuple(new) # tuple to be hashable and call set in set_unique_orbits()


#############################
## RECOMBINATION (haploids) #
#############################

def recombination(population):
# note: 1) if you want to keep old population call with a copy, e.g. pop[:]
    new_pop = []
    while population:
        new1, new2 = recombine(population.pop(), population.pop())
        new_pop.extend([new1, new2])
    new_pop.reverse() # reverse is optional, but easier to debug
    return new_pop

# NOTE: adapted for class individual
def recombine(ind1, ind2, rate = r):
    # segregate entire rows with probability r = 0.5 (mother+father)
    segregation = permutation(ind1.size)[:rate*ind1.size]
    # genotypes
    g1 = ind1.genotype; g2 = ind2.genotype
    g1[segregation], g2[segregation] = g2[segregation], g1[segregation]
    # Modifier allele
    #c1 = ind1.connectivity; c2 = ind2.connectivity
    #c1[segregation], c2[segregation] = c2[segregation], c1[segregation]
    #return ind1, ind2

def flat_recombine(ind1, ind2, rate = r):
    size = ind1.size**2
    segregation = permutation(size)[:rint(rate*size).astype(int)]
    g1 = ind1.genotype.flatten(); g2 = ind2.genotype.flatten()
    g1[segregation], g2[segregation] = g2[segregation], g1[segregation]
    ind1.genotype = g1.reshape((ind1.size, ind1.size))
    ind2.genotype = g2.reshape((ind2.size, ind2.size))

##############
## SELECTION #
##############

def selection(population):
# note: it randomizes the order of the new population
    viable = len(population)
    max_fitness = max(population, key = operator.itemgetter(1))[1]
    new_pop = []
    while (len(new_pop) < P):
        individual = randint(viable)
        genotype = population[individual][0]
        fitness  = population[individual][1]
        fitness /= max_fitness
        r = uniform()
        if (r < fitness):
            new_pop.append(genotype)
    return new_pop


##############
## EVOLUTION #
##############

def evolution(
        population, initial, optimum, gpmap, path_length=array(0),
        generations=G):
# note: 1) if you want to keep old population call with a copy, e.g. pop[:]
    for i in range(1, generations+1):
        #print("generation " + str(i))
        recombined = recombination(population)
        mutated    = mutation(recombined)
        # development returns [tuples (individual, fitness)], length
        developed, length = development(mutated, initial, optimum, gpmap)
        #developed = [(ind, 1) for ind in mutated] # debugging
        path_length = append(path_length, length)
        population  = selection(developed) # note: substituting population
    return population, path_length, developed


############################
## FOUNDER FUNCTIONS (OLD) #
############################

# A. Wagner's method for producing founder: integer phenotype
# (a,c) must be large: (10,0.7) ou (100,0.4)
def wagners_founder(initial=None, optimum=None):
    while(True):
        if not initial:
            initial = spin_phenotype()
        if not optimum:
            optimum = spin_phenotype()

        print(initial, optimum)

        for n in range(F):
            genotype = generate_genotype()
            phenotype, length = develop(genotype, initial, sign)
            if (normed_squared_vector_distance(phenotype, optimum) < converge):
                # NOTE: should change returning gpmap to returning period,
                # for uniformity
                return initial, optimum, genotype, length, sign

        initial, optimum = None, None

    # NOTE: develop returns period = length = inf if no convergence.
    # Should uniformize
    return 0,0,0,0,0

# Siegal's method for producing founder
# note: can return unstable individuals as well.
# Test length - NOTE: changed develop() so that it returns
# period = length = inf if no convergence
def siegals_founder(dtype=Dtype, gpmap=sigmoid2, density=c, verbose=0):

    initial = spin_phenotype(dtype=dtype)
    genotype = generate_genotype(density=density, dtype=dtype)
    period, length, phenotype = develop(
        genotype, initial, gpmap, verbose=verbose)

    # NOTE: should change returning gpmap to returning period, for uniformity
    return initial, phenotype, genotype, length, gpmap

#
def siegals_founder_decimal(
        dtype=Dtype, density=c, basename=def_basename, gpmap=def_gpmap,
        base=def_base, bits=N):

    initial  = generate_decimal_phenotype(base, bits)
    genotype = generate_genotype(bits, density, dtype=dtype)
    period, length, phenotype = brent(
        devo, genotype, initial, basename, gpmap, base, bits)

    # NOTE: returns period insted of gpmap
    return initial, phenotype, genotype, length, period


#############################
## GENERATE ENSEMBLES (OLD) #
#############################

# Wagner's ensemble generation
def generate_ensemble(founder_function=siegals_founder, size=n_runs):
    ensemble = []
    while len(ensemble) < size:
        initial, optimum, founder, length, gpmap = founder_function()
        if length:
            ensemble.append((initial, optimum, [founder]))
    #model_name = founder_function.__name__
    model_name = get_model_name(gpmap)
    middle_name = parameters
    filename = ('ensemble%(size)s%(model_name)s%(middle_name)s_before.npy' %
                locals())
    save_ensemble(ensembles_dir + filename, ensemble)
    return ensemble

#
def generate_unstable_ensemble(founder_function=siegals_founder, size=n_runs):
    ensemble = []
    while len(ensemble) < size:
        initial, optimum, founder, length, gpmap = founder_function()
        if not length:
            ensemble.append((initial, optimum, [founder]))
    #model_name = founder_function.__name__
    model_name = get_model_name(gpmap)
    middle_name = parameters
    filename = ('unstable%(size)s%(model_name)s%(middle_name)s_before.npy' %
                locals())
    save_ensemble(ensembles_dir + filename, ensemble)
    return ensemble

# evolve each population founded by each member of the initial ensemble
def evolve_ensemble(ensemble, filename, gpmap = sigmoid2):
    new_ensemble = []
    for triplet in ensemble:
        initial, optimum, founder = triplet
        ## all individuals are identical in the starting population
        population = [founder[0] for i in range(P)]
        developed  = evolution(population, initial, optimum, gpmap)[2]
        # developed is a list of tuples (ind, fitness): we just want ind
        population = [pair[0] for pair in developed]
        new_ensemble.append((initial,optimum,population))
    # save file
    save_ensemble(filename, ensemble)
    return new_ensemble

#
def compare_ensembles(sample_a, sample_b):
    equal = True
    for a, b in zip(sample_a, sample_b):
        if not len(a[2]) == len(b[2]): return 0
        equal *= a[0].all() == b[0].all() ######## ATENCAO!!! ISTO ESTA' MAL!!!! #####################
        equal *= a[1].all() == b[1].all() ######## ATENCAO!!! ISTO ESTA' MAL!!!! #####################
        for ind1, ind2 in zip(a[2], b[2]):
            equal *= ind1.all() == ind2.all() ######## ATENCAO!!! ISTO ESTA' MAL!!!! #####################
    return equal


################
## NOISE (OLD) #
################

def noise_phenotypes(
        triplet, noise, gpmap=sigmoid2, samples=n_perturb, noise_size=N):
    initial, optimum, population = triplet
    individual = population[0]
    # develop the same genotype n times and return each (phenotype, length) in
    # a list of lists
    phenotypes = [develop(
        individual, initial, gpmap, noise, noise_size) for i in range(samples)]
    #print '', optimum, '\n', array(phenotypes)
    return phenotypes

def noise_ensemble(
        ensemble, noise, gpmap=sigmoid2, samples=n_perturb, noise_size=N,
        stable='s'):
    ensbl_phenotypes = [noise_phenotypes(
        triplet, noise, gpmap, samples, noise_size) for triplet in ensemble]
    # save file
    mod_name = get_model_name(gpmap)
    size = len(ensemble)
    mid_name = parameters
    filename = ('%(stable)s-phenotypes%(samples)s_ensbl%(size)s%(mod_name)s'
                + '%(mid_name)s_%(noise_size)snoise%(noise)s.npz' % locals())
    file = open(data_dir + filename, 'wb')
    # python object format: list of lists of arrays - too big!
    #cPickle.dump(ensbl_phenotypes, file)
    # compressed numpy.array format: 3D array of shape (size, samples, N)
    #savez(file, ensbl_phenotypes)
    file.close()
    return ensbl_phenotypes

def noise_interval(
        ensemble, max_, step, gpmap=sigmoid2, samples=n_perturb, stable='s'):

    noises = arange(step, max_, step)
    n = len(noises)
    size = len(ensemble)

    if samples * n > 3000:
        # using a Generator expression instead of a List comprehension because
        # of memory issues: NOT returning list
        all_phenotypes = (noise_ensemble(
            ensemble, noise, gpmap, samples, stable) for noise in noises )
        # not saving file because it's too BIG! (use independent files for each
        # noise instead (ensbl_phenotypes))

    else:
        all_phenotypes = [noise_ensemble(
            ensemble, noise, gpmap, samples, stable) for noise in noises ]
        # save file
        mod_name = get_model_name(gpmap)
        mid_name = parameters
        filename = ('%(stable)s-phenotypes%(samples)sx%(n)s_ensbl%(size)s'
                    + '%(mod_name)s%(mid_name)s_noise-step%(step)s.npz' %
                    locals())
        file = open(data_dir + filename, 'wb')
        # python object format: list of lists of arrays - too big!
        #cPickle.dump(all_phenotypes, file)
        # compressed numpy.array format: 3D array of shape (n, size, samples, N)
        savez(file, all_phenotypes)
        file.close()
        return all_phenotypes

# for the whole ensemble for a particular level of noise
def noise_stability(ensbl_phenotypes):
    n_samples = float(len(ensbl_phenotypes[0]))
    return [(nonzero(phenotypes)[0].size / N) / n_samples
            for phenotypes in ensbl_phenotypes]

# for the whole noise interval, print summary statistics
def noise_stable_stats(all_phenotypes, step=0.05):
    stable = [noise_stability(phenotypes) for phenotypes in all_phenotypes]
    x = []; y1 = []; y2 = []
    for i, noise in enumerate(stable):
        x.append (step * (i+1))
        y1.append(mean(noise))
        y2.append(std(noise))
    errorbar(x, y1, y2)
    return stable

# further develop noisy (probabilistic) phenotypes without noise (deterministic)
def develop_phenotypes(ensemble, all_phenotypes, gpmap=sigmoid2, times=1):
# note: 1) if you only have one ensemble (1 noise value)
# you should call with [ensbl_phenotypes]
    genotypes = [triplet[2][0] for triplet in ensemble]
    for ensbl_phenotypes in all_phenotypes:
        for phenotypes, genotype in zip(ensbl_phenotypes, genotypes):
            for phenotype in phenotypes:
                phenotype = develop(
                    genotype, phenotype, gpmap, devo_time=times)[0]
    return all_phenotypes # same shape (n = len(noises), ensbl_size, samples, N)

def compare_phenotypes(phenotypes, optimum):

    # get number of stable phenotypes, i.e., not zeros(N)
    stable = nonzero(phenotypes)[0].size / N
    # Attention: watch out for zeros in stable phenotypes.
    # Could control requiring % N == 0
    # alternatives: len(where(x.sum(axis=1) == 0)[0]) or
    # len(where(x.prod(axis=1) == 0)[0])

    # get number of robust phenotypes, i.e., equal (< converge) to
    # target (optimum) phenotype
    equals = [(normed_squared_vector_distance(phenotype, optimum) < converge)
               for phenotype in phenotypes]
    n_equal = equals.count(1)

    # map to list of lists while turning equal (< converge) to
    # really equal (==) with function sign
    l_phenotypes = map(list, sign(phenotypes))
    # get unique elements of list of list phenotypes
    # ("recursive" list comprehension)
    unique = [x for x in l_phenotypes if x not in locals()["_[2]"]] # or "_[1]"
    # notice: forced to use unique for lists because
    # unique for arrays only works in 1D

    # some controls:
    l_optimum = list(sign(optimum))
    print l_phenotypes.count(l_optimum) == n_equal
    if n_equal: unique.remove(l_optimum)
    '''
    if l_phenotypes.count(l_optimum) == n_equal:
        if n_equal: unique.remove(l_optimum)
    else:
        print 'erro 1'
        return 0
        '''
    unstable = len(phenotypes) - stable
    print l_phenotypes.count(N*[0]) == unstable
    if unstable: unique.remove(N*[0])
    '''    if l_phenotypes.count(N*[0]) == unstable:
        if unstable: unique.remove(N*[0])
    else:
        print 'erro 2'
        return 0
        '''
    # compute number of different phenotypes
    different = len(unique)
    # count the number of occurrences of each phenotyope in
    # 1:1 correspondence with unique
    count = [l_phenotypes.count(x) for x in unique]

    # more controls
    print n_equal + int(sum(count)) == stable
    '''
    if n_equal + int(sum(count)) is stable:
        return stable#, n_equal, unique, count
    else:
        print 'erro 3'
        return 0
        '''
#
def noise_stats():
    # get target phenotypes from original ensemble in 1:1 correspondence with
    # phenotypes
    targets = [triplet[1] for triplet in ensemble]
    stable = [compare_phenotypes(phenotype, optimum)
              for phenotype, optimum in zip(phenotypes, targets)]

    stable = float(n - unstable) # stable + unstable = n
    total_different = stable - equal # equal + total_different = stable
    if different: different /= total_different
    if equal: equal /= stable
    return stable/n, equal, different

    stable = [row[0] for row in stats]
    equal  = [row[1] for row in stats]
    diff   = [row[2] for row in stats]
    return mean(stable), std(stable), mean(equal), std(equal), mean(diff), std(diff)

'''
stable_mean = [row[0] for row in stats]
stable_error = [row[1] for row in stats]
errorbar(x, stable_mean, stable_error)
title('stable / n for independent noise for all genes')
xlabel('gaussian noise width')
'''

#####################
## ROBUSTNESS (OLD) #
#####################

## evaluate fraction of mutations that attain optimum target: Wagner's Fig. 2 a)
def mutational_robustness(ensemble, filename):
    #initial, optimum, founder, length, gpmap = wagners_founder()
    robustness = []
    for triplet in ensemble:
        initial, optimum, population = triplet
        pop_size = len(population)
        n_mutations = n_perturb / pop_size
        n = 0
        for genotype in population:
            for i in range(n_mutations):
                mutant = mutate(genotype, 1)
                phenotype = array(develop(mutant, initial, sign)[0])
                if phenotype.all() == optimum.all():
                    n+=1 ######## ATENCAO!!! ISTO ESTA' MAL!!!! ################
        print n
        robustness.append(float(n)/n_perturb)
    save_list(filename, robustness)
    return robustness

## plot fraction of mutations that attain optimum target: Wagner's Fig. 2 a)
def robustness_evolution(ensemble=[], gpmap=sigmoid2, size=n_runs, before=[]):
    #model_name = founder_function.__name__
    model_name = get_model_name(gpmap)
    filename = '%(size)s%(model_name)s%(parameters)s' % locals()
    y  = 'mutational robustness'
    x0 = 'before evo'; x1 = 's = 0.1'; x2 = 's = 1'; x3 = 's = inf'

    ## before
    if not ensemble:
        file = ensembles_dir + 'ensemble%(filename)s_before.dat' % locals()
        try:
            f = open(file, 'r')
        except IOError:
            ensemble = generate_ensemble(
                get_founder_function(gpmap), size, file)
        else:
            ensemble = load_ensemble(file, size)
    if not before:
        file = data_dir + 'stability%(filename)s_before.csv' % locals()
        try:
            f = open(file, 'r')
        except IOError:
            before = mutational_robustness(ensemble, file)
        else:
            before = load_list(file)

    ## evolution
    global s

    s = 0.1
    file = ensembles_dir
    file += 'ensemble%(filename)s_after%(G)s_P%(P)s_s%(s)s.dat' % locals()
    try:
        f = open(file, 'r')
    except IOError:
        new1 = evolve_ensemble(ensemble, file, gpmap)
    else:
        new1 = load_ensemble(file, size)

    file = data_dir
    file += 'stability%(filename)s_after%(G)s_P%(P)s_s%(s)s.csv' % locals()
    try:
        f = open(file, 'r')
    except IOError:
        after1 = mutational_robustness(new1, file)
    else:
        after1 = load_list(file)

    plot(before, after1, 'p')
    axis([0, 1, 0, 1])

    #r.boxplot(before, after1, ylab=y, names=[x0, x1])
    figure()
    boxplot([before, after1])

    s = 1
    new2 = evolve_ensemble(ensemble)
    save_ensemble(ensembles_dir + 'ensemble' + filename + '_after_s1.dat', new2)
    after2 = mutational_robustness(new2)
    save_list(data_dir + 'stability_' + filename + '_after_s1.csv', after2)

    figure()
    plot(before, after2, 'p')
    axis([0, 1, 0, 1])
    #r.boxplot(before, after1, after2, ylab=y, names=[x0, x1, x2])
    close()
    boxplot([before, after1, after2])

    s = inf
    new3 = evolve_ensemble(ensemble)
    fname = ensembles_dir + 'ensemble' + filename + '_after_sinf.dat'
    save_ensemble(fname, new3)
    after3 = mutational_robustness(new3)
    fname = data_dir + 'stability_' + filename + '_after_sinf.csv'
    save_list(fname, after3)

    #r.boxplot(before, after1, after2, after3, ylab=y, names=[x0, x1, x2, x3])
    close
    boxplot([before, after1, after2, after3])

    return before, after1, after2, after3

## plot mean path length against number of generations: Siegal's Fig. 2
def path_length_evolution(ensemble=[]):
    global T
    T = 2
    print T, s
    path_length = array(0)
    for run in range(n_runs):
        ## Founder
        # Wagner's
        #initial, optimum, founder, length, gpmap = wagners_founder()
        # Siegal's ## ATENCAO: MUDEI O RETORNO DA FUNCAO siegals_founder.
        # Mas ja agora, tb deveria testar length de wagners
        initial, optimum, founder, length, gpmap = siegals_founder()
        print initial, '\n', optimum, '\n', founder
        print(length)
        if not path_length.any():
            path_length = length
        else:
            path_length = append(path_length, length)
        # all individuals are identical in the starting population
        population  = [founder for i in range(P)]
        # or for debugging
        #pop = [wagners_founder(initial,optimum)[2],
         #      wagners_founder(initial,optimum)[2]]

        path_length = evolution(
            population, initial, optimum, gpmap, path_length)[1]
        path_length = path_length.reshape(run+1, G+1)
        #plot(range(G+1), path_length[run], 'p')

    avg = mean(path_length,0)
    plot(range(G+1), avg)
    return avg


### RUN
'''s = 1
avg2 = path_length_evolution()
s = inf
avg3 = path_length_evolution()
'''
### plot average path length against network size N
'''n_runs = 60
pop_range = range(4,21,2)
avg = []
for pop_size in pop_range:
    N = pop_size
    initial_path_length = 0.
    for i in range(n_runs):
    #initial_path_length += wagners_founder()[3] # = 11.92
        initial_path_length += siegals_founder()[3] # = 11.46 ~ 11.82

    avg.append(initial_path_length/n_runs)
plot(pop_range,avg,'p')
'''
#print(init_pheno, opt_pheno, genotype)
#stable = devo(genotype, phenotype)

'''
ipython -pylab
import sys
sys.path.append('~/model')
import model
from model import *

reload(model)
from model import *

## some helpfull functions
for i in range(P):
    print selected[i][0][0]

## Obsoletos

## keep only the first member of each tuple in the list: the genotype (removes fitness)
def parse(list):
    new_list=[]
    for pair in list:
        new_list.append(pair[0])
    return new_list

## alternative
def sigmoid2(x):
    for i in range(N):
        x[i] = 1 / ( 1 + exp(-a*x[i]) )
    return x

# numa so linha
genotype[array(randint(0,N*N,c))]=standard_normal(c)

# para confirmar que temos c elementos non-zero
print len(nonzero(genotype)[0]) is int(c*N*N)

n=0
while n<c:
    i=randint(N*N)
    if(a[i] == 0):
        a[i]=1
        n+=1
print a
'''

def thislist():
    """Return a reference to the list object being constructed by the
    list comprehension from which this function is called. Raises an
    exception if called from anywhere else.
    """
    import sys
    d = sys._getframe(1).f_locals
    nestlevel = 1
    while '_[%d]' % nestlevel in d:
        nestlevel += 1
    return d['_[%d]' % (nestlevel - 1)].__self__


#########
## MAIN #
#########

## for profiling purposes only
def main():
    '''file = open(pop_dir + 'founder.dat')
    initial, optimum, founder, length, gpmap = cPickle.load(file)
    gpmap = cSigmoid2
    population  = [founder for i in range(P)]
    #new_pop, path_length = development(population, initial, optimum, gpmap)
    new_pop, path_length = evolution(population, initial, optimum, gpmap, generations = 50)[:2]
    print new_pop[:1], path_length
    '''
    '''filename  = ensembles_dir + 'unstable500_siegals100_N10_c1_T10_before.npy'
    uensemble = load_ensemble(filename, 500)
    return noise_interval(uensemble, 2, 0.05, sigmoid2, 10, 'u')
    '''
    '''file = open(data_dir + 'phenotypes.dat','rb')
    phenotypes = cPickle.load(file)
    file.close()
    for i, phenotype1 in enumerate(phenotypes):
        for phenotype2 in phenotypes[i+1:]:
            vec_dist3(phenotype1, phenotype2)
    '''
    '''dim =  2
    min = -1.
    max =  1.
    n   = 1e4

    jumps, base   = get_equalprobable_states(dim, min, max)

    basename = get_basename(base)
    if basename not in basemaps: load_base(basename, base)
    #basename = 'gray'
    print(basename)

    #[get_vector(generate_decimal_phenotype(base), basename) for i in range(1e7)]
    [get_decimal(spin_phenotype(), basename) for i in range(1e6)]
    '''
    '''#gpmap = sign
    step = get_step_function(dim, 'cstep')
    def gpmap(x): return step(x, jumps, base)
    #def gpmap(x): return cstep(x, jumps, base)
    mymap = False
    devo_function = devo if not mymap else devo_nodict
    list(initial_sampling_stat(n, basename = basename, gpmap = gpmap, base = base, devo_function = devo_function))
    '''
    '''dim = 2
    bit = 30
    noise = 0 #1e-100
    devo_time = 'mean'
    devo_time = inf#devo_times[devo_time]['%d' %dim]['%d' %bit] #4000#M
    samples = 1e2 #50
    filename = data_dir + 'period_length_N%d_noise%s_M%g.dat' %(bit, noise, devo_time)
    data = [test_devo_funtions2(bit, noise, devo_time) for i in xrange(int(samples))]
    save_file(filename, data)
    '''
    #import runs.run_phenotypes_prob2
    #reload(runs.run_phenotypes_prob2)

    #mutants = get_single_mutants(generate_genotype(), 1e5)

    #bit=100;alpha=.02;gamma=.8;delta_in=.2;delta_out=0;samples=300
    #array([generate_scale_free_genotype(bit, alpha, 1 - alpha - gamma, gamma, delta_in, delta_out).in_degree().values() for i in range(samples)]).mean()

    bits = N
    [[int_decimal_hamming_distance(i, j) for j in range(2**bits)] for i in range(2**bits)]
    #[[int_vector_hamming_distance(get_vector(i, def_basename, def_base, bits), get_vector(j, def_basename, def_base, bits)) for j in range(2**bits)] for i in range(2**bits)]

'''
import cProfile
import pstats
filename_pstats = 'data/profiling/model2.pstats'
cProfile.run('import src.model2; src.model2.main()', filename_pstats)
stats = pstats.Stats(filename_pstats).sort_stats('cumulative').print_stats(15)

cProfile.run('reload(src.model2); src.model2.main()', filename_pstats)
'''


#########################
## EXPERIMENTS AND RUNS #
#########################


#######################################
## EIGENVECTORS AND VALUES EXPERIMENT #
#######################################

## 2d
def phi(vector):
    return arctan(vector[1] / vector[0])

def devo2d(genotype, phenotype):
    global N
    N = 2
    epolar = array([phi(vector) for vector in eig(genotype)[1]])
    for i in range(20):
        phenotype = dot(genotype, phenotype)
        print epolar - phi(phenotype)

## 3d
def theta(vector):
    return arctan( sqrt(vector[0]*vector[0]+vector[1]*vector[1]) / vector[2] )

def epolar(genotype):
    evalues, evectors = eig(genotype)
    epolar = [array([theta(vector), phi(vector)]) for vector in evectors]
    return array(epolar)

def devo3d(genotype = None, phenotype = None):
    global N
    N = 3
    genotype = genotype if genotype else generate_genotype()
    phenotype = phenotype if phenotype else spin_phenotype()
    evectors = epolar(genotype)
    for i in range(50):
        phenotype = dot(genotype, phenotype)
        vector = array([theta(phenotype), phi(phenotype)])
        print evectors - vector

## n-dimensional
def eigen_devo(genotype = None, phenotype = None):
    genotype = genotype if genotype else generate_genotype()
    phenotype = phenotype if phenotype else spin_phenotype()
    initial = phenotype.copy()

    eval, evect = eig(genotype) ;print eval#, evect
    print linalg.solve(evect,phenotype)
    '''
    for i in range(10):
        print i
        x = linalg.solve(evect,phenotype) ;print x
        for val in eval:
            print x/pow(val,i)
        phenotype = dot(genotype, phenotype)
        '''
    return develop(genotype, initial, sign)


######################################
## INITIAL STATE SAMPLING AND STATS ##
######################################

'''def initial_sampling_stats(samples = n_perturb, verbose = 0, founder_function = siegals_founder):
    length = [founder_function(verbose = verbose)[3] for i in range(samples)]
    stable = nonzero(length)[0].size / float(samples)
    print stable
'''
def initial_sampling_stat(samples = n_perturb, density = c, basename = def_basename, gpmap = def_gpmap, base = def_base, devo_function = devo, mymap = False, bits = N):
    return (brent(devo_function, generate_genotype(density = density), generate_decimal_phenotype(base, bits, mymap), basename, gpmap, base) for i in xrange(samples))

def initial_sampling_period(samples = n_perturb, density = c, basename = def_basename, gpmap = def_gpmap, base = def_base, devo_function = devo, zero = 0):
    sampling = initial_sampling_stat(samples, density, basename, gpmap, base, devo_function)
    # discount the all-zero state: phe_decimal == 0
    return [stat[0] for stat in sampling if stat[2] is not zero] # brent[0] == period, brent[1] == length, brent[2] == phe_decimal

def stability_vs_c(samples = n_perturb, basename = def_basename, gpmap = def_gpmap, base = def_base, devo_function = devo, zero = 0, period = 1, less = False, step = 0.1):
    step = 0.1 if step < 0.1 else step
    if less:
        return [where(array(initial_sampling_period(samples, density, basename, gpmap, base, devo_function, zero)) <= period)[0].size / float(samples) for density in arange(step, 1 + step, step)]
    else:
        return [initial_sampling_period(samples, density, basename, gpmap, base, devo_function, zero).count(period) / float(samples) for density in arange(step, 1 + step, step)]

# SAME STABILITY INCLUDES CYCLES!
def stability_vs_c2(
        samples=n_perturb,
        basename=def_basename, gpmap=def_gpmap, base=def_base,
        devo_function=devo, sel_period=1, noise=None, discrete=True):

    # -1 not to filter (there's no negative phe_decimal in any map).
    # Only filter if state = 0 is in base
    zero = get_decimal(zeros(N), basename, base, bits) if 0. in base else -1

    stability = []
    for density in arange(.1, 1.1, .1):
        same_stability = 0
        for i in xrange(samples):
            genotype = generate_genotype(density=density)
            if discrete and not noise:
                period, length, phe_decimal = brent(
                    devo_function, genotype, generate_decimal_phenotype(base),
                    basename, gpmap, base)
                # fixpoints
                if period == sel_period == 1 and phe_decimal is not zero:
                    same_stability += 1
                # cycles
                if (period == sel_period  > 1 and
                    zero not in get_orbit(
                        period, genotype, phe_decimal, basename, gpmap, base)):
                    same_stability += 1
            else:
                period, length, phenotype = develop(
                    genotype, generate_vector_phenotype(base), gpmap, noise)
                # fixpoints
                if period == sel_period == 1 and phenotype.any() != 0:
                    same_stability += 1
        stability.append(same_stability / float(samples))
    return stability


############
## Elhanan #
############

#
def compare_stability_phenotype(
        source_period, target_period, source_phenotype, target_phenotype):

    same_stability = 0; same_phenotype = 0

    # target ("wild-type") and perturbed ("mutant") phenotypes are
    # BOTH Fixpoints
    if target_period == 1 and source_period == 1:

        same_stability = 1

        # Equal
        if target_phenotype == source_phenotype: same_phenotype = 1

    # target ("wild-type") and perturbed ("mutant") phenotypes are BOTH Cycles
    elif target_period > 1 and source_period > 1:

        same_stability = 1

        # Equal
        if target_phenotype == source_phenotype: same_phenotype = 1

    #print target_period, source_period, same_stability
    #print target_phenotype, source_phenotype, same_phenotype

    return same_stability, same_phenotype

# for Elhanan to build a graph like in Sevim's Fig. 3
def graph(samples = n_perturb, step = .2):
    for density in arange(step, 1.01, step):
        genotypes = [generate_genotype(density = density)
                     for i in range(samples)]
        for i, genotype in enumerate(genotypes):
            phenotypes = array([
                (phe_decimal, devo(genotype, phe_decimal))
                for phe_decimal in range(pow(def_base.size, N))])
            savetxt(
                data_dir + 'elhanan-graph/c%.1f-%s.txt' %(density, i+1),
                phenotypes, '%d, %d', '\n')


#############################
## PROFILING DEVO FUNCTIONS #
#############################

#
def test_devo_funtions(
        genotype_function=generate_genotype, phenotype_function=spin_phenotype,
        devo_function=devo, basename=def_basename, gpmap=def_gpmap,
        base=def_base):

    genotype = genotype_function()
    phenotype = phenotype_function()
    vector_f = list(develop(genotype, phenotype, gpmap))
#    decimal_f = list(brent(
 #       devo_function, genotype, get_decimal(phenotype, def_basename),
  #      basename, gpmap, base))
    #vector_f[2] = get_decimal(vector_f[2], basename)
    return vector_f#, decimal_f

#
def test_devo_funtions2(
        bits=N, noise=0, devo_time=M, dim=2, min_=-1., density=c, binary=False,
        converge_time=T, converge_threshold=converge):
    devo_function = devo
    genotype_function  = generate_genotype if not binary else spin_genotype
    if not noise:
        phenotype_function = generate_decimal_phenotype
    else:
        phenotype_function = generate_vector_phenotype
    base, basename, jumps = get_base(dim, bits, min_)
    gpmap = get_gpmap(base, basename, density, bits, binary, jumps)
    initial, genotype = (phenotype_function(base, bits),
                         genotype_function(bits, density))
    period, length, phe_decimal, phenotype = general_develop(
        genotype, initial, basename, gpmap, base, bits, noise, devo_function,
        converge_time, converge_threshold, devo_time)

    return period#, length


################
## PATH LENGTH #
################

# studying a cut-off for devo time
def path_length_vs_N(
        bits=N, dim=2, min_=-1., density=c, noise=0, samples=1e3, binary=False,
        converge_time=T, converge_threshold=converge, devo_time=M):

    devo_function = devo
    genotype_function = generate_genotype if not binary else spin_genotype
    if not noise:
        phenotype_function = generate_decimal_phenotype
    else:
        phenotype_function = generate_vector_phenotype
    base, basename, jumps = get_base(dim, bits, min_)
    gpmap = get_gpmap(base, basename, density, bits, binary, jumps)

    filename = data_dir + 'path_leng_vs_N%d.txt' %bits

    lengths = []
    for i in xrange(int(samples)):
        initial, genotype = (phenotype_function(base, bits),
                             genotype_function(bits, density))
        period, length, phe_decimal, phenotype = general_develop(
            genotype, initial, basename, gpmap, base, bits, noise,
            devo_function, converge_time, converge_threshold, devo_time)
        if period == 1:
            lengths.append(length)
            if bits > N_dict:
                save_txt3(filename, lengths, 'w')

    save_txt3(filename, lengths, 'w')

    return array(lengths)

#
def path_length_vs_N2(
        individuals, basename, gpmap, base, bits=N, noise=0,
        devo_function=devo, converge_time=T, converge_threshold=converge,
        devo_time=M, filter=None, sampling=False, datapoints=None,
        filename=None):

    samples = 0
    lengths = []
    for initial, genotype in individuals:

        # develop and test stability: period = 1 if fixed point else 0
        period, length, phenotype = test_stability(
            genotype, initial, basename, gpmap, base, bits, noise,
            devo_function, converge_time, converge_threshold, devo_time, filter)

        if period:
            lengths.append(length)

        # save (i, n) to file every data point
        samples +=1
        if sampling and samples in datapoints:
            save_txt3(filename, lengths, 'w')

    return array(lengths)

# studying a cut-off for devo time in brent:
# the relationship between (maximum) power and length (transient time)
def power_vs_mu(
        bits=N, dim=2, min_=-1., density=c, samples=1e3, binary=False):

    devo_function = devo
    genotype_function = generate_genotype if not binary else spin_genotype
    phenotype_function = generate_decimal_phenotype
    base, basename, jumps = get_base(dim, bits, min_)
    gpmap = get_gpmap(base, basename, density, bits, binary, jumps)

    filename = data_dir + 'power_vs_mu_N%d.txt' %bits

    for i in xrange(int(samples)):

        initial, genotype = (
            phenotype_function(base, bits), genotype_function(bits, density))
        period, length, phenotype, power = brent(
            devo_function, genotype, initial, basename, gpmap, base, bits,
            sampling_power = True)

        if period == 1:
            save_txt2(filename, (power, length))


##############
## STABILITY #
##############

# develop one individual and test period for fixed point or cycle
def test_stability(
        genotype, initial, basename, gpmap, base, bits=N, noise=0,
        devo_function=devo, converge_time=T, converge_threshold=converge,
        devo_time=M, _filter=None, fixed_states=None):

    period, length, phe_decimal, phenotype = general_develop(
        genotype, initial, basename, gpmap, base, bits,
        noise, devo_function, converge_time, converge_threshold,
        devo_time, fixed_states)

    # fixpoints
    if not noise:
        # _filter: 0 > None is True, 0 > 0 is False
        if period == 1 and phe_decimal != _filter:
            return period, length, phenotype

    else:
        if period == 1 and linalg.norm(phenotype) > _filter:
            return period, length, phenotype

    # cycles
    return 0, 0, 0

def get_stability_and_phenotypes(
        n, phenotypes, samples, genotype, initials, basename, gpmap, base,
        bits=N, noise=0,
        devo_function=devo, converge_time=T, converge_threshold=converge,
        devo_time= M, _filter=None, experiment='stability', sampling=False,
        datapoints=None, filename=None, fixed_states=None):

    for initial in initials:

        # develop and test stability: period = 1 if fixed point else 0
        period, length, phenotype = test_stability(
            genotype, initial, basename, gpmap, base, bits,
            noise, devo_function, converge_time, converge_threshold,
            devo_time, _filter, fixed_states)
        n += period

        # add phenotypes to stack
        if period and experiment == 'phenotypes':
            if noise and 'cSigmoid' not in (gpmap.__name__ and gpmap.func_doc):
                phenotype = get_decimal(phenotype, basename, base, bits)
            phenotypes.append(phenotype)

        # save (i, n) to file every data point
        samples +=1
        if sampling and samples in datapoints and experiment != 'distribution':
            save_txt2(filename, (samples,n))

    return n, phenotypes, samples

# sample a population of individuals; from run_stability_sample.py
def stability_sample(
        individuals, initials, basename, gpmap, base, bits=N, noise=0,
        devo_function=devo, converge_time=T, converge_threshold=converge,
        devo_time=M, filter=None, get_initials=True, experiment='stability',
        sampling=False, datapoints=None, filename=None, fixed_states=None):

    # load functions
    a = get_activating_fraction
    d = diag
    t = trace

    n = 0
    samples = 0
    phenotypes = []
    stability_dist = []

    for individual in individuals:

        if get_initials:
            initial, genotype = individual
            initials = [initial]
        else:
            genotype = individual

        n, phenotypes, samples = get_stability_and_phenotypes(
            n, phenotypes, samples, genotype, initials,
            basename, gpmap, base, bits, noise,
            devo_function, converge_time, converge_threshold,
            devo_time, filter, experiment, sampling, datapoints,
            filename, fixed_states)

        if experiment == 'distribution':
            stability_dist.append((t(genotype), a(d(genotype)), n))
            n = 0

    return n, phenotypes, stability_dist

# sample the phenotypes (including orbits) of a population of individuals
def phenotype_sample(
        individuals, initials, basename, gpmap, base, bits=N, noise=0,
        devo_function=devo, converge_time=T, converge_threshold=converge,
        devo_time=M, get_initials=True, experiment='phenotypes', sampling=False,
        datapoints=None, filename=None, fixed_states=None, _get_orbit=True,
        only_unique=False, noise_function=random_noise, noise_time='before'):

    L = base.size**bits #sizes[str(bits)][str(base.size)]['L']
    if not base[0]:
        # all zero state is impossible to converge with def gpmap for min = 0.
        L -= 1
    n = 0 # stability
    periods = [] # experiment == 'distribution'

    # try to load saved files and pick up from there
    try:
        # filename is myrun.temp_file
        # skips first line ('# %s\n'%samples)
        phenotypes = loadtxt(filename.replace('.tsv', '.txt'), int)
    except IOError:
        phenotypes = []
        samples = 0
    else:
        phenotypes = list(phenotypes)
        samples = int(load_last_line_of_txt(filename)[0])

    for individual in individuals:

        if get_initials:
            initial, genotype = individual
            initials = [initial]
        else:
            genotype = individual

        for initial in initials:

            # develop and get period
            period, length, phe_decimal, phenotype = general_develop(
                genotype, initial, basename, gpmap, base, bits, noise,
                devo_function, converge_time, converge_threshold, devo_time,
                fixed_states, noise_function, noise_time)
            periods.append(period)

            # get fixed points
            if period == 1:
                if (noise and 'cSigmoid' not in
                    (gpmap.__name__ and gpmap.__doc__)):
                    phenotype = get_decimal(phenotype, basename, base, bits)
                else:
                    phenotype = phe_decimal
                # for uniform output (get_orbit returns tuples)
                if _get_orbit and not noise:
                    phenotype = (phenotype,)
                n +=1

            # get orbits
            else:
                if _get_orbit and not noise:
                    phenotype = get_orbit(
                        period, genotype, phe_decimal, basename, gpmap, base,
                        bits, devo_function, fixed_states)
                else:
                    phenotype = None#()

            if phenotype is not None:
                phenotypes.append(phenotype)

            samples +=1
            unique_phenotypes = unique(phenotypes)

            if sampling and samples in datapoints:

                # save period distribution to file every data point
                # develop(noise) returns period = 1 or inf, so cannot
                # study period distribution
                if 'distribution' in experiment and not noise:
                    #data = append(
                     #   [unique_hist(periods)], unique(periods).reshape(1,-1),
                      #  axis=0)
                    data = unique(periods), bincount(periods)[1:]
                    savetxt(filename, data, '%d')

                if experiment == 'phenotypes':
                    save_phenotype_sample(
                        filename, unique_phenotypes, phenotypes, samples,
                        only_unique)

            if unique_phenotypes.size == L and experiment == 'phenotypes':
                save_phenotype_sample(
                    filename, unique_phenotypes, phenotypes, samples,
                    only_unique)
                return n, phenotypes

    return n, phenotypes


# map = unique(x), unique(y)
def map_x_to_y(data, map):
    x_map = list(map[0])
    y_map = map[1]
    return array([y_map[x_map.index(x)] for x in data])

def order_x_by_y(x, y):
    x = array(x)
    y = array(y)
    return x[lexsort((x,y))[::-1]]

# int:   numpy.bincount
# float: scipy.stats.itemfreq
def unique_hist(data):
    count = list(data).count
    return [count(x) for x in unique(data)]

def unique_hist_array(data):
    return 0

def get_entropy(x):
    return x * log(x)

def evaluate_entropy(phenotypes, normalize=True):
    unique_phenotypes = unique(phenotypes)
    n = float(len(phenotypes))
    count = phenotypes.count
    entropy = sum((get_entropy(count(x) / n) for x in unique_phenotypes))
    norm = -log(unique_phenotypes.size) if normalize and entropy else 1
    return entropy / norm

# general funtion
def evaluate_entropy2(data, normalize='unique'):
    _unique = set(data)
    n = float(len(data))
    count = list(data).count
    entropy = sum((get_entropy(count(x) / n) for x in _unique))
    #pk = [count(x) / n for x in _unique]
    #entropy = stats.entropy(pk) #from scipy import stats
    if entropy:
        if isinstance(normalize, int):
            return entropy / -log(normalize)
        if normalize == 'all':
            return entropy / -log(n)
        if normalize == 'unique':
            return entropy / -log(len(_unique))
        else:
            return entropy
    else:
        return 0


def get_number_cycling_genes(
        samples, basename=def_basename, gpmap=csign, base=def_base, bits=N):

    cycling_genes = [[] for i in range(100)]

    for i in range(samples):
        genotype = spin_genotype()
        initial =  generate_decimal_phenotype(base)

        period, length, phe_decimal = general_develop(
            genotype, initial, basename, gpmap, base)[:3]

        if period == 1:
            cycling_genes[period].append(0)

        else:
            orbit = get_orbit(
                period, genotype, phe_decimal, basename, gpmap, base)

            all_diffs = [
                binary_repr(
                    bitwise_xor(decimal1, decimal2), width=bits)
                for decimal1, decimal2 in zip(orbit, orbit[1:])]

            all_diffs = array(
                [decimal_str_to_vector(decimal_str, (0,1))
                 for decimal_str in all_diffs])

            cycling_genes[period].append(
                all_diffs.sum(axis=0).nonzero()[0].size)

    return cycling_genes

def find_bias(p, samples=1e6):
    def f(x): return p - get_activating_fraction(sign(normal(x, 1, samples)))
    return optimize.brentq(f, -3, 3)


############################
## RegulonDB's E. coli TRN #
############################

def get_network(
        short_name, filters=[], knockouts=[], random=False, verbose=False):

    fnames = net_short_names[short_name]

    # read first network
    network = read_network(fnames[0])

    # convert TF name to gene if merging tf-tf with other files
    # ATTENTION: 'network_tf_tf' has always to be the first file if multiple
    # fnames if 'meta' in short_name:
    network = parse_network(network, tf_to_gene)

    # read other networks, if any
    for fname in fnames[1:]:
        network = append(network, read_network(fname), axis=0)

    # parse filters = ['(obsolete)', '?', '+-']
    for filter in filters:
        network = filter_network(network, filter)

    if 'network_tf_gene7.0' in fnames  and '(obsolete)' not in filters:
        network = filter_network(network, '(obsolete)')

    if 'sRNADataSet' in fnames and 'antisense' not in filters:
        network = filter_network(network, 'antisense')

    # for verbose purposes
    n_comments = sum([network_files[fname]['comments'] for fname in fnames])

    # build the network matrix/genotype; nodes are tf/gene names
    genotype, nodes, n, r, existing = build_network(
        network, fnames, n_comments, verbose)

    # silence a gene
    # ATTENTION: nodes is list!
    for tf in knockouts:
        try:
            knockout = nodes.index(tf)
        except: #IndexError ?
            knockout = nodes.index(tf_to_gene(tf))
        genotype = gene_knockout(genotype, knockout)
        del nodes[knockout]

    if verbose:
        verbose_network(
            short_name, genotype, network, nodes, n, r, knockouts, existing)

    if random: genotype = randomize_genotype(genotype)

    return genotype, network, nodes

def gene_knockout(genotype, knockout):
    genotype = delete(genotype, knockout, axis=0)
    return     delete(genotype, knockout, axis=1)

def verbose_network(
        short_name, genotype, network, nodes, n, r, knockouts, existing):

    # get number of genes coding for TFs
    if 'tf_tf' in short_name:
        n_tfs = len(genotype)
    else:
        n_tfs = nonzero([node in tf_gene_map['gene_to_tf'].keys()
                         for node in nodes])[0].size

    # get number of genes that are not regulated
    not_regulated = where(get_in_degrees(genotype) == 0)[0].size

    # get genes that are only self-regulated and regulate at least another gene
    candidates = get_candidates_list(genotype, network, nodes)
    candidates.sort(key=lambda t: t[2], reverse=True)
    #candidates = array(sorted(candidates.items(), key=lambda t: t[1])[::-1])

    # filename
    fname  = short_name + '_N%d' %len(nodes)
    suffix = '_%s' %knockouts if knockouts else ''

    # save list of candidates, translating genes to TFs when possible
    candidates = array(candidates)
    if candidates.size:
        candidates[:,0] = parse_list(candidates[:,0], gene_to_tf)
        file = nets_dir + fname + '_candidates' + suffix + '.txt'
        savetxt(file, candidates, '%s')

    # save tf or gene names with out degree
    nodes = parse_list(nodes, gene_to_tf)
    out_degrees = get_out_degrees(genotype)
    file = nets_dir + fname + '_names' + suffix + '.txt'
    savetxt(file, array((nodes, out_degrees)).T, '%s')

    save_network(fname, existing)

    # write some statistics to a log file
    print_some_network_statistics(
        short_name, genotype, network, n_tfs, not_regulated, candidates, n, r,
        knockouts)

def save_network(fname, network):
    file = open(nets_dir + fname + '.txt',  'w')
    new = []
    for connection, link in network.iteritems():
        source, target = parse_list(connection.rsplit('-', 1), gene_to_tf)
        new.append((source, target, link.strip('?')))
    savetxt(file, new, '%s')
    file.close()

def build_network(
        network, fnames, n_comments=0, verbose=False, base=def_base):

    sources  = [edge[0] for edge in network]
    targets  = [edge[1] for edge in network]
    nodes    = list(unique(sources + targets))
    size     = len(nodes)
    genotype = zeros((size, size))
    del sources, targets

    # for verbose purposes
    sources = []
    targets = []
    n = {'-':0, '+':0, '+-':0, '-?':0, '+?':0, '?': 0, 'errors':0, 'repeats':0,
         'lines':len(network)} # total
    r = {'-':0, '+':0, '+-':0, '-?':0, '+?':0, '?': 0} # self-regulation
    existing = {}
    if verbose:
        file = open(nets_dir + regulon_err_file,  'a')
        print >> file, network_files[fnames[0]]['release']
        print >> file, "printing to log (from %s)" %fnames

    for i, edge in enumerate(network):
        source = nodes.index(edge[0])
        target = nodes.index(edge[1])
        link   = edge[2]

        sources.append(edge[0])
        #targets.append(edge[1])

        # repeated interactions in the data file
        if verbose:
            connection = '-'.join([edge[0], edge[1]])
            if connection in existing.keys():
                n['repeats'] += 1
                if existing[connection] != link and link != '(obsolete)':
                    print >> file, connection,
                    (' link already exists with different sign between lines' +
                     '%s and %s!' %(
                         n_comments + i, n_comments + i + n['repeats']))
            else:
                existing[connection] = link

        # negative
        if   link == '-':
            genotype[source, target] = base[0]
            n['-']  += 1
            if source == target:  r['-'] += 1
        elif link == '-?':
            genotype[source, target] = base[0]
            n['-?'] += 1
            if source == target: r['-?'] += 1

        # positive
        elif link == '+':
            genotype[source, target] = base[1]
            n['+']  += 1
            if source == target:  r['+'] += 1
        elif link == '+?':
            genotype[source, target] = base[1]
            n['+?'] += 1
            if source == target: r['+?'] += 1

        # random
        elif link == '+-':
            genotype[source, target] = base[random_integers(0,1)]
            n['+-'] += 1
            if source == target: r['+-'] += 1
        elif link == '?':
            genotype[source, target] = base[random_integers(0,1)]
            n['?']  += 1
            if source == target:  r['?'] += 1

        # not attributed
        else:
            n['errors'] += 1
            sources.pop()
            #targets.pop()
            if verbose: print >> file, n_comments + i, link

    # verbose
    #nodes = unique(sources + targets)
    n['regulators']   = unique(sources).size
    if verbose: file.close()

    return genotype.T, nodes, n, r, existing

def build_all_networks(
        short_names=net_short_names.keys(), filters=[[], ['?', '+-']],
        verbose=True):

    for short_name in short_names:
        for filter in filters:
            genotype, network, nodes = get_network(
                short_name, filter, verbose = verbose)
    #return genotype.T, network, nodes

# get genes that are only self-regulated and regulate at least another gene
def get_candidates_i(genotype):
    candidates = []
    for i in range(len(genotype)):
        row = nonzero(genotype[i,:])[0]
        # regulated by only one gene
        if row.size == 1:
            # only self-regulation
            if row == [i]:
                col = nonzero(genotype[:,i])[0]
                # regulates at least one more gene other than itself
                if col.size > 1:
                    candidates.append(i)
    return candidates

# get their out-degree and sign of self-regulation as well
def get_candidates_list(genotype, network, nodes):
    self_reg     = where(network[:,0] == network[:,1])[0]
    #print network, self_reg
    self_reg     = network[self_reg]
    candidates   = []
    candidates_i = get_candidates_i(genotype)
    for gene in candidates_i:
        degree = nonzero(genotype.T[gene])[0].size
        gene   = nodes[gene]
        sgn    = self_reg[where(self_reg == gene)[0][0]][2]
        candidates.append((gene, sgn, degree))
    return candidates
