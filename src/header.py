import sys
#from pprint import pprint as pp; pp(sys.path)
# evo
if sys.platform == 'linux2' and sys.hexversion == 34013936:
    sys.path.append("/home/mmclaren/pyenv/lib/python2.7/site-packages")

#from pylab import *
from numpy import *
from numpy.random import *
from scipy import stats, weave, misc, ndimage, optimize
from scipy.stats import rv_discrete
#import scipy.weave as weave
#import networkx as nx

import os
import operator
import csv
import cPickle
import copy
import itertools
from types import *
from random import choice
#import psyco

# dirs
workdir       = ''.join(os.getcwd().partition('model')[:2])
maps_dir      = workdir + '/maps/'
data_dir      = workdir + '/data/'
ensembles_dir = workdir + '/ensembles/'
analysis_dir  = workdir + '/analysis/'
pop_dir       = workdir + '/populations/'
figs_dir      = workdir + '/figs/'
nets_dir      = workdir + '/networks/'
pa_dir        = workdir + '/perturbations/' # perturb_and_analyse
epi_dir       = workdir + '/epistasis/'
surv_dir      = workdir + '/survival/'
paper_dir     = 'figs/paper/'
report_dir    = 'report/'
conf_dir      = 'conferences/'
logs_dir =  workdir + '/logs_and_runs/'

# mut_bias
mut_biases = [2.33, 1.64, 1.28, 1.037, .84, .674, .583, .52, .385, .25,
              .228, .2, .176, .151, .126, .1, .076, .05, .025, 0,
              -.025, -.126, -.25, -.385, -.52, -.674, -.84, -1.037, -1.28,
              -1.64, 2.33]
qp_biases = [.99, .95, .9, .85, .8, .75, .72, .7, .65, .6, .59, .58, .57,
             .56, .55, .54, .53, .52, .51, .5, .49, .45, .4, .35, .3, .25,
             .2, .15, .1, .05, .01]
#append([.05], arange(.1, 1, .1))
mut_biases_dict = dict(zip(qp_biases, mut_biases))
# std of the normal distribution if not deterministic in get_pos_neg_zeros()
qp_scale = .145

## PARAMETERS 0
n_states = 2
min_state = -1.
max_state =  1.
# already built maps. Note: There are more!
bases = ['gray', 'binary', 'decimal', '[-1.  1.]']
# default base map. ATTENTION: model.py may change this!
# Note: '[-1.  1.]' is the FASTEST!
def_basename = bases[-1]
# default normalization function. ATTENTION: model.py may change this!
# Note: sign is the FASTEST!
def_gpmap = sign
# for putting an upper limit on full enumeration of (binary) genotype space
enum_bits = 5
# for putting a lower limit on the use of the sign function to correct for
# sign(0) = 0
sign_bits = 5
# saved dictionary/maps are N = 2, 4 <= N <= 11.
# For bigger N use decimal_to_vector and vector_to_decimal (strings and ints)
N_dict = 11
#
Dtype = float64

## PARAMETERS 1

# A. Wagner: 4 < N < 10. Ciliberti: 5 < N < 40. GP Wagner: 4
#N = [4, 6, 8, 10, 10, 10, 11, 10, 10, 4, 12, 20, 40, 4, 7, 10,
 #    8, 25, 10, 10, 100, 10, 20, 4, 12, 20, 10, 10, 20, 50]
# in the literature
N = 10

# graph or network density (aka connectivity coeficient)
# (fraction of links): http://en.wikipedia.org/wiki/Graph_density
# A. Wagner: 0.4 < c < 1. Siegal: 0.144 < c < 0.75. Ciliberti: 0.1 < c < 0.5
#c = [ 0.4, 0.6, 0.8, 1, 0.144, 0.4, 0.75, 0.75, 0.26, 0.5, 0.75, 0.25, 0.5,
 #     0.25, 0.5, 0.08, 0.25, 0.75, 1, 0.25] in the literature
c = 1.

# population size (Wagner and Siegal: 500. GP Wagner: 1000)
P = 500

# recombination (0 == no recombination)
r = .5

# mutation rate (per genome)
# Wagner: 1. Siegal: 0.1. Mike: 0.1, 0.01, 0.003. GP Wagner: 0.00317 per edge
# (N=4, c(initial)=.5) i.e. u ~ 0.025 (initial)
u = .1

# insertion rate (per genome), i.e., the number of connections (edges) created
# GP Wagner: 2*uA = 2*uD = 3.17e-5 per edge (N=4) i.e. uA ~ 2.5e-4
uA = 0 #u/2

# deletion rate (per genome), i.e., the number of connections (edges) erased
# GP Wagner: 2*uA = 2*uD = 3.17e-5 per edge (N=4) i.e. uD ~ 2.5e-4
uD = 0 #uA

# selection strength
# Wagner: 0.1. Siegal: 0.1 strong, 1 intermediate, inf no stabilizing selection
s = .1

# the fitness of individuals with period (stability) NOT selected
# (period > 1 == cycles)
# Wagner = exp(-1/s) ~ 0, Siegal = 0, random = 1
zero_fit = 0

# number of generations (Wagner and Siegal: 400)
G = 400

# Wagner: 150 to 300 sample (ensemble) size
n_runs = 50

# number of perturbations (mutations or noise) to the ensemble members when
# evaluating genome stability/robustness
# mutations: Wagner: 5000/~P ~ 10
# 10 < n_noise < 1000
n_perturb = 100

# 1 or 2% of the population, i.e., low frequency
n_invaders = int(.02 * P)

# making sure it's even not odd (because of recombination)
n_invaders = n_invaders if not bool(n_invaders % 2) else n_invaders + 1

## PARAMETERS 2
# steepness of the sigmoid function (a = 10 steeper)
## NOTE: I wonder if I can compare this to n_states
# A. Wagner: uses sign. Siegal & Bergman a = 100. Mike a = 1, 10
a = 100

# maximum development time (to reach equilibrium):
# depends on N - with noise you have big path lengths! >100
# A. Wagner: 3N. Usually takes < 20 or it doesn't converge (N=10)
# Siegal and Mike: 100
# with brent algorithm you define it in powers of 2:
# 64 is a good choice (no noise)
M = 100

# maximum trials for the founder genotype to converge
# A. Wagner: 5*pow(2,2*N)
F = 10000

# average over last T phenotypes: T = 2
# Siegal and Mike: T = 10
T = 10

# fraction of 1s in the integer-valued (spin) phenotype
# A. Wagner: 0.1 < p < 0.9
# ATTENTION: I'm using it as activating fraction, like McDonald et al.,
# as opposed to Wagner's repressing fraction (?)
p = 0.5

# stability convergence threshold value (epsilon)
## NOTE: I wonder if I can compare this to the width of my steps
converge = 1e-4

ps_all = arange(11)/10.
qs_all = arange(91)/90.

# saved dictionary/maps are 2 <= k <= 4 for N=10
# For bigger k use decimal_to_vector and vector_to_decimal (strings and ints)
def k_dict(bits):
    if bits <= 4:
        return 9
    elif bits == 5:
        return 7
    elif bits == 6:
        return 6
    elif bits == 10:
        return 3
    else:
        return 2

parameters = '_N%(N)s_c%(c)s_T%(T)s' % locals()

def get_model(number, gen=G):
    model_names = ['random: no evolution',
                   'random evolution: no selection, just sampling',
                   'Siegal: evolution with selection for fixed points',
                   ('Wagner: evolution with selection for fixed points and'
                    + 'target phenotype')]
    model = []
    loop = zip(
        [0, gen, gen, gen], [None, inf, inf, 0.1], [None, 1, 0, 0], model_names)
    for g, s, min_fit, name in loop:
        model.append(dict(G=g, s=s, min_fit=min_fit, name=name))
    return model[number-1]
