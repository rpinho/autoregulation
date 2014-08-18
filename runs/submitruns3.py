#!/usr/bin/env python
#!/usr/bin/python

#import subprocess
import sys
import os
import time
# barley
if sys.platform == 'linux2' and '+' in sys.version:
    import submit_barley as submit
# evo or mac
else:
    import submit
sys.path.append(os.getcwd())
import runs.run3
reload(runs.run3)
from runs.run3 import *
script = '/runs/run3.py'

# networks
bits = [4]#range(4,20)
#list(get_int_datapoints(log10(3000), log10(10000), 4))[1:]
densities = [c]#c#None
degrees = [None]#None#'N-1'#2
act_fractions = [p]#[None] #[0., .1]#arange(0, 1.1, .1)
#ps = [p]#[.5407885]#, .53403611]#[.9, 1.]
#ps = ['ps_all']#arange(.2, 1.1, .1)
#qs = [.63]#'qs_all']#arange(91)/90.#arange(0, 1.01, .01)
#delete(qs, range(0,101,10))
#   = delete(qs, [5, 16, 25, 36, 45, 56, 65, 76, 85, 96])
# repeated int elements of arange(0, 1.01, .01)
#   = ps[::-1]
#   = qs[::-1]

# development process parameters
#phe_decimal   = 1023
mins = [0.]#, -1.]
dims = [2]#range(3, 10)#k_dict(bits[0]) + 1)#[2]#[3]
noises = [1e-100, .01, .05, .1, .15, .2, .25]#, .3, .35, .4, .45]
#noise_function = 'random_noise'
noise_function = 'random_flip'#force_different_flip'#neighbour_flip'#
noise_time = 'before'#shmulevich'#after'#
if 'flip' in noise_function and noise_time != 'shmulevich':
    noises = [1e-100, 1, 2, 3, 4]
_devo_times = [inf]#'maximum'#'mean'#inf ## NOTE: inf is not a String!
if any(noises):
    devo_time = 'maximum'
    if dims == devo_times[devo_time].keys():
        _devo_times = [devo_time]#'mean'
    else:
        _devo_times = [100]#[devo_times[devo_time]['2']['1.0']['%d' %bits[0]]]
gpmap = None#'csign'#cHeaviside0#
steps = 'uniform'#'quantiles'#
steepness = a
converge_time = T
converge_threshold = converge
robustness_threshold = .1

# genotype/initials sampling
deterministic = True#False
# NOTE: i don't think random is used in get_individuals()/genotypes
random = True#False
stable = False#True
#if stable and any(noises) and _devo_times == [100]:
 #   _devo_times = [200]
sel_period = 1
max_trials = F*10 # for generate_stable_genotype()
binary = True#False
full_enum = False#True
if full_enum:
    binary = True
symm = False#True
shuffling = False # does NOT increase the search speed
genotype_func = 'generate_genotype'
#'diagonal_p_genotype'#'diagonal_0_genotype'#'generate_whole_genotype'
ini_random = False#True
get_initials = True#False
## ATTENTION: if myrun.experiment == 'distribution': myrun.get_initials = False
initials = None
#array([0, 2, 6, 8, 18, 20, 24, 26, 54, 56, 60, 62, 72, 74, 78, 80])#None
## ATTENTION: if not get_initials and not initials:
# initials = phenotypic_space_function and n_samples *= initials.size
graph = False#True
scale_free = False#True

# mutants
n_perturb = 100
change_sign = False#True
deletion = False#True
n_mutants = None
all_mutants = None#1
n_neighbours = 0 # degree of neighborhood of perturbed initials

# evolution
models = [4] # 1, 2, 3 or 4
periods = [1]#, 2, 3, 4, 5]
generations = 1e6#G
sel_strengths = [s]#*10, s/10]
mut_rates = [u]#*10, u/10]
mut_biases = [None]
rec_rates = [r]#, 0]
selection = False # save pop before (False) or after (True) selection

# autoregulation
pop_size = P#1e3
pop_stat_func = 'get_ind_stat'#'get_pop_stat'
pop_stats = ['2p']#, '1-p']#sign']#np', 'nq']#p']#, 's', '1-p']#, 'q', 'P', 'c']
#'p', '1-p', 'qp', 's' + ['l', 'f']
#pop_stats     = ['robustness']#, 'survivability']#diff']#initials_all
#pop_stats     = ['conservation']
# 'equal' stats have to be separated because of:
# pop_function = p.pop_robustness if 'equal' in self.pop_stats else p.pop_evolution
#pop_stats     = ['equal', 'equal_fix', 'intersect']
#pop_stats     = ['phe_decimal']#notype']#diversity']
#q+q']#n+clusters', 'n-clusters']#'q+clusters', 'q-clusters',
#'a+clusters', 'a-clusters']
#g1, g2 = 1e5, 4.9e5

# RegulonDB's E. coli TRN
network       = None#'meta3'#'tf_gene'#None
filters       = [['?', '+-']]#, []]
knockouts     = []#'ArgR']#'CRP'#'Fis']#'CpxR#'LexA
fixed_states  = None#{'LexA':1}
#random        = True

# nx.DiGraph(nx.scale_free_graph) paramaters
delta         = .1
end           = .3
alphas        = arange(delta, end, delta)
delta         = .01
gammas        = arange(delta, .21, delta)
delta         = .1
deltas_in     = arange(.7, 1+delta, delta)[:-1]
deltas_out    = array([.2])#arange(delta, 1+delta, delta)

# new Evolution
shortname = 'm3r05b'#None#'m3r05p7'
if shortname:
    (models, periods, generations, sel_strengths, mut_rates, rec_rates, bits,
     densities, degrees, n_runs, precision, _filter, suffix, binary, graph,
     rec_function, ps, qs, mut_biases, mins, _devo_times) = map(
         atleast_1d, pop_model(shortname)[:-2])
    gpmap, steepness = pop_model(shortname)[-2:]
# output
i = 0; j = 0; n_nodes = 16#8
start = 1#j + n_nodes*i + 1
n_runs = 1#j + n_nodes*(i+1)
runs = arange(start, n_runs + 1)
delete_runs = []#[51, 72, 85, 87, 90]
runs = delete(runs, [where(runs == run)[0] for run in delete_runs])
if not shortname:
    samples = 1e8#0#65536#None
else:
    samples = None
#generations = 1e6 + 10
#mut_rates = rec_rates = [0]
#suffix = ['-u0g%d' %(generations-1e6)]
n_datapoints = 100 #'all1000'
if shortname:
    datapoints = get_pop_datapoints(
        shortname, generations, _filter, n_datapoints)[:samples]
else:
    datapoints = None #array([1e6])
experiment = 'pop.stats'
#phenotypes'#-period_distribution'#stability'#'#path_length'#
#pop.evolution'#robustness_and_survivability'#qp.
#initials'#stable.pop.evaluate'#final_vs_p'#.u0'#stable.pop.superfrankenmatrix'
#survivability'#robustness'#autoregulation'#epistasis'#perturb_and_analyse'
#'stability'#'scale_free_graph'#'distribution'
#'viability'
extension = '.dat'
reshape = True#False
normalize = False#True # normalize entropy by -log(unique_phenotypes.size)
sampling = False#True
if 'stability' in experiment or 'phenotypes' in experiment:
    sampling = True # because there's no save
only_unique = False#True
_get_orbit = False#True
verbose = True#False
rm_file = False#True

# run/cluster
submitrun = True#False # submit or just create and save myrun file? (dry-run)
setup = False # pre-setup
if not submitrun:
    setup = True
save = True # pickle the Run instance to a file after setup() in run.py
interactive = False
job = ''
jobname = None
#qname = 'high'#'low'#'med'#'main.q'#'long.q'#'precise.q'#'precise-long.q'

# barley
#if sys.platform == 'linux2' and '+' in sys.version:#sys.hexversion == 34013680:
#    if qname == 'low': qname = 'long.q'
#    if qname == 'med': qname = 'main.q'
#    qname = 'main.q'

# default behaviors
if binary:
    if not any(mut_biases):
        change_sign   = True

    if 'generate' in genotype_func:
        genotype_func = 'spin_genotype'

if len(filters) == 1:
    filter = filters[0]

# experiment
if 'epistasis' in experiment:
    n_perturb = 16

if ('autoregulation' in experiment or
    'stable.pop' in experiment or
    'qp.pop' in experiment):

    if 'qp.pop' in experiment:
        genotype_func = 'qp_genotype'
    else:
        genotype_func = 'diagonal_p_genotype'

    if 'evolution' not in experiment:
        generations = None
        pop_size    = 1e3# if not random else 1e5
        ini_random  = True
        periods     = atleast_1d(sel_period)

if 'stats' in experiment:
    runs = atleast_1d(n_runs)
    if 'qp' in pop_stats and 'm3r05' in shortname and n_datapoints == 'all1000':
        runs = datapoints
    if 'robustness' in pop_stats:
        pop_stat_func = 'get_pop_stat'
    #datapoints = [g1, g2]

elif shortname:
    extension = suffix[0] + extension

if 'final_vs_p' in experiment:
    genotype_func = 'qp_genotype'
    ps            = [None]
    pop_size      = None

if genotype_func != 'qp_genotype':
    if not shortname:
        qs = [None] # for get_suffix() in myrun
    if genotype_func != 'diagonal_p_genotype':
        ps = [None] # for get_suffix() in myrun

# run jobs interactively/locally on corn or evo
if (len(sys.argv) == 2 and
    (sys.argv[1] == 'corn' or sys.argv[1] == 'evo' or sys.argv[1] == 'nohup')):
    interactive = True

# mac and mpinho-laptop2
elif sys.platform != 'darwin' and sys.hexversion != 33949168:
    # parse command line
    if (len(sys.argv) < 2):
        print "usage: submitruns3.py <dir> <qname>"
        sys.exit()
    # dir is the name of a particular set of runs
    dir = sys.argv[1]
    qname = sys.argv[2] if (len(sys.argv) == 3) else ""
    # e.g. '-l hostname=evo14' or '-l mem_free=1G'
    qsub_extras = sys.argv[3] if (len(sys.argv) == 4) else ""

    # This makes sure <workdir>/<dir> exists
    # All output files go into <dir>/*
    # stdout and stderr for into <dir>/outs and <dir>/errs, respectively
    submit.ensure_dir_exists(dir)
    submit.ensure_dir_exists(dir+"/outs")
    submit.ensure_dir_exists(dir+"/errs")

    # do everything in <dir>
    os.chdir(dir)

# Maybe you want to put your own loops for parameter ranges, in
# this example we just have one loop 'idx'.
for run in runs:
    x = itertools.product(
        bits, densities, degrees, mins, dims, _devo_times, noises,
        act_fractions, ps, qs)
    for bit, density, degree, min, dim, devo_time, noise, act_fraction, p, q in x:
        #y = itertools.product(alphas, gammas, deltas_in, deltas_out)
        #for alpha, gamma, delta_in, delta_out in y:
        z = itertools.product(
            models, periods, sel_strengths, mut_rates, rec_rates, mut_biases)
        for n_model, sel_period, sel_strength, mut_rate, rec_rate, mut_bias in z:
        #for filter in filters:

            if degree == 'N-1':
                degree = bit-1

            # 'prefix' is the prefix of the output files
            prefix = 'N%d' %bit
            if density and density < c:
                prefix += 'c%g' %density
            elif degree and degree < bits:
                prefix += 'k%s' %degree
            if not min:
                prefix += 'm%d' %min
            if dim > def_base.size:
                prefix += 'k%s' %dim
            if not full_enum:
                if binary:
                    prefix += 'b'
                else:
                    prefix += 'r'
                if stable:
                    prefix += 's'
                else:
                    prefix += 'u'
            if noise:
                prefix += get_noise_run_prefix(
                    noise_function, noise, noise_time)
            if samples and not full_enum:
                prefix += 'n%d' %log10(samples)
            if steps == 'quantiles':
                prefix += 'q'

            if experiment == 'viability':
                if deletion:
                    prefix += 'u0'
                elif change_sign:
                    prefix += 'u-'
                else:
                    prefix += 'uu'
                prefix += 'all' if all_mutants else '%d' %n_mutants

            elif 'pop.' in experiment and generations:
                #prefix = get_evolution_shortname(
                 #   n_model, generations, sel_period, sel_strength, mut_rate,
                  #  rec_rate, run, g, precision)
                prefix = get_evolution_shortname(
                    n_model, None, None, sel_strength, mut_rate, rec_rate,
                    shortname, run)
                if 'stats' in experiment:
                    prefix = prefix.replace(print_run(run)[1:], pop_stats[0])
                elif 'robustness' in experiment:
                    prefix += 'r'
                elif len(prefix) > 10:
                    prefix = prefix.translate(None,'n')

            elif experiment is 'scale_free_graph':
                prefix = ('N%sa%.2fg%.2fi%.2fo%.2f' %
                          (bit, alpha, gamma, delta_in, delta_out))

            elif network:
                prefix = network

            elif experiment == 'autoregulation' or 'stable.pop' in experiment:
                prefix = 'p%gn%dn%d' %(act_fraction, log10(pop_size), run)
                if 'superfrankenmatrix' in experiment or p is not None:
                    prefix += 'p%g' %p

            elif 'qp.pop' in experiment:
                prefix = 'S' if stable else 'R'
                prefix += 'q%.2g' %q if type(q) != str else q

                if 'stats' not in experiment:
                    prefix += 'p%.2g' %p if type(p) != str else p
                    prefix += 'n%d' %run

                if sel_period  > 1:
                    prefix += 'p%s' %sel_period

                if len(prefix) > 12:
                    prefix = prefix.translate(None,'SRpn')

                if len(prefix) > 11:
                    prefix = prefix.translate(None,'SRp')

                if len(prefix) > 10:
                    prefix = prefix.translate(None,'SR')

            elif experiment == 'final_vs_p':
                prefix = 'n%d' %run

            else:
                # 'maximum' or 'mean'
                if isinstance(devo_time, str):
                    prefix += '%s' %str(devo_time)[:4]
                # inf, 100, etc. NOTE: type(inf) == float
                else:
                    prefix += '%g' %devo_time

                if experiment == 'stability':
                    prefix += 's'

                #elif experiment == 'phenotypes':
                 #   prefix += 'p'

            if not jobname:
                jobname = job + prefix

            # here is where you set up the parameters for a particular run
            # in a Run instance
            myrun               = Run()
            myrun.bits          = bit
            myrun.density       = density
            myrun.degree        = degree
            myrun.samples       = samples
            myrun.stable        = stable
            myrun.binary        = binary
            myrun.full_enum     = full_enum
            myrun.symm          = symm
            myrun.random        = random
            myrun.initials      = initials
            #myrun.phe_decimal   = phe_decimal
            myrun.get_initials  = get_initials
            myrun.min           = min
            myrun.dim           = dim
            myrun.noise         = noise
            myrun.noise_function = noise_function
            myrun.noise_time    = noise_time
            myrun.run           = run
            myrun.experiment    = experiment
            myrun.sampling      = sampling
            myrun.n_datapoints  = n_datapoints
            myrun.network       = network
            myrun.only_unique   = only_unique
            myrun.get_orbit     = _get_orbit
            # files
            myrun.prefix        = prefix
            myrun.extension     = extension
            myrun.verbose       = verbose
            myrun.rm_file       = rm_file
            myrun.precision     = precision
            myrun.reshape       = reshape
            myrun.dir           = get_experiment_dir(myrun.experiment)
            myrun.save          = save
            myrun.submitrun     = submitrun
            # genotype
            myrun.filters       = filter
            myrun.knockouts     = knockouts
            myrun.graph         = graph
            myrun.scale_free    = scale_free
            myrun.genotype_func = genotype_func
            myrun.act_fraction  = act_fraction
            myrun.p             = p
            myrun.q             = q
            myrun.deterministic = deterministic
            myrun.max_trials    = max_trials
            myrun.shuffling     = shuffling
            # devo
            myrun.devo_time     = devo_time
            myrun.gpmap         = gpmap
            myrun.steepness     = steepness
            myrun.converge_time = converge_time
            myrun.converge_threshold   = converge_threshold
            myrun.robustness_threshold = robustness_threshold
            myrun.steps         = steps
            myrun.ini_random    = ini_random
            myrun.pop_size      = int(pop_size)
            myrun.pop_stats     = pop_stats
            myrun.pop_stat_func = pop_stat_func
            myrun.pop_start     = start
            myrun.sel_period    = sel_period
            # mutants
            myrun.mut_bias      = mut_bias
            myrun.change_sign   = change_sign
            myrun.deletion      = deletion
            myrun.n_perturb     = n_perturb
            myrun.n_mutants     = n_mutants
            myrun.all_mutants   = all_mutants
            myrun.n_neighbours  = n_neighbours
            myrun.shortname     = shortname

            # nx.DiGraph(nx.scale_free_graph) paramaters
            if experiment is 'scale_free_graph':
                myrun.alpha = alpha
                myrun.gamma = gamma
                myrun.delta_in = delta_in
                myrun.delta_out = delta_out

            # e.coli network
            if network:
                test_genotype, net, nodes = get_network(network, filter)
                myrun.bits = len(test_genotype) - len(knockouts)
                myrun.degree = get_in_degrees(test_genotype).mean()
                myrun.density = myrun.n_zeros = None
                myrun.binary = True
                myrun.prefix += 'N%dm%d' %(myrun.bits, min)
                del test_genotype

                # fixing expression state levels of some genes
                if fixed_states:
                    myrun.fixed_states = fixed_states
                    myrun.devo_function = 'devo_fixed_state'
                    myrun.nodes = parse_list(nodes, gene_to_tf)
                    myrun.prefix += str(fixed_states)

                if knockouts:
                    myrun.prefix += str(knockouts)

                if random:
                    myrun.prefix += '-random'

            if gpmap:
                myrun.prefix += myrun.gpmap

            # pop.evolution experiments
            if 'pop.' in experiment and generations:
                myrun.model        = n_model
                model              = get_model(myrun.model, generations)
                myrun.sel_strength = sel_strength if n_model == 4 else model['s']
                myrun.min_fit      = model['min_fit']
                myrun.generations  = generations#model['G']
                myrun.sel_period   = sel_period
                myrun.mut_rate     = mut_rate
                myrun.rec_rate     = rec_rate
                if rec_function: myrun.rec_function = rec_function
                myrun.selection    = selection
                myrun.datapoints   = datapoints
                myrun.dir          = get_pop_dir(shortname, run, myrun.dir)
                if 'stats' not in experiment:
                    submit.ensure_dir_exists(myrun.dir)

                if experiment == 'perturb_and_analyse':
                    myrun.normalize    = normalize

                if generations:
                    myrun.suffix = '-' + get_evolution_filename(
                        myrun.model, myrun.generations, myrun.sel_period,
                        sel_strength=myrun.sel_strength,
                        mut_rate=myrun.mut_rate,
                        rec_rate=myrun.rec_rate, run=myrun.run,
                        mut_suffix = get_mut_suffix(
                            myrun.deletion, myrun.change_sign, myrun.mut_bias),
                        rec_function = rec_function)

            # set it up before saving
            if setup:
                myrun.setup()

            # pickle the Run instance to a file
            fname = myrun.prefix + ".run"
            file  = open(fname,  "w")
            cPickle.dump(myrun, file)
            file.close()

            # submit a job that reads the pickled Run into the run.py script
            # to run it
            line = 'time ' + workdir + script + ' < %s.run' % myrun.prefix

            # mac and mpinho-laptop2
            if (sys.platform != 'darwin' and sys.hexversion != 33949168 and
                not interactive):
                submit.qsub_launch(
                    dir, myrun.prefix, line, qname, jobname, qsub_extras)
                print line

            elif submitrun:
                if sys.platform != 'darwin' and sys.argv[1] == 'nohup':
                    os.system(('nohup ' + line + '> %s.out 2> %s.err &' %
                               (myrun.prefix, myrun.prefix)))
                else:
                    os.system(line)
                #subprocess.Popen('./runs/run.py < N4k2n3infcsign.run')
                #subprocess.call(['./runs/run.py', '< N4k2n3infcsign.run'])
                print line

            else:
                print ('Job not submitted. Pickled the run to file as %s.run'%
                       myrun.prefix)

            jobname = None
