#!/usr/bin/env python
#!/usr/bin/python

import sys
import os
workdir = ''.join(os.getcwd().partition('model')[:2])
sys.path.append(workdir)
import src.classes3
reload(src.classes3)
from src.classes3 import *

# filename_exceptions
density_precisions = {300:3, 1000:3, 2000:3, 3000:4, 5000:3, 7000:3, 10000:4}
#density_precisions = {174:3, 300:3, 418:3, 559:3, 747:3, 1000:3, 2000:3,
 #                     3000:4, 5000:3, 7000:3, 10000:4}
# if N > 1000 then precision > 2
precision = None

class Run:
    def __repr__(self):
        return str(self.__dict__)

    def __init__(self):

        # Network parameters
        self.bits = N
        self.density = None#c
        self.n_zeros = None#0
        self.degree = None#self.bits
        # fraction of activating (+1) interactions in the network
        self.act_fraction = p

        # State space parameters
        self.dim = 2
        self.min = -1.

        # Development process parameters
        self.devo_time = inf
        self.noise = 0
        self.gpmap = None
        self.steepness = a
        self.converge_time = T
        self.converge_threshold = converge
        self.robustness_threshold = .1
        self.steps = 'uniform' #'quantiles'

        # Network sampling parameters
        self.samples = 1e4
        # full enumeration of genotype space (binary matrices)
        self.full_enum = False
        # with or without symmetry constraints
        self.symm = False
        # conditional sampling of genotype space (real matrices)
        self.stable = False
        # random binary matrices
        self.binary = False
        # samples = number of individuals (pairs matrix x initial) if get_initials
        # else just number of matrices, and initials = phenotypic_space_function
        self.get_initials = True
        self.initials = None
        # file with saved genotypes
        self.genotypes = None
        self.network = None
        self.knockouts = []
        self.fixed_states = None

        # Output parameters
        # ATTENTION: for decimals use any True value if you want to filter.
        self.filter = None
        # either just count # fixpoints (stability) or actually store the list
        # of fixpoint phenotypes (unique(phenotypes) = #stable states)
        self.experiment = 'stability'
        # save just unique(phenotypes) or all phenotypes
        self.only_unique = True
        # save (i, n) to file every data point
        self.sampling = True
        # maximum number of data points (i, n) to be stored for each set of
        # parameters (min, dim, bits)
        self.n_datapoints = 20
        self.datapoints = None
        self.filename = None
        self.temp_file = None
        self.dir = data_dir
        self.suffix = ''
        self.extension = '.dat'
        self.verbose = False

        # devo function
        self.devo_function = 'devo'

        # Mutations
        self.n_mutants = None
        self.all_mutants = None
        self.deletion = False
        self.change_sign = False
        self.mut_bias = False

        # Evolution
        self.mut_rate = u
        self.rec_rate = r
        self.rec_function = 'recombine'


    def get_prefix(self):

        if 'pop.' in self.experiment and 'evolution' not in self.experiment:
            prefix = self.experiment.split('.')[-1]#.strip('pop.')
        else:
            prefix = self.experiment

        if 'stable.pop' in self.experiment:
            prefix = 'stable.pop'

        if 'qp.pop.initials' in self.experiment:
            prefix = 'qp.pop'

        return prefix


    def get_filename(self, precision):
        return ''.join((print_N(self.bits),
                        print_c(self.density, precision),
                        print_min(self.min),
                        print_dim(self.dim),
                        print_noise(
                            self.noise_function, self.noise, self.noise_time)))


    def get_suffix(self, gpmap, devo_time):
        return get_suffix(
            self.samples, self.full_enum, self.symm, self.stable, self.binary,
            self.filter, self.fix_gpmap, gpmap, self.steepness, self.steps,
            self.converge_time, self.converge_threshold, devo_time,
            self.act_fraction, self.p, self.q, self.sel_period,
            self.n_mutants, self.all_mutants, self.deletion, self.change_sign,
            self.mut_bias, self.network, self.knockouts, self.fixed_states,
            self.random, self.get_initials, self.graph, self.scale_free,
            self.genotype_func, self.pop_size, self.run, self.experiment,
            self.deterministic)


    def get_individuals(
            self, genotype_function, genotypic_space_function,
            phenotypic_space_function):

        return get_individuals(
            self.full_enum, self.stable, self.samples, self.density,
            self.basename, self.base, self.bits,
            genotype_function, genotypic_space_function, self.get_initials,
            phenotypic_space_function,
            self.noise, self.p, self.q, self.act_fraction, self.genotypes,
            self.network, self.filters,
            self.knockouts, self.random,
            self.initials, self.binary, self.graph, self.scale_free,
            #self.sel_period, self.max_trials, self.deterministic,
            shuffling=self.shuffling)

    def setup(self, gpmap=None, precision=2):

        ## Default behaviors
        if self.full_enum:
            self.binary = True

        if self.experiment == 'distribution':
            self.get_initials = False
        # This guarantees a full enumeration of phenotype space in
        # model.get_individuals()

        ## Set devo time (depends on N)
        key = self.degree if self.degree else self.density
        devo_time = self.devo_time # save str for suffix in filename
        if (type(self.devo_time) is str and
            str(self.dim) in devo_times[self.devo_time].keys()):
            self.devo_time = (devo_times[self.devo_time]['%d' %self.dim]
                              ['%s' %key]['%d' %self.bits])

        ## Get density, number of zeros and degree
        self.density, self.n_zeros, self.degree = get_density_degree_nzeros(
            self.density, self.n_zeros, self.degree, self.bits)

        ## Get base
        self.base, self.basename, jumps = get_base(
            self.dim, self.bits, self.min, self.steps)

        ## Get gpmap
        if not gpmap and self.gpmap:
            gpmap = globals()[self.gpmap]

        if gpmap:
            # constant gpmap independent of base
            self.fix_gpmap = True

            # continous function (with zeros)
            if gpmap.__name__ == 'cSigmoid':

                # force the use of development and not brent
                self.noise = 1e-100

                # gpmap only takes 1 argument; also sets gpmap.__doc__
                def gpmap(x):
                    """cSigmoid"""
                    return cSigmoid(x, self.steepness)

            elif gpmap.__name__ == 'cSigmoid2':
                self.noise = 1e-100
                def gpmap(x):
                    """cSigmoid2"""
                    return cSigmoid2(x, self.steepness)
        else:
            self.fix_gpmap = False
            gpmap = get_gpmap(
                self.base, self.basename, self.density, self.bits, self.binary,
                jumps)

        # save gpmap name only, not the function itself
        if type(gpmap) is str:
            self.gpmap = gpmap
        else:
            if gpmap.__name__ != 'gpmap':
                self.gpmap = gpmap.__name__
            else:
                self.gpmap = gpmap.__doc__

        self.stepfunction = gpmap.__doc__
        # gpmap = globals()[individual.gpmap] # from classes.py L224

        ## Generating functions

        # generating function for sampling genotypes
        genotype_function = globals()[self.genotype_func]
        #self.real_genotype] if not self.binary else spin_genotype

        # generating function for full enumeration of genotypes
        if self.symm:
            genotypic_space_function = generate_all_symm_genotypes
        else:
            genotypic_space_function = generate_full_genotype_space

        # generating function for full enumeration of phenotypes
        phenotypic_space_function = generate_full_phenotype_space
        #generate_all_symm_phenotypes # if noise return vectors, else decimals

        ## Get individuals/samples

        if ('pop' not in self.experiment and
            self.experiment != 'autoregulation' and
            self.experiment != 'final_vs_p'):
            individuals, initials, n_samples = self.get_individuals(
                genotype_function, genotypic_space_function,
                phenotypic_space_function)

        ## Set filter

        # no need to filter if you have the phenotypes.
        # You can apply the filter afterwards
        if self.experiment is 'phenotypes':
            self.filter = None

        # if not noise, we use decimals and devo, so can only filter state = 0
        if self.filter and not self.noise:
             if 0. in self.base:
                 self.filter = get_decimal(
                     zeros(self.bits), self.basename, self.base)
             else:
                 self.filter = None # decimals

        ## Prepare output files

        # filename
        if not self.filename:

            if precision and self.bits in density_precisions:
                precision = density_precisions[self.bits] # filename_exceptions

            prefix = self.get_prefix()

            filename = self.get_filename(precision)

            suffix = self.get_suffix(gpmap, devo_time)

            self.filename = self.dir + prefix + filename + suffix + self.suffix

        # temporary file
        if self.sampling:

            # save temp file at these equally spaced (in log-scale) datapoints
            self.datapoints = get_int_datapoints(
                0, log10(n_samples), self.n_datapoints)

            # overwrite existing file.
            # Note: path_length_vs_N2 (in model.py) saves with:
            # save_txt3(filename, lengths, 'w') so no need..
            if self.experiment != 'path_length':
                self.temp_file = self.filename + '.tsv'
                if self.submitrun and self.experiment != 'phenotypes':
                    file = open(self.temp_file, 'w')
                    file.close()

        #else: self.datapoints = None

        if ('pop' not in self.experiment and
            self.experiment != 'autoregulation' and
            self.experiment != 'final_vs_p'):
            return individuals, initials, gpmap
        else:
            return 1, 1, gpmap


    def get_population(self):

        if 'autoregulation' in self.experiment:
            return Population(
                size=self.pop_size, random=True, ini_random=self.ini_random,
                binary=self.binary, density=self.density,
                basename=self.basename, gpmap=self.gpmap, base=self.base,
                genotype_function=globals()[self.genotype_func],
                bits=self.bits,
                act_fraction=self.act_fraction, p=self.p, q=self.q)

        elif 'stable.pop' in self.experiment:
            return Population(
                size=self.pop_size,
                binary=self.binary, conditional=True, density=self.density,
                basename=self.basename, gpmap=self.gpmap, base=self.base,
                genotype_function=globals()[self.genotype_func],
                bits=self.bits,
                act_fraction=self.act_fraction, p=self.p, q=self.q)

        elif 'qp' in self.experiment:
            return Population(
                size=self.pop_size, random=self.random,
                ini_random=self.ini_random,
                binary=self.binary, conditional=self.stable,
                density=self.density,
                basename=self.basename, gpmap=self.gpmap, base=self.base,
                genotype_function=globals()[self.genotype_func],
                bits=self.bits, sel_period=self.sel_period,
                act_fraction=self.act_fraction, p=self.p, q=self.q,
                max_trials=self.max_trials, deterministic=self.deterministic)

        # population is random:
        # (different genotypes, not conditional on stability)
        else:
            return Population(
                random=True, binary=self.binary, density=self.density,
                basename=self.basename, gpmap=self.gpmap, base=self.base,
                genotype_function=globals()[self.genotype_func],
                bits=self.bits, deterministic=self.deterministic,
                devo_time=self.devo_time, noise=self.noise,
                steepness=self.steepness, converge_time=self.converge_time,
                converge_threshold=self.converge_threshold)


    ##
    def copy_old_pop_to_new_pop(self, pop):
        # set Evolution parameters
        pop.sel_strength = self.sel_strength
        pop.min_fit = self.min_fit
        pop.sel_period = self.sel_period
        pop.rec_rate = self.rec_rate
        pop.rec_function = self.rec_function
        pop.mut_rate = self.mut_rate
        pop.change_sign = self.change_sign
        pop.deletion = self.deletion
        pop.mut_bias = self.mut_bias
        pop.noise = self.noise
        pop.converge_time = self.converge_time
        pop.converge_threshold = self.converge_threshold
        pop.devo_time = self.devo_time

        # legacy code: copy old class to new class
        if pop.__module__ != 'src.classes3':
            return src.classes3.Population(pop = pop)

        else:
            return pop

    ## The experiment: Random Sampling of both genotype and initial state for
    # estimating stability (= i.e. prob(fixpoints))
    def stability_sample(self, individuals, initials, gpmap):

        devo_function = globals()[self.devo_function]

        sampling = False if self.experiment == 'viability' else self.sampling

        self.stability, phenotypes, stability_dist = stability_sample(
            individuals, initials, self.basename, gpmap, self.base, self.bits,
            self.noise, devo_function, self.converge_time,
            self.converge_threshold,
            self.devo_time, self.filter, self.get_initials, self.experiment,
            sampling, self.datapoints, self.temp_file)

        if self.experiment == 'phenotypes':
            return phenotypes

        elif self.experiment == 'distribution':
            return array(stability_dist)

        else:
            return self.stability


    ##
    def phenotype_sample(self, individuals, initials, gpmap):

        devo_function = globals()[self.devo_function]
        noise_function = globals()[self.noise_function]

        return phenotype_sample(
            individuals, initials, self.basename, gpmap, self.base, self.bits,
            self.noise, devo_function, self.converge_time,
            self.converge_threshold,
            self.devo_time, self.get_initials, self.experiment, self.sampling,
            self.datapoints, self.temp_file, self.fixed_states, self.get_orbit,
            self.only_unique, noise_function, self.noise_time)


    ## Generating mutants and evaluating their stability
    def viability_sample(self, individuals, gpmap):

        n = 0
        samples = 0
        for initial, genotype in individuals:

            mutants = get_mutants(
                initial, genotype, self.n_mutants, self.all_mutants,
                self.change_sign, self.deletion)

            n += self.stability_sample(mutants, None, gpmap)

            del mutants

            # save (i, n) to file every data point
            samples +=1
            if self.sampling and samples in self.datapoints:
                save_txt2(self.temp_file, (samples, n))

        del individuals

        return n


    ## Studying a cut-off for devo time
    def path_length_vs_N(self, individuals, gpmap):

        devo_function = globals()[self.devo_function]

        # save temporary file at equally spaced (in log-scale) datapoints
        self.sampling = True
        self.datapoints = get_int_datapoints(
            0, log10(self.samples), self.n_datapoints)

        # filename
        if self.degree < self.bits:
            key = self.degree
        else:#if self.dim > def_base.size:
            key = self.dim
        self.temp_file = (self.dir + 'path_leng_vs_N%dk%d.txt' %
                          (self.bits, key))

        return path_length_vs_N2(
            individuals, self.basename, gpmap, self.base, self.bits, self.noise,
            devo_function,
            self.converge_time, self.converge_threshold, self.devo_time,
            self.filter, self.sampling, self.datapoints, self.temp_file)


    # samples = 1e+04, n_runs = 20
    def phenotype_vs_p(self):
        ps = [[] for i in range(2**self.bits)] #1024
        for i in xrange(int(self.samples)):
            p = randint(0,11)/10.
            initial, genotype, phe_decimal = generate_stable_genotype(
                self.bits, self.density, self.basename, self.base,
                globals()[self.genotype_func], p, self.q, self.act_fraction,
                self.binary, self.stable, self.sel_period, self.max_trials)
            ps[phe_decimal].append(p)
        save_file(self.filename + self.extension, ps)
        return ps

    ## Adjusting nx.DiGraph(nx.scale_free_graph) paramaters
    def degree_parameters_sample(self):
        return mean([
            generate_scale_free_genotype(
                self.bits, self.alpha, self.gamma, self.delta_in,
                self.delta_out).in_degree().values()
            for i in arange(self.samples)])


    ## ATTENTION: mutations only! See below for initials
    def perturb_and_analyse(self, generations, pop = None):
        if not pop:
            pop = get_pop(self.shortname, self.run, generations)[0]
        return perturb_and_analyse2(pop, None, self.n_perturb, self.normalize)[0]
        '''
        # all individuals in the population have the same initial state
        initial  = self.pop.ini_decimal
        # Perturbed initials are nearest neighbours
        initials = get_nearest_neighbours(initial, self.n_neighbours)
        # Generate perturbed ensembles
        perturb_and_analyse2(pop, initials, None, self.normalize)[0]
        '''
    ##
    def epistasis(self, generations, pop = None):
        if not pop:
            pop = get_pop(
                self.shortname, self.run, generations,
                extension=self.extension)[0]
        return evaluate_epistasis(pop, self.n_perturb)[0]

    ##
    def robustness(self, generations, pop = None):
        if not pop:
            pop = get_pop(
                self.shortname, self.run, generations,
                extension=self.extension)[0]
        return evaluate_robustness(pop, None, self.n_perturb)

    ##
    def survivability(self, generations, pop = None):
        if not pop:
            pop = get_pop(
                self.shortname, self.run, generations,
                extension=self.extension)[0]
        return evaluate_survival(pop)

    ##
    def autoregulation(self):
        pop = self.get_population()
        autoregulation(pop, self.pop_stats, self.n_perturb)
        save_file(self.filename + self.extension, pop)

    ##
    def generate_stable_pop(self):
        pop = self.get_population()
        save_file(self.filename + self.extension, pop)
        return pop

    ##
    def stable_pop_robustness(self):
        pop = load_file(self.filename + self.extension, self.verbose)
        evaluate_robustness_and_survivability_matrix(pop, 'robustness')
        save_file(self.filename + self.extension, pop)
        return pop

    ##
    def stable_pop_survivability(self):
        pop = load_file(self.filename + self.extension, self.verbose)
        evaluate_robustness_and_survivability_matrix(pop, 'survivability')
        save_file(self.filename + self.extension, pop)
        return pop

    ##
    def superfrankenmatrix(self):
        pop = load_file(self.filename + self.extension, self.verbose)
        mutant_fname = self.filename.replace(
            get_suffix(binary=self.binary,act_fraction=self.act_fraction),
            get_suffix(binary=self.binary,act_fraction=self.p))
        mutant_pop = load_file(mutant_fname + self.extension, self.verbose)
        pop.mutate_diagonal(mutant_pop)
        pop.cycle_develop()
        evaluate_robustness_and_survivability_matrix(
            pop, 'robustness_and_survivability')
        self.filename = self.filename.replace(
            get_suffix(binary=self.binary, act_fraction=self.act_fraction),
            get_suffix(binary=self.binary, act_fraction=self.act_fraction,
                       p=self.p))
        save_file(self.filename + self.extension, pop)
        return pop

    ##
    def qp(self):
        pop = self.get_population()
        pop.cycle_develop()
        if not self.stable:
            pop.s = pop.get_stability_matrix()
        #evaluate_robustness_and_survivability_matrix2(pop)
        save_file(self.filename + self.extension, pop)
        return pop

    ##
    def robustness_and_survivability_matrix(self, gpmap):
        for generation in self.datapoints[1:]:
            pop, fname = get_pop(
                self.shortname, self.run, generation, '', self.extension,
                self.verbose)
            if (pop and pop.size and
                (not (hasattr(pop, 'robustness') and
                      hasattr(pop, 'survivability')) or
                 self.rm_file)):
                pop = self.copy_old_pop_to_new_pop(pop)
                if pop.noise:
                    # for continuous (cSigmoid) phenotypes
                    pop.robustness_threshold = self.robustness_threshold
                evaluate_robustness_and_survivability_matrix2(pop, gpmap=gpmap)
                save_file(fname + self.extension, pop)
        return pop

    ##
    def qp_pop_robustness_and_survivability_matrix(self):

        p = self.p
        if p == .5:
            ps = [p]
        elif p == 'ps_all':
            ps = ps_all
        else:
            ps = arange(0, 1 + p, p)

        q = self.q
        if q == .5:
            qs = [q]
        elif q == 'qs_all':
            qs = qs_all
        else:
            qs = arange(0, 1 + q, q)

        for p, q in itertools.product(ps, qs):

            pop, fname = get_non_evolved_pop(
                'qp.pop', self.stable, self.run, None, p, q, self.sel_period,
                self.pop_size, self.binary, self.genotype_func, self.extension,
                self.verbose)

            if (pop and pop.size and
                (not (hasattr(pop, 'robustness') and
                      hasattr(pop, 'survivability')) or self.rm_file)):

                evaluate_robustness_and_survivability_matrix2(pop)
                save_file(fname + self.extension, pop)

        return 1

    ##
    def qp_initials_robustness_and_survivability(self):
        pop = load_file(self.filename + self.extension, self.verbose)
        evaluate_robustness_and_survivability_matrix2(
            pop, False, False, self.n_neighbours)
        save_file(self.filename + self.extension, pop)
        return pop

    ##
    def stats(self):
        import src.plotting_tools2 as p
        pop_function = p.get_pop_function(
            self.pop_stats, self.experiment, self.shortname, self.n_datapoints)
        return pop_function(
            shortname=self.shortname, stats=self.pop_stats,
            pop_stat_func=self.pop_stat_func,
            datapoints=self.datapoints, start=self.pop_start,
            end=self.run, samples=self.samples, g=self.run,
            n_runs=self.run, period=self.sel_period,
            act_fraction=self.act_fraction, p=self.p, q=self.q,
            pop_size=self.pop_size, stable=self.stable,
            experiment=self.experiment, genotype_func=self.genotype_func,
            verbose=self.verbose, _plotting=False,
            extension=self.extension, _reshape=self.reshape)

    ## Evolution
    def evolution(self, logfile, gpmap):

        # try to read saved pop, starting from the last generation, to 0
        for g in self.datapoints[::-1]:

            pop, fname = get_pop(self.shortname, self.run, g)
            #, '', self.extension, self.verbose, self.experiment)

            if pop:

                # pop was saved before selection if not self.selection
                if not self.selection:
                    pop.selection()

                # del vectors to save space
                for ind in pop.individuals:
                    ind.connectivity = ind.optimum = None

                break

        # if no pop
        if not pop:
            pop = self.get_population()
            fname = self.filename + print_generations(g, self.precision)
            save_file(fname + self.extension, pop)

        # legacy code: copy old class to new class and update pop parameters
        pop = self.copy_old_pop_to_new_pop(pop)

        # Evolution
        for g in range(pop.generation + 1, int(self.generations) + 1):
            # log
            print >> logfile, 'generation', g

            # recombine, mutate, develop and evaluate fitness.
            # selection if myrun.selection
            pop.evolution(self.selection, gpmap)

            # save every x generations
            if g in self.datapoints:
                fname = self.filename + print_generations(g, self.precision)
                save_file(fname + self.extension, pop)
            #print >> logfile, pop

            # if not myrun.selection (do selection if saved pop before selection)
            if not self.selection:
                pop.selection()

        return pop, fname


    ## The output
    def save_file(self, data=None, suffix='', extension='.dat'):

        # save just unique(phenotypes) or all phenotypes
        if 'phenotypes' in self.experiment:
            output = data if not self.only_unique else list(set(data))

        elif self.experiment == 'distribution':
            output = array(data)

        elif 'pop' in self.experiment:
            output = data

        else:
            output = self.stability

        _save = save if extension == '.npy' else save_file
        _save(self.filename + suffix + extension, output)

        del output


if __name__ == "__main__":
#def main():
    #myrun = runs()

    # load the pickled Run instance from stdin
    try:
        myrun = cPickle.load(sys.stdin)
    except:
        myrun = cPickle.load(open('N4c1inf.run'))

    if myrun.verbose:
        print myrun

    #if myrun.experiment == 'stability' or myrun.experiment == 'phenotypes':
     #   myrun.sampling = True

    # set it up
    individuals, initials, gpmap = myrun.setup()

    if myrun.verbose:
        print myrun

    if myrun.save:
        fname = myrun.prefix + ".run"
        file = open(fname,  "w")
        cPickle.dump(myrun, file)
        file.close()

    # write some things to a log file
    fname = myrun.prefix + ".log"
    logfile = open(fname,  "a")
    print >> logfile, 'printing to log (from %s)' % myrun.prefix
    print >> logfile, myrun.filename
    print >> logfile, myrun

    # extra set up
    if myrun.network and type(myrun.fixed_states) == dict:
        myrun.fixed_states = array(
            [(list(myrun.nodes).index(gene), state)
             for gene, state in myrun.fixed_states.iteritems()]).T

    # run it
    if 'stability' in myrun.experiment:
        n = myrun.stability_sample(individuals, initials, gpmap)
        # no need to save: save_txt2() in model.get_stability_and_phenotypes()
    elif myrun.experiment == 'distribution':
        n = myrun.stability_sample(individuals, initials, gpmap)
        savetxt(myrun.filename + '.txt', n, '%f %f %d')
    elif 'phenotypes' in myrun.experiment:
        n, phenotypes = myrun.phenotype_sample(individuals, initials, gpmap)
        #myrun.save_file(phenotypes, extension='.npy')
    elif myrun.experiment == 'path_length':
        lengths = myrun.path_length_vs_N(individuals, gpmap)
        n = lengths.size
        # no need to save: save_txt3() in model.path_length_vs_N2()
        #myrun.save_file(lengths, extension='.npy')
        #savetxt(myrun.filename + '.txt', lengths, '%d ')
    elif myrun.experiment == 'viability':
        n = myrun.viability_sample(individuals, gpmap) / float(myrun.n_mutants)
    elif myrun.experiment == 'final_vs_p':
        data = myrun.phenotype_vs_p()
    elif myrun.experiment == 'autoregulation':
        myrun.autoregulation()
    elif 'stats' in myrun.experiment:
        data = myrun.stats()
    elif myrun.experiment == 'stable.pop':
        pop = myrun.generate_stable_pop()
    elif myrun.experiment == 'stable.pop.evaluate':
        pop = myrun.generate_stable_pop()
        evaluate_robustness_and_survivability_matrix2(pop)
        save_file(myrun.filename + myrun.extension, pop)
    elif myrun.experiment == 'stable.pop.robustness':
        pop = myrun.stable_pop_robustness()
    elif myrun.experiment == 'stable.pop.survivability':
        pop = myrun.stable_pop_survivability()
    elif myrun.experiment == 'stable.pop.superfrankenmatrix':
        pop = myrun.superfrankenmatrix()
    elif myrun.experiment == 'qp.pop':
        pop = myrun.qp()
    elif myrun.experiment == 'qp.pop.initials':
        pop = myrun.qp_initials_robustness_and_survivability()
    elif myrun.experiment == 'qp.pop.evolution':
        pop = myrun.evolution(logfile, gpmap)[0]
    elif myrun.experiment == 'pop.evolution':
        pop = myrun.evolution(logfile, gpmap)[0]
    elif myrun.experiment == 'pop.evolution.u0':
        pop, fname = myrun.evolution(logfile, gpmap)
        save_file(
            fname + '-u0g%d' %(myrun.generations-1e6) + myrun.extension, pop)
    elif myrun.experiment == 'pop.robustness_and_survivability':
        pop = myrun.robustness_and_survivability_matrix(gpmap)
    elif myrun.experiment == 'qp.pop.robustness_and_survivability':
        myrun.qp_pop_robustness_and_survivability_matrix()
    elif 'pop' in myrun.experiment:
        for generation in myrun.datapoints[1:]:

            pop, pop_fname = get_pop(
                myrun.shortname, myrun.run, generation, '', myrun.extension,
                myrun.verbose)
            data, data_fname = get_pop(
                myrun.shortname, myrun.run, generation, '', myrun.extension,
                myrun.verbose, myrun.experiment.split('.')[-1])

            # pop exists but there is still no data: run
            if pop and not any(data):
                if myrun.experiment == 'pop.perturb_and_analyse':
                    data = myrun.perturb_and_analyse(generation, pop)
                elif myrun.experiment == 'pop.epistasis':
                    data = myrun.epistasis(generation, pop)
                elif myrun.experiment == 'pop.robustness':
                    data = myrun.robustness(generation, pop)
                elif myrun.experiment == 'pop.survivability':
                    data = myrun.survivability(generation, pop)
                if any(data):
                    myrun.save_file(
                        data, print_generations(generation, myrun.precision),
                        myrun.extension)

            # pop doesnt exist but there's data:
            # rm data (m4r05b corn/evo1 rsync fiasco, rm pop..)
            elif not pop and any(data) and myrun.rm_file:
                rm_file(data_fname, myrun.verbose)

    elif myrun.experiment == 'scale_free_graph':
        if ((1-myrun.alpha-myrun.gamma) > 0 and
            (1-myrun.alpha-myrun.gamma+myrun.alpha+myrun.gamma) == 1):
            n = myrun.degree_parameters_sample()
            file = open(data_dir + myrun.prefix + '.txt', 'w')
            print >> file, n
            file.close()
        else:
            n = 0

    # write to log file
    if 'pop.evolution' in myrun.experiment:
        print >> logfile, 'generations', pop.generation
    elif 'stats' in myrun.experiment:
        if any(data) and type(data) == dict:
            print [(key, ma.masked_invalid(value).mean())
                   for key, value in data.items()]
    elif myrun.experiment == 'pop.epistasis':
        if any(data):
            print >> logfile, ma.masked_invalid(data).mean()
    elif myrun.experiment == 'pop.robustness':
        if any(data):
            print >> logfile, ma.masked_invalid(data.equal + data.intersect).mean()
    elif myrun.experiment == 'pop.survivability':
        if any(data):
            print >> logfile, ma.masked_invalid(data).mean()
    elif myrun.experiment == 'pop.perturb_and_analyse':
        pass #print >> logfile, zip(data.keys(), map(mean, data.values()))
    elif myrun.experiment == 'stable.pop':
        if pop:
            print >> logfile, pop.size
    elif myrun.experiment == 'stable.pop.robustness':
        if pop:
            print >> logfile, pop.robustness.mean()
    elif myrun.experiment == 'stable.pop.survivability':
        if pop:
            print >> logfile, pop.survivability.mean()
    elif (myrun.experiment == 'stable.pop.superfrankenmatrix' or
          myrun.experiment == 'stable.pop.evaluate'):
        if pop:
            print >> logfile, pop.robustness.mean(), pop.survivability.mean()
    elif (myrun.experiment == 'qp.pop' or
          myrun.experiment == 'qp.pop.robustness_and_survivability'):
        pass
        #if pop: print >> logfile, pop.survivability.mean() #pop.robustness.mean()
    elif myrun.experiment == 'pop.robustness_and_survivability':
        if pop:
            print >> logfile, pop.robustness.mean(), pop.survivability.mean()
    elif myrun.experiment == 'qp.pop.initials':
        if pop:
            print >> logfile, pop.robustness.mean(), pop.survivability.mean()
    elif myrun.experiment == 'final_vs_p':
        print >> logfile, ma.masked_invalid(data).mean()
    elif myrun.experiment != 'autoregulation':
        print >> logfile, n, myrun.samples, n/float(myrun.samples)

    logfile.close()
