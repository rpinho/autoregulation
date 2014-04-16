import src.model
reload(src.model)
from src.model import *

pa_stats = ('stable',
            'equal', 'equal fix', 'equal cycle',
            'unique', 'unique fix', 'unique cycle',
            'entropy', 'entropy fix', 'entropy cycle',
            'hamming', 'epistasis')

#psyco.full() # compile as much code as possible to speedup simulation

class Individual:
    def __repr__(self):
        return str(self.__dict__)

    def __init__ (
            self, genotype = None, initial = None,
            basename = def_basename, gpmap = def_gpmap, base = def_base,
            bits = N, density = c, phe_decimal = None, act_fraction = None,
            binary = None, stable = None, period = None,
            genotype_function = generate_genotype,
            random = False, conditional = False, ind = None,
            max_trials = F, get_orbits = False, dtype = float32):

        # some default behaviors: see each case for variations from this
        self.phenotype    = None # vector
        self.phe_decimal  = phe_decimal
        self.optimum      = None # vector
        self.opt_decimal  = None
        self.length       = None
        self.stable       = stable
        self.period       = period
        self.orbit        = None
        self.fitness      = None
        self.dtype        = None
        self.connectivity = None
        self.p            = act_fraction
        self.binary       = binary

        # default/blank
        # given (default is None)
        self.genotype     = genotype
        # given (default is None). Note: initial should come in decimal code
        self.ini_decimal  = initial
        # Note (see below): if any(self.genotype): self.size = len(self.genotype)
        self.size         = bits

        # base maps / state space
        self.base         = base
        self.basename     = basename
        self.gpmap        = gpmap if type(gpmap) == str else gpmap.__name__
        self.stepfunction = None#gpmap.__doc__

        if random:
            self.genotype = genotype_function(self.size, density, dtype = dtype)
            self.ini_decimal = generate_decimal_phenotype(self.base, self.size)

            self.period, self.length, self.phe_decimal = brent(
                devo, self.genotype, self.ini_decimal,
                self.basename, gpmap, self.base, self.size)

            if self.period > 1:
                if get_orbits:
                    self.orbit = get_orbit(
                        self.period, self.genotype, self.phe_decimal,
                        self.basename, gpmap, self.base, self.size)
                self.phe_decimal = None
                self.stable      = False
            else:
                self.stable      = True

        elif conditional:
            self.period = period
            self.ini_decimal, self.genotype, self.phe_decimal = generate_stable_genotype(
                self.size, density, self.basename, gpmap, self.base,
                genotype_function, p, q, act_fraction, binary, stable,
                self.period, max_trials, dtype)
            if self.period > 1:
                self.stable = False
            else:
                self.stable = True

        # Common to all constructs
        if any(self.genotype):
            self.size         = len(self.genotype)              # bits
            self.dtype        = self.genotype.dtype
            #self.connectivity = get_connectivity(self.genotype) # Modifier allele

        # copy another instance
        elif ind:
            # new pop has the same number of keys/attrib as old pop
            #self.__dict__ = pop.__dict__
            # new pop has the same number of keys/attrib as a new instance
            # could also use getattr() and setattr()
            for key in ind.__dict__.keys():
                # values/variables
                self.__dict__[key] = ind.__dict__[key]
            for key in ind.__class__.__dict__.keys():
                if '__' not in key:
                    # methods
                    self.__class__.__dict__[key] = ind.__class__.__dict__[key]

            del ind
            #return None

    def cycle_develop(self, pop, gpmap = None, devo_function = devo):

        if self.stepfunction:
            step  = globals()[self.stepfunction]
            jumps = get_jumps(self.base.size)
            def gpmap(x): return step(x, jumps, self.base)
        else:
            gpmap = globals()[self.gpmap] if self.gpmap != 'gpmap' else gpmap

        self.period, self.length, self.phe_decimal, self.phenotype = general_develop(
            self.genotype, self.ini_decimal, self.basename, gpmap, self.base,
            self.size, pop.noise, devo_function, pop.converge_time,
            pop.converge_threshold, pop.devo_time)

#    def get_stability(self):
#        self.stable = True if self.period == 1 else False
#        return self.stable

    # ind_functions

    # returns True or False
    def get_stability(self, sel_period = 1):
        return self.period == sel_period

    # returns True or False or NaN
    # self is mutant
    def get_robustness(self, wildtype, gpmap, mutant_pop, sel_period=1):
        # contingent on self and wildtype stability (sel_period = 1)
        #return (self.phe_decimal == wildtype.phe_decimal)*self.period*wildtype.period == 1

        if not (self.period and
                (self.phe_decimal != None or any(self.phenotype))):
            self.cycle_develop(mutant_pop, gpmap)

        # wildtype is stable (sel_period = 1)
        if wildtype.period == sel_period:

            # mutant and wildtype are stable and equal
            if self.period == wildtype.period:

                if (self.phe_decimal != None and wildtype.phe_decimal != None):
                    return self.phe_decimal == wildtype.phe_decimal

                elif (any(self.phenotype) and any(wildtype.phenotype)):
                    return convergence([self.phenotype, wildtype.phenotype],
                                       mutant_pop.converge_threshold)
            else: return False

        # wildtype is not stable
        else: return NaN

    def get_connectivity(self):
        return get_connectivity(self.genotype)

    def get_activating_fraction(self):
        return get_activating_fraction(self.genotype)

    def get_activating_fraction_diagonal(self):
        return get_activating_fraction_diagonal(self.genotype)

    def get_activating_fraction_non_diagonal(self):
        return get_activating_fraction_non_diagonal(self.genotype)

    def get_number_of_positive_diagonal(self):
        return get_number_of_positive_diagonal(self.genotype)

    def get_number_of_positive_non_diagonal(self):
        return get_number_of_positive_non_diagonal(self.genotype)

    # label returns labels, n
    def get_number_pos_clusters(self):
        return label_clusters(self.genotype,-1)[1]

    def get_number_neg_clusters(self):
        return label_clusters(self.genotype, 1)[1]

    # bin[0] is size 0; bin[1:].size is # clusters
    def get_size_pos_clusters(self):
        return bincount(label_clusters(self.genotype,-1)[0].flatten())[1:]

    def get_size_neg_clusters(self):
        return bincount(label_clusters(self.genotype, 1)[0].flatten())[1:]


class Population:
    def __repr__(self):
        return str(self.__dict__)

    def __init__ (
            self, founder_function = siegals_founder, size = P,
            stable = True, dtype = float32,
            founder = None, initials = None, mutations = None,
            random = False, ini_random = False, initial = None,
            binary = False, conditional = False, density = c,
            basename = def_basename, gpmap = None, base = def_base,
            genotype_function = None,
            genotypic_space_function = generate_full_genotype_space,
            bits = N, _empty = False, full_enum = False,
            sel_strength = s, min_fit = zero_fit, sel_period = 1,
            rec_rate = r, mut_rate = u,
            ins_rate = uA, del_rate = uD, all_mutants = False,
            change_sign = False, deletion = False, mut_bias = False,
            rec_function = 'recombine', act_fraction = p,
            p = p, q = p, max_trials = F, deterministic = True,
            pop = None,
            devo_time = inf, noise = 0, steepness = a, converge_time = T,
            converge_threshold = converge):

        # main attribute of this class
        self.individuals   = []

        # some default behaviors: see each case for variations from this
        self.initial       = None # vector
        self.ini_decimal   = None
        self.optimum       = None # vector
        self.opt_decimal   = None
        self.stable        = None
        self.length        = None
        self.dtype         = None

        # default
        self.genotype      = founder # given (default is None)
        self.bits          = bits
        self.density       = density

        # random or conditional (stable xor unstable) individuals
        self.random        = random
        self.binary        = binary
        self.full_enum     = full_enum

        # development process parameters
        self.devo_time          = devo_time
        self.noise              = noise
        self.steepness          = steepness
        self.converge_time      = converge_time
        self.converge_threshold = converge_threshold

        # generating genotypes parameters
        self.act_fraction  = act_fraction
        self.p             = p
        self.q             = q
        self.deterministic = deterministic

        # generating function for sampling genotypes
        if not genotype_function:

            if not self.binary:
                genotype_function = generate_genotype

            else:
                genotype_function = spin_genotype

        # base maps / state space
        # base should be None if you want real-valued phenotypes
        # between min_state and max_state
        self.base         = base
        self.basename     = basename

        # gpmap
        if gpmap:
            if type(gpmap) == str or NoneType:
                self.gpmap = gpmap
            else:
                self.gpmap = gpmap.__name__
        else:
            gpmap = get_gpmap(
                self.base, self.basename, self.density, self.bits, self.binary,
                get_jumps(self.base.size))

            self.gpmap = gpmap.__name__

        ## Evolution
        self.generation   = 0
        # selection
        self.sel_strength = sel_strength
        self.sel_period   = sel_period
        self.min_fit      = min_fit
        # recombination
        self.rec_rate     = rec_rate
        self.rec_function = rec_function
        # mutations
        self.mut_rate     = mut_rate
        self.change_sign  = change_sign
        #if self.binary: self.change_sign = True
        self.mut_bias     = mut_bias
        self.deletion     = deletion
        self.ins_rate     = ins_rate
        self.del_rate     = del_rate
        self.all_mutants  = all_mutants

        # empty instance
        if _empty:
            #return None
            self.individuals = []
            self.size = 0

        # copy another instance
        elif pop:
            # new pop has the same number of keys/attrib as old pop
            #self.__dict__ = pop.__dict__
            # new pop has the same number of keys/attrib as a new instance
            # could also use getattr() and setattr()
            for key in pop.__dict__.keys():
                self.__dict__[key] = pop.__dict__[key] # values/variables

            # methods
            for key in pop.__class__.__dict__.keys():
                if '__' not in key:
                    self.__class__.__dict__[key] = pop.__class__.__dict__[key]

            #for ind in pop.individuals: ind = Individual(ind = ind)
            del pop
            #return None

        ## Founder is a given genotype
        # if mutations: ALL members of the same population are single-mutants
        # from the wildtype founder
        # else: ALL members of the same population have the SAME genotype
        # (isogenic population)
        elif any(self.genotype):

            self.dtype = self.genotype.dtype

            self.bits = len(self.genotype) # N

            if mutations:
                initials = [initial]

            self.set_perturbed_individuals(initials, mutations)

            # if mutations: size = mutations, else: size L = k^N or size given
            self.size = len(self.individuals)
            del initials, founder


        ## Null model: not conditional on initial conditions
        #  (i.e. population founder).
        #  ALL members of the same population have DIFFERENT genotypes
        elif self.random:

            self.size  = int(size)
            self.dtype = dtype

            # all individuals in the population have different initial states
            if ini_random:
                initials = (generate_decimal_phenotype(self.base, self.bits)
                            for i in range(self.size))

            # all individuals in the population have the same initial state
            else:
                if initial:
                    self.ini_decimal = initial
                else:
                    self.ini_decimal = generate_decimal_phenotype(
                        self.base, self.bits)
                initials = (self.ini_decimal for i in range(self.size))

            self.set_random_individuals(genotype_function, initials)

            del initials

            # optimum for fitness (model 4)
            if not self.noise:
                self.opt_decimal = generate_decimal_phenotype(
                    self.base, self.bits)
                self.optimum = get_vector(
                    self.opt_decimal, self.basename, self.base, self.bits)

            else:
                self.optimum = generate_vector_phenotype(
                    self.base, self.bits, self.noise) # real_phenotype()

            # Notes: the optimum or opt_decimal are random because it's not
            # conditional on the founder
            # They do not necessarilly relate to the phenotype of the individual
            # who later founds an ensemble's population.
            # It's not the wild-type phenotype (although it could be with strong
            # selection). It's USED ONLY FOR SELECTION


        ## Full enumeration:
        #  ALL members of the population have DIFFERENT genotypes
        elif self.full_enum:

            # Binary networks (-1, 1)
            self.individuals = [
                Individual(
                    genotype, None, self.basename, gpmap, self.base)
                for genotype in genotypic_space_function(self.bits)]
            self.size   = len(self.individuals)
            self.binary = True


        ## All individuals (size N) in the population have DIFFERENT genotypes
        #  and initial states
        elif conditional:
            self.stable = stable
            self.size   = size
            self.dtype  = dtype
            self.max_trials = max_trials # for generate_stable_genotype()

            self.set_stable_individuals(genotype_function)
            self.size = len(self.individuals)


        # Founder is an individual conditional on stable
        else:
            ## ATTENTION: this might be outdated/broken
            founder = Individual(
                founder_function, stable, dtype,
                basename = self.basename, gpmap = gpmap, base = self.base)
            self.size = int(size)
            self.individuals = [copy.copy(founder) for i in range(self.size)]

            # ALL members of the same population have the SAME genotype
            # and initial and optimum phenotypes
            self.genotype    = founder.genotype # the founder genotype
            self.initial     = founder.initial
            self.ini_decimal = founder.ini_decimal
            self.optimum     = founder.optimum  # Note: could be None
            self.opt_decimal = founder.opt_decimal

            self.stable      = stable         # population mean stability
            self.length      = founder.length # population mean path length
            self.dtype       = dtype

            del founder

        # Auxiliary variables for statistical analysis and evolutionary processes
        self.phenotypes   = None # auxiliary list for statistical analysis
        self.equal        = None
        self.entropy      = None
        self.connectivity = None
        self.M1           = None
        self.c            = None
        self.densities    = None
        self.hamming      = None # mean hamming distance to self.optimum
        self.fitness      = None # mean fitness in the population
        self.max_fitness  = None # max fitness in the population
        self.stables      = None # evolutionary record of mean stability

    def set_perturbed_individuals(self, initials, mutations):

        if mutations:
            # returns ((initial, mutant))
            mutants = get_mutants(
                initials[0], self.genotype, mutations, self.all_mutants,
                self.change_sign, self.deletion)

            self.individuals = [
                Individual(
                    mutant, initial, self.basename, self.gpmap, self.base)
                for initial, mutant in mutants]
        else:
            self.individuals = [
                Individual(
                    self.genotype, ini_decimal, self.basename, self.gpmap,
                    self.base)
                for ini_decimal in initials]

    def set_stable_individuals(self, genotype_function):

        individuals = get_individuals(
            self.full_enum, self.stable, self.size, self.density,
            self.basename, self.base, self.bits, genotype_function,
            p=self.p, q=self.q, act_fraction=self.act_fraction,
            binary=self.binary, sel_period=self.sel_period,
            max_trials=self.max_trials)[0]

        self.individuals = [
            Individual(
                genotype, initial, self.basename, self.gpmap, self.base,
                self.bits, self.density, phe_decimal, self.act_fraction,
                self.binary, self.stable, self.sel_period)
            for initial, genotype, phe_decimal in individuals]

    def set_random_individuals(self, genotype_function, initials):

        self.individuals = [
            Individual(
                genotype_function(
                    self.bits, self.density, p = self.p, q = self.q,
                    act_fraction = self.act_fraction,
                    deterministic = self.deterministic, dtype = self.dtype),
                initial, self.basename, self.gpmap, self.base)
            for initial in initials]

    # EVOLUTION #

    def recombination(self):
        recombine = globals()[self.rec_function]

        for ind1, ind2 in zip(self.individuals[0::2], self.individuals[1::2]):
            recombine(ind1, ind2, self.rec_rate)

    def mutation(self):

        for individual in self.individuals:
            individual.genotype = mutate(
                individual.genotype, rate = self.mut_rate,
                change_sign = self.change_sign, deletion = self.deletion,
                bias = self.mut_bias)

    def insertion_and_deletion(self):

        for individual in self.individuals:
            individual.genotype = insert_and_delete(
                individual.genotype, self.ins_rate, self.del_rate)
            individual.connectivity = get_connectivity(individual.genotype)

    def cycle_develop(self, gpmap=None, devo_function=devo):

        for individual in self.individuals:
            individual.cycle_develop(self, gpmap, devo_function)

    def evaluate_stability(self):

        self.stable = array([
            individual.get_stability()
            for individual in self.individuals])

    def evaluate_fitness(
            self, distance_function=normed_vector_hamming_distance,
            gpmap=None, fitness=gaussian_fitness):

        self.max_fitness = 0
        for individual in self.individuals:

            # the period NOT selected
            if self.sel_period and individual.period != self.sel_period:
                # Wagner = exp(-1/s) ~ 0, Siegal = 0, random = 1
                individual.fitness = self.min_fit

            # the period selected
            else:
                if self.sel_strength == inf:
                    individual.fitness = 1
                else:

                    if not any(individual.phenotype):

                        if not self.noise:
                            if individual.period == 1:
                                individual.phenotype = (
                                    get_vector(
                                        individual.phe_decimal,
                                        individual.basename,
                                        individual.base, individual.size),)

                            else:
                                orbit =  get_orbit(
                                    individual.period, individual.genotype,
                                    individual.phe_decimal,
                                    individual.basename,
                                    globals()[individual.gpmap],
                                    individual.base, individual.size)

                                individual.phenotype = (
                                    get_vector(
                                        phe_decimal, individual.basename,
                                        individual.base, individual.size)
                                    for phe_decimal in orbit)

                                del orbit

                        else:
                            individual.cycle_develop(self, gpmap)

                    # list for cycle in max(fitness)
                    #individual.phenotype = individual.phenotype.reshape(
                     #   -1, individual.size)

                    if not any(individual.optimum):
                        individual.optimum = self.optimum

                    # Generalization of Wagner's model to include cycles
                    individual.fitness = max([
                        fitness(
                            phenotype, individual.optimum, distance_function,
                            self.sel_strength)
                        for phenotype in individual.phenotype])

            # max fitness in the population for normalization
            self.max_fitness = max(self.max_fitness, individual.fitness)

            ## ATTENTION: del vectors to save space
            individual.phenotype = individual.optimum = None

    def evaluate_stability_fitness_hamming_and_M1(
            self, distance_function=int_vector_hamming_distance,
            correct4zeros=True):

        # arbitrary individual of the ensemble
        size = self.individuals[0].size # == N

        self.stable       = 0
        self.fitness      = 0
        self.hamming      = 0
        self.M1           = 0
        self.c            = 0
        self.connectivity = zeros(size)
        self.densities    = zeros(size + 1)

        for individual in self.individuals:

            if not any(individual.phenotype):
                individual.phenotype = get_vector(
                    individual.phe_decimal, individual.basename,
                    individual.base, individual.size)

            # fixpoints only
            if (individual.period == 1 and not
                (correct4zeros and not len(nonzero(individual.phenotype)[0]))):

                individual.stable = True

                self.stable  += 1
                self.hamming += distance_function(
                    individual.phenotype, individual.optimum)

            # (fixpoints and cycles)
            if individual.fitness:
                self.fitness += individual.fitness

            ## Modifier allele

            # M1 is the frequency of the wildtype allele (self.density)
            self.M1 += len(where(individual.connectivity == self.density)[0])

            # computes the connectivity vector of the average individual
            self.connectivity += individual.connectivity

            # computes the frequency distribution of alleles c
            for density in rint(individual.connectivity*individual.size):
                self.densities[density] += 1

            # computes the density of the average network in the population
            self.c += individual.connectivity.mean()

            individual.phenotype = None

        # conditional on stable (fixpoints only)
        if self.stable:
            self.hamming /= float(self.stable)
            if self.fitness: self.fitness /= float(self.stable)

        # absolute in the whole population (fixpoints and cycles)
        self.stable       /= float(self.size)
        self.connectivity /= float(self.size)
        self.c            /= float(self.size)
        self.M1           /= float(self.size * individual.size)
        self.densities    /= float(self.size * individual.size)

    def selection(self):

        new_pop = []

        # if there's any individual with fitness not zero, do selection;
        # otherwise: extinction
        if self.max_fitness:
            while len(new_pop) < self.size:
                individual = self.individuals[randint(self.size)] # sampling
                fitness  = individual.fitness
                fitness /= self.max_fitness
                if uniform() < fitness:
                    new_pop.append(copy.copy(individual))

        del self.individuals
        self.individuals = new_pop
        self.size = len(self.individuals)
        del new_pop

    def evolution(
            self, selection=True, gpmap=None, devo_function=devo,
            distance_function=normed_vector_hamming_distance,
            computer_other_properties=False):

        self.generation += 1

        # the if statements are just to save speed; they are not necessary
        if self.rec_rate:
            self.recombination()

        self.mutation()

        if self.ins_rate or self.del_rate:
            self.insertion_and_deletion()

        self.cycle_develop(gpmap, devo_function)

        self.evaluate_fitness(distance_function, gpmap)

        if computer_other_properties:
            self.evaluate_stability_fitness_hamming_and_M1()

        if selection:
            self.selection()

        #if not self.size:
         #   sys.exit("Population has gone extint at generation %s." %self.generation)

    # mutant is diagonal
    def mutate_diagonal(self, pop):

        for individual, mutant in zip(self.individuals, pop.individuals):
            individual.genotype = ldu_composition_genotype(
                individual.genotype, mutant.genotype)

    def get_genotypes(self, flat=False, _tuple=False):

        if flat:

            if _tuple:
                return [tuple(
                    individual.genotype.flatten())
                        for individual in self.individuals]

            else:
                return array([
                    individual.genotype.flatten()
                    for individual in self.individuals])

        else:
            return array([
                individual.genotype
                for individual in self.individuals])

    def get_connectivities(self):

        return array([individual.get_connectivity()
                      for individual in self.individuals])

    def get_activating_fraction_diagonal(self):

        return array([individual.get_activating_fraction_diagonal()
                      for individual in self.individuals])

    # get_ind_stat

    def get_stability_matrix(self, gpmap = None):

        if not self.individuals[0].period:
            self.cycle_develop(gpmap) # needed for individual.get_stability()

        return array([individual.get_stability(self.sel_period)
                      for individual in self.individuals])

    # self is mutant: if wildtype is not stable, returns a vector of NaNs
    def get_robustness_matrix(self, wildtype, gpmap = None):

        return array([individual.get_robustness(wildtype, gpmap, self)
                      for individual in self.individuals])

    def set_survivability_mean(self, prefix = ''):

        setattr(self, prefix + 'survivability',
                array([getattr(individual, prefix + 'survivability').mean()
                       for individual in self.individuals]))

    # individual.robustness is either a boolean vector or a vector of all NaNs,
    # so mean() is fine
    def set_robustness_mean(self, prefix = ''):

        setattr(self, prefix + 'robustness',
                array([getattr(individual, prefix + 'robustness').mean()
                       for individual in self.individuals]))

class Ensemble:
    def __init__ (
            self, founder_function=siegals_founder, size=n_runs,
            samples=n_perturb, stable=1, dtype=float32, founder=None,
            initials=None, population=None, mutations=None,
            all_mutants=False, n_neighbours=0,
            basename=def_basename, gpmap=def_gpmap, base=def_base):

        self.size   = size
        self.random = None

        # receives ONE founder GENOTYPE and creates SINGLE-MUTANTS
        if any(founder):
            # full enumeration of all initial phenotypes = k^N +1
            genotype = founder
            self.dtype = genotype.dtype
            self.populations = [
                Population(
                    founder=genotype, initials=initials,
                    basename=basename, gpmap=gpmap, base=base)]

            for i in range(1, self.size):
                # each ensemble member genotype differs from the previous by
                # one-point mutation
                genotype = mutate(genotype, 1)
                self.populations.append(
                    Population(
                        founder=genotype, initials=initials,
                        basename=basename, gpmap=gpmap, base=base))

        # receives a POPULATION of founder INDIVIDUALS
        # if mutations: each individual wildtype founds a population of
        # single-mutants with the same initial phenotype
        # else: each individual founds a population of isogenic individuals but
        # with different initial phenotypes
        elif any(population):
            self.populations = []
            for individual in population.individuals:
                # Perturbed initials are nearest neighbours or
                # generate_full_phenotype_space()
                if n_neighbours:
                    initials = get_nearest_neighbours(
                        individual.ini_decimal, n_neighbours)
                else:
                    initials = generate_full_phenotype_space(
                        0, individual.basename, individual.base,
                        individual.size)

                perturbed_pop = Population(
                    founder=individual.genotype, initials=initials,
                    mutations=mutations,
                    initial=individual.ini_decimal, binary=population.binary,
                    basename= population.basename, gpmap=population.gpmap,
                    base=population.base,
                    sel_period=population.sel_period, all_mutants=all_mutants,
                    devo_time=population.devo_time, noise=population.noise,
                    steepness=population.steepness,
                    converge_time=population.converge_time,
                    converge_threshold=population.converge_threshold)

                perturbed_pop.ini_decimal = individual.ini_decimal
                # == None (ind) if population.random
                perturbed_pop.opt_decimal = individual.opt_decimal

                self.populations.append(perturbed_pop)
                del perturbed_pop

            # the target used for selection,
            # not necessarilly the wild-type (Elhanan's original) phenotype
            self.opt_decimal = population.opt_decimal
            self.random      = population.random
            self.size        = population.size
            del population

        else:
            self.populations = [
                Population(
                    founder_function, samples, stable, dtype)
                for i in range(self.size)]

            self.samples = samples
            self.stable  = stable

    def develop(
            self, start_over=False, noise=1e-100, sample_init=False,
            devo_time=M, window=T, gpmap=def_gpmap):

        time = devo_time if devo_time > window else window
        for population in self.populations:
            #initial = population.initial if start_over else None
            for individual in population.individuals:
                initial = individual.initial if start_over else None
                if sample_init:
                    initial = spin_phenotype()
                individual.develop(initial, gpmap, noise, time, window)

    def cycle_develop(self, gpmap = None):
        for population in self.populations: population.cycle_develop(gpmap)

    # Note: NO test for stability, just successive iterations with noise
    def noise_develop(self, noise, devo_time=M):
        for population in self.populations:
            for individual in population.individuals:
                for t in range(devo_time):
                    individual.ini_decimal = devo_noise(
                        individual.genotype, individual.ini_decimal, noise,
                        individual.basename, globals()[individual.gpmap],
                        individual.base)

    '''def print_phenotypes(self, i = None, j = None):
        i = i if i else self.size
        j = j if j else self.samples
        for population in self.populations[:i]:
            print population.optimum, population.opt_decimal, '\n'
            for individual in population.individuals[:j]:
                print individual.phenotype, individual.phe_decimal, individual.stable
            print '\n'
            '''
    # self is mutant; pop is wildtype
    def set_ind_survivability_matrix(self, pop, prefix = ''):
        # wildtype is Individual of pop; mutants is Population of self
        for wildtype, mutants in zip(pop.individuals, self.populations):
            setattr(
                wildtype, prefix + 'survivability',
                mutants.get_stability_matrix(gpmap))

    # self is mutant; pop is wildtype
    def set_ind_robustness_and_survivability_matrix(
            self, pop, prefix='', gpmap=None):

        # wildtype is Individual of pop; mutants is Population of self
        for wildtype, mutants in zip(pop.individuals, self.populations):

            if not (wildtype.period and
                    (wildtype.phe_decimal or any(wildtype.phenotype))):
                wildtype.cycle_develop(pop, gpmap)

            if pop.noise:
                mutants.robustness_threshold = pop.robustness_threshold

            setattr(
                wildtype, prefix + 'robustness',
                mutants.get_robustness_matrix(wildtype, gpmap))

            setattr(
                wildtype, prefix + 'survivability',
                mutants.get_stability_matrix(gpmap))



class PopulationAnalysis():

    def __init__ (self, population):

        self.population = population # main attribute of this class
        self.entropy    = None

    def evaluate_genotype_diversity(self):

        genotypes = [
            tuple(
                reshape(
                    individual.genotype, individual.size * individual.size))
            for individual in self.population.individuals]
        unique_genotypes = set(genotypes)
        #[x for x in genotypes if x not in locals()["_[2]"]]

        entropy = 0.
        for genotype in unique_genotypes:
            frequency = genotypes.count(genotype) / float(self.population.size)
            entropy  += frequency * log(frequency)

        if entropy: entropy /= -log(self.population.size)
        #population.entropy = entropy
        #self.entropy.append( entropy)
        return entropy

    def evaluate_phenotype_diversity(self):

        phenotypes = [individual.phe_decimal
                      for individual in self.population.individuals]

        entropy = 0.
        for phenotype in unique(phenotypes):
            frequency = phenotypes.count(phenotype) / float(self.population.size)
            entropy  += frequency * log(frequency)

        if entropy: entropy /= -log(self.population.size)
        #population.entropy = entropy
        #self.entropy.append( entropy)
        return entropy


class EnsembleAnalysis:
    def __init__ (self, ensemble):

        # main attribute of this class
        self.ensemble      = ensemble

        # main average measures on phenotypes
        self.length        = None
        self.stable        = None
        self.entropy       = None
        self.entropy_fix   = None
        self.entropy_cycle = None
        self.hamming       = None
        self.gain          = None

        # integer measures on phenotypes
        self.unique        = None

        # for noise interval
        self.noises        = None
        self.stable_mean   = None
        self.stable_std    = None

        # helper lists to store past information (evolutionary records)
        self.lengths       = None
        self.stables       = None
        #self.equals        = None
        self.entropies     = None
        self.uniques       = None

        # pixel map, all against all, 1 point per individual:
        # shape (ensemble.size, population.size)
        self.period        = None
        self.phenotypes    = None
        self.robust        = None

        # list of unique orbits in the phenotypic (state) space
        self.orbits        = None

    def evaluate_stability(self):

        self.stable = []
        for population in self.ensemble.populations:
            stable = 0
            for individual in population.individuals:
                if individual.stable: stable +=1
            stable /= float(population.size)
            population.stable = stable
            self.stable.append(stable)

        if not any(self.stables): self.stables = []
        self.stables.append(self.stable)

    def evaluate_path_length(self):

        self.length = []
        for population in self.ensemble.populations:
            length = []
            for individual in population.individuals:
                if individual.stable:
                    length.append(individual.length)
            if length:
                length = mean(length)
                population.length = length
                self.length.append(length)

        if not any(self.lengths): self.lengths = []
        self.lengths.append(self.length)

    def noise_interval(self, max, step, gpmap = def_gpmap, stable = ''):

        self.stable_mean = []
        self.stable_std  = []
        self.noises = arange(step, max, step)

        mod_name = get_model_name(gpmap)
        size = self.ensemble.size
        samples = self.ensemble.populations[0].size
        mid_name = parameters
        name = ensembles_dir + '%(stable)sensemble%(size)s_population%(samples)s%(mod_name)s%(mid_name)s_10noise' % locals()

        for noise in self.noises:
            self.ensemble.develop(True, noise, gpmap = gpmap)

            filename = '%s%s.dat' %(name, noise)
            file     = open(filename, 'wb')
            cPickle.dump(self.ensemble, file)
            file.close()

            self.evaluate_stability()
            self.stable_mean.append( mean(self.stable))
            self.stable_std.append(  std(self.stable))

    # get fixed points
    def get_stable_phenotypes(self):

        for population in self.ensemble.populations:
            population.phenotypes = [
                individual.phe_decimal
                for individual in population.individuals if individual.stable]

    def compare_phenotypes(self, normed = True, cycles = True, random = False):

        # if only fixed points, create the lists population.phenotypes
        # get arbitrary population of the ensemble
        test_pop = self.ensemble.populations[0]
        if not cycles and not test_pop.phenotypes:
            self.get_stable_phenotypes()

        # if cycles are included in the analysis,
        # create the array self.phenotypes
        if cycles and not any(self.phenotypes):
            self.pixel() # warning: slow

        # test for full enumeration of initial phenotypes
        # get arbitrary individual of the ensemble
        test_ind  = test_pop.individuals[0]
        full_enum = True if test_pop.size == pow(test_ind.base.size, test_ind.size) else False

        self.equal   = []
        self.entropy = []
        self.unique  = []
        for i, population in enumerate(self.ensemble.populations):

            if cycles:
                phenotypes = list(self.phenotypes[i])
            else:
                phenotypes = population.phenotypes

            n = float(len(phenotypes))

            if n:
                if normed:
                    # maximum entropy possible in the population
                    # makes entropy independent of stability
                    normalize = log(n)

                unique_phenotypes = unique(phenotypes)
                n_unique = len(unique_phenotypes)
                self.unique.append(n_unique)

                # when no selection (s = inf and min_fit = 1 or random)
                # then there's no optimum. choose one randomly
                if random or population.random:
                    target = unique_phenotypes[randint(n_unique)]
                else:
                    target = population.opt_decimal

                equal = phenotypes.count(target)
                # don't count twice if target was chosen from phenotypes
                # or full enumeration includes pop.ini_decimal
                if random or population.random or full_enum:
                    equal -= 1
                equal /= n
                population.equal = equal
                self.equal.append(equal)

                entropy = 0.
                for phenotype in unique_phenotypes:
                    frequency = phenotypes.count(phenotype) / n
                    entropy  += frequency * log(frequency)

                if entropy: entropy /= -normalize
                population.entropy = entropy
                self.entropy.append( entropy)

        # for evolutionary records
        '''if not any(self.equals): self.equals = []
        self.equals.append(self.equal)

        if not any(self.entropies): self.entropies = []
        self.entropies.append(self.entropy)

        if not any(self.uniques): self.uniques = []
        self.uniques.append(self.unique)
        '''

    def pixel(self):

        self.period = []
        self.length = []
        self.phenotypes = []
        if not self.orbits: self.set_unique_orbits()
        # get arbitrary individual of the ensemble
        test_ind = self.ensemble.populations[0].individuals[0]
        max_fixpoint = pow(test_ind.base.size, test_ind.size) # pow(k, N)

        for population in self.ensemble.populations:
            for individual in population.individuals:
                self.length.append(individual.length)
                self.period.append(individual.period)
                if individual.period > 1:
                    if any(individual.orbit):
                        orbit = individual.orbit
                    else:
                        orbit = get_orbit(
                            individual.period, individual.genotype,
                            individual.phe_decimal,
                            individual.basename, globals()[individual.gpmap],
                            individual.base, individual.size)
                    phe_decimal = max_fixpoint + self.orbits.index(orbit)
                else:
                    phe_decimal = individual.phe_decimal
                self.phenotypes.append(phe_decimal)

        self.period = reshape(self.period, (self.ensemble.size, population.size))
        self.length = reshape(self.length, (self.ensemble.size, population.size))
        self.phenotypes = reshape(self.phenotypes, (self.ensemble.size, population.size))

    def set_unique_orbits(self):
        orbits = []
        for population in self.ensemble.populations:
            for individual in population.individuals:
                if individual.period > 1:
                    individual.orbit = get_orbit(
                        individual.period, individual.genotype,
                        individual.phe_decimal,
                        individual.basename, globals()[individual.gpmap],
                        individual.base, individual.size)
                    orbits.append(individual.orbit)
        # list to be indexable (sets are not: but set gets unique)
        self.orbits = list(set(orbits))
        # shuffle to destroy orbits' auto-correlation generated by looping
        shuffle(self.orbits)
        del orbits

    def table(self, initials, distance_function = int_vector_hamming_distance):

        # arbitrary individual of the ensemble
        size = self.ensemble.populations[0].individuals[0].size # == N
        self.stable      = []
        self.equal       = []
        self.equal_fix   = []
        self.equal_cycle = []
        self.hamming     = []
        self.gain        = zeros((size+1, size+1))#, dtype = int)
        for population in self.ensemble.populations:

            # get Target (wild-type phenotype)
            gpmap = population.gpmap
            if type(gpmap) == str:
                gpmap = globals()[gpmap]

            period, length, target_decimal = brent(
                devo, population.genotype, population.ini_decimal,
                population.basename, gpmap, population.base, size)

            if period > 1:
                target_orbit = get_orbit(
                    period, population.genotype, target_decimal,
                    population.basename, gpmap, population.base, size)

            stable      = 0
            equal_fix   = 0
            equal_cycle = 0

            for individual in population.individuals:

                # target ("wild-type") and perturbed ("mutant") phenotypes are Fixpoints
                if period == individual.period == 1:

                    # number of Fixed points produced by fixpoint wildtypes
                    stable += 1

                    # Equal
                    if individual.phe_decimal == target_decimal: equal_fix += 1

                    # Hamming distances
                    target_hamming = distance_function(
                        get_vector(
                            individual.phe_decimal, individual.basename, individual.base),
                        get_vector(
                            target_decimal, individual.basename, individual.base))

                    self.hamming.append(target_hamming)

                    # perturbation is in initial state
                    if initials:
                        source_hamming = distance_function(
                            get_vector(
                                individual.ini_decimal, individual.basename, individual.base),
                            get_vector(
                                population.ini_decimal, individual.basename, individual.base))

                    # perturbation is mutation
                    else:
                        source_hamming = genotype_distance(
                            individual.genotype, population.genotype)

                    # Gain function
                    self.gain[source_hamming, target_hamming] += 1

                # target ("wild-type") and perturbed ("mutant") phenotypes are Cycles
                elif period == individual.period > 1:

                    # Equal:
                    if individual.orbit == target_orbit: equal_cycle += 1

            # Normalize counts by population size
            self.stable.append(stable)
            self.equal.append(equal_fix + equal_cycle)
            self.equal_fix.append(equal_fix)
            self.equal_cycle.append(equal_cycle)


    def set_phenotypes(self):

        if not self.orbits: self.set_unique_orbits()
        # get arbitrary individual of the ensemble
        test_ind = self.ensemble.populations[0].individuals[0]
        max_fixpoint = pow(test_ind.base.size, test_ind.size) # pow(k, N)

        for population in self.ensemble.populations:
            population.phenotypes       = []
            population.phenotypes_fix   = []
            population.phenotypes_cycle = []
            for individual in population.individuals:
                if individual.period > 1:
                    phe_decimal = max_fixpoint + self.orbits.index(individual.orbit)
                    population.phenotypes_cycle.append(phe_decimal)
                else:
                    phe_decimal = individual.phe_decimal
                    population.phenotypes_fix.append(phe_decimal)
                population.phenotypes.append(phe_decimal)


    def evaluate_diversity_entropy(self, normalize = True):

        if not self.orbits:         self.set_unique_orbits()
        # get arbitrary population of the ensemble
        test_pop = self.ensemble.populations[0]
        if not test_pop.phenotypes: self.set_phenotypes()

        self.entropy       = []
        self.entropy_fix   = []
        self.entropy_cycle = []

        self.unique        = []
        self.unique_fix    = []
        self.unique_cycle  = []

        for population in self.ensemble.populations:

            self.entropy.append(
                evaluate_entropy(
                    population.phenotypes, normalize))

            self.entropy_fix.append(
                evaluate_entropy(
                    population.phenotypes_fix, normalize))

            self.entropy_cycle.append(
                evaluate_entropy(
                    population.phenotypes_cycle, normalize))

            self.unique.append(
                unique(population.phenotypes).size)

            self.unique_fix.append(
                unique(population.phenotypes_fix).size)

            self.unique_cycle.append(
                unique(population.phenotypes_cycle).size)


    def evaluate_robustness(self):

        self.equal     = []
        self.equal_fix = []
        self.intersect = []
        for population in self.ensemble.populations:

            # get wildtype phenotype
            gpmap = population.gpmap

            if type(gpmap) == str:
                gpmap = globals()[gpmap]

            period, length, population.phe_decimal = brent(
                devo, population.genotype, population.ini_decimal,
                population.basename, gpmap, population.base, population.bits)

            # set as the uniform output
            if period > 1:
                population.orbit = set(
                    get_orbit(
                        period, population.genotype, population.phe_decimal,
                        population.basename, gpmap, population.base,
                        population.bits))

            else:
                population.orbit = set((population.phe_decimal,))

            equal     = 0
            equal_fix = 0
            intersect = 0
            for mutant in population.individuals:

                if mutant.period > 1:
                    mutant.orbit = set(
                        get_orbit(
                            mutant.period, mutant.genotype, mutant.phe_decimal,
                            mutant.basename, globals()[mutant.gpmap], mutant.base,
                            mutant.size))

                else:
                    mutant.orbit = set((mutant.phe_decimal,))

                if population.orbit == mutant.orbit:
                    equal += 1

                elif population.orbit.intersection(mutant.orbit):
                    intersect += 1

                # old code/legacy: wildtype and mutant phenotypes are equal fixpoints
                if period == mutant.period == 1 and population.phe_decimal == mutant.phe_decimal:
                    equal_fix += 1

                del mutant.orbit

            del population.orbit

            # Store
            self.equal.append(equal)
            self.equal_fix.append(equal_fix)
            self.intersect.append(intersect)


    def evaluate_robustness_matrix(self, pop):
        for individual, population in zip(pop.individuals, self.ensemble.populations):
            individual.robustness = get_ind_robustness_matrix(
                individual, population.individuals)

        pop.robustness = array([
            individual.robustness.mean()
            for individual in pop.individuals])

    def evaluate_survivability_matrix(self, pop):
        for individual, population in zip(pop.individuals, self.ensemble.populations):
            individual.survivability = get_ind_survivability_matrix(
                population.individuals)

        pop.survivability = array([
            individual.survivability.mean()
            for individual in pop.individuals])



## Aliases
individual = Individual
population = Population


###############
## ROBUSTNESS #
###############

def evaluate_robustness(pop, initials=None, mutations=None):

    if pop and pop.size:
        ensemble = Ensemble(
            population=pop, initials=initials, mutations=mutations)
        ensemble.cycle_develop()
        analysis = EnsembleAnalysis(ensemble)
        analysis.evaluate_robustness()
        del pop, ensemble, analysis.ensemble
        return analysis
    else: return None

def get_ind_robustness_matrix(individual, mutants):
    return array([(individual.phe_decimal == mutant.phe_decimal)
                  *individual.period*mutant.period == 1
                  for mutant in mutants])

def get_ind_survivability_matrix(mutants):
    return array([mutant.period == 1 for mutant in mutants])

def evaluate_robustness_and_survivability_matrix(pop, stat):

    if pop and pop.size:
        ensemble = Ensemble(population=pop, mutations=True, all_mutants=True)
        ensemble.cycle_develop()
        analysis = EnsembleAnalysis(ensemble)
        if stat == 'robustness':
            analysis.evaluate_robustness_matrix(pop)
        if stat == 'survivability':
            analysis.evaluate_survivability_matrix(pop)
        if stat == 'robustness_and_survivability':
            ensemble.evaluate_robustness_and_survivability_matrix(pop)

def evaluate_robustness_and_survivability_matrix2(
        pop, mutations=True, all_mutants=True, n_neighbours=0, gpmap=None):

    if not mutations:
        prefix = 'initials_'
        if n_neighbours:
            prefix += '%d_' %n_neighbours
        else:
            prefix += 'all_'
    else:
        prefix = ''

    if pop and pop.size:

        # ensemble is mutant
        ensemble = Ensemble(
            population=pop, mutations=mutations, all_mutants=all_mutants,
            n_neighbours=n_neighbours)
        ensemble.cycle_develop(gpmap)

        # pop is wildtype
        if pop.sel_period == 1:
            ensemble.set_ind_robustness_and_survivability_matrix(
                pop, prefix, gpmap)
            pop.set_survivability_mean(prefix)
            pop.set_robustness_mean(prefix)

        else:
            ensemble.set_ind_survivability_matrix(pop, prefix)
            pop.set_survivability_mean(prefix)

def evaluate_survival(pop):

    if pop and pop.size:
        ensemble = Ensemble(population=pop, mutations=True, all_mutants=True)
        ensemble.cycle_develop()
        return array([
            get_ind_stat(population, 's')
            for population in ensemble.populations])
    else:
        return None #[NaN]*pop_size*bits**2



####################
## AUTO-REGULATION #
####################

# x is an instance of src.classes.individual
ind_functions = {'a':
                 lambda x: x.get_activating_fraction(),

                 'p':
                 lambda x: x.get_activating_fraction_diagonal(),

                 '1-p':
                 lambda x: x.get_activating_fraction_non_diagonal(),

                 'P':
                 lambda x: x.get_number_of_positive_diagonal(),

                 'q':
                 lambda x: x.get_number_of_positive_non_diagonal(),

                 'qp':
                 lambda x:(
                     x.get_number_of_positive_diagonal(),
                     x.get_number_of_positive_non_diagonal()),

                 'np':
                 lambda x: nonzero(diag(x.genotype))[0].size,

                 'nq':
                 lambda x: nonzero(nondiag(x.genotype))[0].size,

                 'sign':
                 lambda x: x.genotype.prod(),

                 's':
                 lambda x: x.period == 1,

                 'l':
                 lambda x: x.period,

                 'f':
                 lambda x: x.fitness,

                 'phe_decimal':
                 lambda x: x.phe_decimal,

                 'phenotype':
                 lambda x: x.phenotype,

                 'c':
                 lambda x: x.get_connectivity().mean(),

                 'kin':
                 lambda x: get_in_degrees(x.genotype),

                 'kout':
                 lambda x: get_out_degrees(x.genotype),

                 'kin+':
                 lambda x: get_in_number_of_positive(x.genotype),

                 'kout+':
                 lambda x: get_out_number_of_positive(x.genotype),

                 'initials_1_robustness.matrix':
                 lambda x: x.initials_1_robustness,

                 'initials_1_survivability.matrix':
                 lambda x: x.initials_1_survivability,

                 'initials_all_robustness.matrix':
                 lambda x: x.initials_all_robustness,

                 'initials_all_survivability.matrix':
                 lambda x: x.initials_all_survivability,

                 'robustness.matrix':
                 lambda x: x.robustness,

                 'survivability.matrix':
                 lambda x: x.survivability,

                 'robustness.diff':
                 lambda x: difference_diag_nondiag(x.robustness, mean),

                 'survivability.diff':
                 lambda x: difference_diag_nondiag(x.survivability, mean),

                 'q+clusters':
                 lambda x: bincount(
                     label_clusters(
                         nondiag(
                             x.genotype),-1)[0].flatten())[1:],

                 'q-clusters':
                 lambda x: bincount(
                     label_clusters(
                         nondiag(
                             x.genotype), 1)[0].flatten())[1:],

                 'a+clusters':
                 lambda x: bincount(
                     label_clusters(
                         x.genotype,-1)[0].flatten())[1:],

                 'a-clusters':
                 lambda x: bincount(
                     label_clusters(
                         x.genotype, 1)[0].flatten())[1:],

                 'n+clusters':
                 lambda x: [label_clusters(x.genotype, -1)[1]],

                 'n-clusters':
                 lambda x: [label_clusters(x.genotype,  1)[1]],

                 'q+slices':
                 lambda x: find_objects(nondiag(x.genotype), -1),

                 'q-slices':
                 lambda x: find_objects(nondiag(x.genotype),  1),

                 'q+q':
                 lambda x: bincount(
                     label_clusters(
                         nondiag(x.genotype), -1)[0].flatten())[1:] / float(
                             get_number_of_positive_non_diagonal( x.genotype))}

# from itertools import chain; list(chain.from_iterable(q))
def get_ind_stat(pop, stat, pop_size=P, bits=N, **args):

    if pop and pop.size:
        return array([ind_functions[stat](individual)
                      for individual in pop.individuals])

    elif 'matrix' in stat:
        return array([NaN]*pop_size*bits**2).reshape((pop_size, -1))

    elif stat == 'qp':
        return array([NaN]*pop_size*2).reshape((pop_size, -1))

    elif stat == 'phenotype':
        return array([NaN]*pop_size*bits).reshape((pop_size, -1))

    else:
        return array([NaN]*pop_size)

def get_pop_stat(pop, stat, pop_size=1000, experiment='', p=p, **args):

    if  'p_vs_' in experiment:
        return ([p]*pop_size, getattr(pop, stat))

    if not (pop and pop.size):
        return array([NaN]*pop_size)

    if 'qp.pop' in experiment or experiment == 'autoregulation':
        if hasattr(pop, stat):
            return getattr(pop, stat)
        else:
            return array([NaN]*pop_size)

    else:
        return [getattr(individual, stat) for individual in pop.individuals]

def autoregulation(pop, stats, mutations=n_perturb, initials=None):
    pop.cycle_develop() # for period in stats
    for stat in stats:
        setattr(pop, stat, get_ind_stat(pop, stat, pop.size))
    analysis = evaluate_robustness(pop, initials, mutations)
    pop.equal, pop.equal_fix, pop.intersect = array(
        analysis.equal), array(analysis.equal_fix), array(analysis.intersect)
    del analysis, pop.individuals, pop.optimum, pop.base
    #return pop


## Auxiliary functions

# Generate and analyse perturbed ensembles
def perturb_and_analyse(pop, data, initials=None, mutations=None):

    ## Pertubations: mutations or initial conditions
    ensemble = Ensemble(population=pop, initials=initials, mutations=mutations)
    ensemble.cycle_develop()

    ## Analysis
    analysis = EnsembleAnalysis(ensemble)
    # entropy and set both individual and unique orbits
    analysis.evaluate_entropy_and_set_orbits()
    # stability, robustness, hamming distance and gain function (if initials)
    analysis.table(any(initials))

    ## Append analysed data
    for stat, name in zip([analysis.stable, analysis.equal, analysis.hamming, analysis.entropy, analysis.equal_fix, analysis.equal_cycle, analysis.entropy_fix, analysis.entropy_cycle], ['stable', 'equal', 'hamming', 'entropy', 'equal fix', 'equal cycle', 'entropy fix', 'entropy cycle']):
        data[name].append(stat)
    if any(initials): data['gain'] += analysis.gain

    del ensemble, analysis
    #return data

def perturb_and_analyse2(pop, initials = None, mutations = None, normalize = True,
                         distance_function = int_vector_hamming_distance, n_bins = 20):

    if not pop: return {}, 0, 0

    ## Pertubations: initial conditions or mutations
    # if initials:  each individual of the original population (pop) generates an isogenic population with different initial conditions
    # if mutations: each individual of the original population (pop) generates a population of single mutants with the same initial state
    ensemble = Ensemble(population=pop, initials=initials, mutations=mutations)
    ensemble.cycle_develop()

    ## Analysis
    analysis = EnsembleAnalysis(ensemble)
    # evaluate number of unique phenotypes and entropy for each population
    # it also sets both individual and unique orbits and lists of phenotypes
    analysis.evaluate_diversity_entropy(normalize)
    # stability, robustness, hamming distance and gain function
    analysis.table(any(initials), distance_function)

    ## Epistasis
    all_epistasis    = evaluate_epistasis(pop, distance_function = distance_function)
    # mask Inf and NaNs
    all_epistasis    = ma.masked_invalid(all_epistasis)
    epistasis        = Epistasis()
    epistasis.mean   = all_epistasis.mean()
    epistasis.std    = all_epistasis.std()
    epistasis.hist   = histogram(all_epistasis.compressed(), n_bins)
    epistasis.unique = unique(all_epistasis.compressed()), array(unique_hist(all_epistasis.compressed()))

    ## Append analysed data
    data = dict(zip(pa_stats,
                    (analysis.stable, analysis.equal, analysis.equal_fix, analysis.equal_cycle, analysis.unique, analysis.unique_fix,
                     analysis.unique_cycle, analysis.entropy, analysis.entropy_fix, analysis.entropy_cycle, analysis.hamming, epistasis)))

    return data, ensemble, analysis



##############
## EPISTASIS #
##############

class Epistasis():
    def __repr__(self):
        return str(self.__dict__)

# individual is class instance; mutant is genotype array
def get_hamming(individual, mutant, distance_function = int_vector_hamming_distance):
    # get mutant phenotype
    period, length, phe_decimal = test_stability(mutant, individual.ini_decimal, individual.basename,
                                                 globals()[individual.gpmap], individual.base, individual.size)
    # if mutant is still a fixed point
    if period == 1:
        return distance_function(get_vector(phe_decimal,            individual.basename, individual.base, individual.size),
                                 get_vector(individual.phe_decimal, individual.basename, individual.base, individual.size))
    else:
        return inf

def evaluate_epistasis(pop, n_mutations = None, change_sign = False, deletion = False,
                       distance_function = int_vector_hamming_distance):
    if pop and pop.size:
        # init
        if pop.binary: change_sign = True
        # get number of non-zero mutational targets
        size = n_mutations if n_mutations else rint(pop.individuals[0].genotype.size*pop.density).astype(int)
        # number of epistatic interactions (i.e. pairs of nodes or edges)
        n = misc.comb(size, 2, 1)
        epistasis = []
        genes     = []

        # for all individuals in the population
        for individual in pop.individuals:

            # only fixed-point wildtypes
            if individual.period == 1:

                # only mutate non-zero elements, randomly
                mutating = get_mutating(individual.genotype, n_mutations)

                # first single mutant; gene_i is index of non-zero element of flattened genotype
                for i, gene_i in enumerate(mutating):
                    mutant_i  = get_single_mutant(individual.genotype, gene_i, change_sign, deletion)
                    hamming_i = get_hamming(individual, mutant_i, distance_function)

                    # second single mutant
                    for gene_j in mutating[i+1:]:
                        mutant_j  = get_single_mutant(individual.genotype, gene_j, change_sign, deletion)
                        hamming_j = get_hamming(individual, mutant_j, distance_function)

                        # double mutant
                        mutant_ij  = get_double_mutant(mutant_i, mutant_j, gene_j)
                        hamming_ij = get_hamming(individual, mutant_ij, distance_function)

                        epistasis.append(get_epistasis(hamming_i, hamming_j, hamming_ij, individual.size))
                        genes    .append((gene_i, gene_j))

                if mutating.size < size:
                    n_nan = (n-misc.comb(mutating.size, 2, 1))
                    epistasis.extend([NaN]*n_nan)
                    genes    .extend([NaN]*n_nan)

            # wildtpe is a cycle
            else:
                epistasis.extend([NaN]*n)
                genes    .extend([NaN]*n)

        #if len(epistasis) == pop.size*n: epistasis = reshape(epistasis, (pop.size, n))
        return epistasis, genes
    else: return [], []

#    get_epistasis(mutant_i, mutant_j, get_double_mutant(mutant_i, mutant_j, gene_j), individual) for gene_j in mutations[i+1:])

#
def filter_ensemble_populations(
        model_number, experiment='initials', start=0, end=100, threshold=1):
    filename = ensembles_dir + 'model%s_%s_ensemble' %(model_number, experiment)

    for i in range(start, end):
        file = open(filename + '%s.dat' %i)
        ensemb = cPickle.load(file)
        file.close()

        new_ensemb = []
        for population in ensemb.populations:
            new_pop = []
            for individual in population.individuals:
                if int_vector_hamming_distance(
                        get_vector(
                            individual.ini_decimal, individual.basename,
                            individual.base),
                        get_vector(population.ini_decimal, individual.basename,
                                   individual.base)) <= threshold:
                    new_pop.append(individual)
            population.individuals = new_pop
            population.size = len(new_pop)
            new_ensemb.append(population)
        ensemb.populations = new_ensemb

        file = open(filename + '%s-filtered.dat' %i, 'wb')
        cPickle.dump(ensemb,file)
        file.close()

        del ensemb, new_ensemb, population, new_pop

def load_noise_interval(filename, develop=False, save=False):
    file   = open(filename + '.dat', 'rb')
    ensemb = cPickle.load(file)
    file.close()

    if develop: ensemb.develop()
    if save:
        file = open(filename + '-dev_aft.dat', 'wb')
        cPickle.dump(ensemb, file)
        file.close()

    analysis = EnsembleAnalysis(ensemb)
    analysis.evaluate_stability()
    stable_mean = mean(analysis.stable)
    stable_std  = std( analysis.stable)

    # old saved ensembles don't have decimal attributes
    '''analysis.compare_phenotypes()
    stable_mean = mean(analysis.stable)
    stable_std  = std( analysis.stable)
    '''
    del ensemb, analysis
    return stable_mean, stable_std

# noises = append(array([0.05,0.1,0.15]), arange(0.2,1.5,0.1))
def reconstruct_noise_interval(filename, noises, develop=False, save=False):
    stable_mean = []; stable_std = []
    for i in noises:
        file = filename + str(i)
        y1, y2 = load_noise_interval(file, develop, save)
        stable_mean.append(y1); stable_std.append(y2)
    #fig = figure()
    #errorbar(x, stable_mean, stable_std)
    return noises, stable_mean, stable_std #fig

def stability_fitness_and_hamming_evolution(pop, generations = G, sel_strength = s, min_fit = zero_fit, sel_period = 1, basename = def_basename, gpmap = def_gpmap, devo_function = devo, rec_rate = r, mut_rate = u, ins_rate = uA, del_rate = uD, keys = ['stable', 'fitness', 'hamming', 'size', 'M2', 'connectivity', 'densities', 'c']):

    data = {}

    # generation 0
    pop.cycle_develop(gpmap, devo_function)
    pop.evaluate_fitness()
    pop.evaluate_stability_fitness_hamming_and_M1()
    for key, stat in zip(keys, [pop.stable, pop.fitness, pop.hamming, pop.size, 1-pop.M1, pop.connectivity, pop.densities, pop.c]):
        data[key] = [stat]

    # Evolution
    pop.sel_strength = s
    for i in range(1, generations + 1):
        pop.evolution(True, gpmap, devo_function, computer_other_properties = True)
        for key, stat in zip(keys, [pop.stable, pop.fitness, pop.hamming, pop.size, 1-pop.M1, pop.connectivity, pop.densities, pop.c]):
            data[key].append(stat)

        # if population has gone extint
        if not pop.size: break

    return data

def elhanan(stable):

    while True:
        wildtype = Individual(random = True)
        if wildtype.stable == stable: break

    if wildtype.period > 1: wildtype.phenotype = wildtype.orbit
    else:                   wildtype.phenotype = wildtype.phe_decimal

    ## Perturbed initials (in decimal code) are nearest neighbours
    initials = get_nearest_neighbours(wildtype.ini_decimal, wildtype.size)
    init_same_stability = 0
    init_same_phenotype = 0
    mut_same_stability  = 0
    mut_same_phenotype  = 0

    for ini_decimal in initials:

        ## Perturbing initial conditions or ER
        mutant_period, length, mutant_phe_decimal = brent(devo, wildtype.genotype, ini_decimal)
        if mutant_period > 1: mutant_phenotype = get_orbit(mutant_period, wildtype.genotype, mutant_phe_decimal)
        else:                 mutant_phenotype = mutant_phe_decimal

        a, b = compare_stability_phenotype(wildtype.period, mutant_period, wildtype.phenotype, mutant_phenotype)
        init_same_stability += a;  init_same_phenotype += b

        ## Perturbing genotype or GR
        mutant_genotype = mutate(wildtype.genotype, 1)
        mutant_period, length, mutant_phe_decimal = brent(devo, mutant_genotype, wildtype.ini_decimal)
        if mutant_period > 1: mutant_phenotype = get_orbit(mutant_period, mutant_genotype, mutant_phe_decimal)
        else:                 mutant_phenotype = mutant_phe_decimal

        a, b = compare_stability_phenotype(wildtype.period, mutant_period, wildtype.phenotype, mutant_phenotype)
        mut_same_stability += a; mut_same_phenotype += b

    n = len(initials)
    mut_diff_stability = n - mut_same_stability
    mut_diff_phenotype = mut_same_stability - mut_same_phenotype
    init_diff_stability = n - init_same_stability
    init_diff_phenotype = init_same_stability - init_same_phenotype

    return mut_diff_phenotype, mut_diff_stability, init_diff_phenotype, init_diff_stability # GR_x, GR_y, ER_x, ER_y


## for profiling purposes only
def main():
    '''file = open(pop_dir + 'individual.dat', 'rb')
    ind = cPickle.load(file)
    file.close()
    population  = [ind for i in range(10000)]

    new = []
    for ind in population:
        x = ind.develop()
        new.append(x)
   
    new = [ ind.develop(cSigmoid2) for ind in population ]
    
    file = open(ensembles_dir + 'ensemble2.dat', 'rb')
    ensemb = cPickle.load(file)
    file.close()
    ensemb.develop(gpmap = sign)
    analysis = EnsembleAnalysis(ensemb)
    analysis.binary_to_decimal()
    
    return Ensemble(founder = generate_genotype())
    
    file   = open(ensembles_dir + 'ensemble.dat', 'rb')
    ensemb = cPickle.load(file)
    file.close()
    ensemb.cycle_develop()
    #ensemb.develop(start_over = True, gpmap = sign)
    '''
    #pop = Population(size = 5000, random = True) # random population

    '''filename = 'random_pop5000.dat' #'population.dat'
    file     = open(pop_dir + filename , 'rb')
    pop      = cPickle.load(file)
    file.close()
    
    #[pop.evolution() for i in range(3)]
    
    #gpmap = sign0
    dim =  2
    min = -1.
    max =  1.
    jumps, base = get_equalprobable_states(dim, min, max)
    basename    = array2string(base, precision = 2)
    step        = get_step_function(dim, 'cstep')
    #def gpmap(x): return step(x, jumps, base)
    #gpmap = csign
    gpmap = sign
    #gpmap = sign0
    [pop.cycle_develop(gpmap, devo) for i in range(4)]
    '''
    #stability_fitness_and_hamming_evolution(50)

    '''import runs.run_stability_vs_c_evolved
    reload(runs.run_stability_vs_c_evolved)
    '''
    #import runs.run_insertion_deletion2

    '''binary = True
    gpmap = get_gpmap(binary = binary)
    pop = Population(random = True, binary = binary, gpmap = gpmap)
    generations = 2e2
    [pop.evolution(change_sign = binary) for g in range(1, int(generations) + 1)]
    return pop
    '''
'''
import cProfile
import pstats
filename = data_dir + 'profiling/classes.pstats'
cProfile.run('import src.classes; src.classes.main()', filename)
stats = pstats.Stats(filename).sort_stats('cumulative').print_stats(15)

cProfile.run('reload(src.classes); src.classes.main()', filename)

fig.canvas.draw() 

bins = arange(0,1.1,0.1)
label_hist(fig, ['% stable runs', '% stable runs equal to deterministic', 'entropy of stable runs'], 3, )
ylabel('# genomes')
xlabel('0 <= x <=1 (x see legend)')
title('deterministic stable ensemble 300x100x10\nsampling initial phenotype')

label_hist(fig, ['sampling initial', 'noise = 0.05', 'noise = 0.15', 'det. after noise = 0.05', 'det. after noise = 0.15'], 5, code=0)

get_printoptions()
set_printoptions(threshold=nan)
'''
