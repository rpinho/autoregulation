def get_model_label_mut(
        shortname, sel_strength='\sigma', sel_period='l',
        inf='\infty', space_char=',\/', mut_rate='\mu'):

    if shortname == 'm4r05b':
        return r'$%s=%s%s%s=%g$' %(sel_strength, s,   space_char, mut_rate,   u)

    if shortname == 'm4r05u1':
        return r'$%s=%s%s%s=%g$' %(sel_strength, s,   space_char, mut_rate,   1)

    if shortname == 'm4r05u001':
        return r'$%s=%s%s%s=%g$' %(sel_strength, s,   space_char, mut_rate, .01)

    # default
    return shortname

def get_model_label(
        shortname, sel_strength='\sigma', sel_period='l',
        inf='\infty', space_char=',\/', mut_rate='\mu'):

    if shortname == 'm4r05p0':
        return r'$%s=%s$'        %(sel_strength, s)

    if shortname == 'm4r05b':
        return r'$%s=%s%s%s=%d$' %(sel_strength, s,   space_char, sel_period, 1)

    if shortname == 'm3r05p2':
        return r'$%s=%s%s%s=%d$' %(sel_strength, inf, space_char, sel_period, 2)

    if shortname == 'm3r05p3':
        return r'$%s=%s%s%s=%d$' %(sel_strength, inf, space_char, sel_period, 3)

    if shortname == 'm3r05p4':
        return r'$%s=%s%s%s=%d$' %(sel_strength, inf, space_char, sel_period, 4)

    if shortname == 'm3r05p5':
        return r'$%s=%s%s%s=%d$' %(sel_strength, inf, space_char, sel_period, 5)

    if shortname == 'm3r05p6':
        return r'$%s=%s%s%s=%d$' %(sel_strength, inf, space_char, sel_period, 6)

    if shortname == 'm3r05p7':
        return r'$%s=%s%s%s=%d$' %(sel_strength, inf, space_char, sel_period, 7)

    if shortname == 'm3r05b':
        return r'$%s=%s%s%s=%d$' %(sel_strength, inf, space_char, sel_period, 1)

    if shortname == 'm2r05b':
        return r'$%s=%s$'        %(sel_strength, inf)

    if shortname == 'b05p05q05':
        return r'$p_0=%.1f %s q_0=%.1f$' %(.5, space_char, .5)

    if shortname == 'm3p05q09':
        return r'$p_0=%.1f %s q_0=%.1f$' %(.5, space_char, .9)

    if shortname == 'm3p05q01':
        return r'$p_0=%.1f %s q_0=%.1f$' %(.5, space_char, .1)

    if shortname == 'm3p09q09':
        return r'$p_0=%.1f %s q_0=%.1f$' %(.9, space_char, .9)

    if shortname == 'm3p09q01':
        return r'$p_0=%.1f %s q_0=%.1f$' %(.9, space_char, .1)

    # default
    return shortname

def get_model_label_mut_bias(shortname, mut_bias='neutrality = '):

    if shortname == 'b09p05q05':
        return '%s%.1f' %(mut_bias, .9)

    if shortname == 'b08p05q05':
        return '%s%.1f' %(mut_bias, .8)

    if shortname == 'b07p05q05':
        return '%s%.1f' %(mut_bias, .7)

    if shortname == 'b06p05q05':
        return '%s%.1f' %(mut_bias, .6)

    if shortname == 'b05p05q05':
        return '%s%.1f' %(mut_bias, .5)

    if shortname == 'b04p05q05':
        return '%s%.1f' %(mut_bias, .4)

    if shortname == 'b03p05q05':
        return '%s%.1f' %(mut_bias, .3)

    if shortname == 'b02p05q05':
        return '%s%.1f' %(mut_bias, .2)

    if shortname == 'b01p05q05':
        return '%s%.1f' %(mut_bias, .1)

    if shortname == 'b005p05q05':
        return '%s%.2f' %(mut_bias, .05)

def get_model_label_init_var(shortname):

    if shortname == 'm3r05b':
        return '$var(q_0) < var(p_0)$'

    if shortname == 'b05p05q05':
        return '$var(q_0) = var(p_0) = 0$'

    if shortname == 'm3r05d':
        return '$var(q_0) = var(p_0) > 0$'

    return shortname

def get_model_label_mut_text(shortname):

    if 'r05b' in shortname:
        return '$\mu$'

    if 'r05u1' in shortname:
        return 'higher mutation rate'

    if 'r05u001' in shortname:
        return 'lower mutation rate'

    return shortname

def get_model_label_sel_text(shortname):

    if 'r05b' in shortname:
        return '$\sigma$'

    if 'r05s1' in shortname:
        return 'weaker selection'

    if 'r05s001' in shortname:
        return 'stronger selection'

    return shortname

def get_model_label_rec_text(shortname):

    if 'r05b' in shortname:
        return 'recombination'

    if 'r0' in shortname:
        return 'no recombination'

    return shortname

def get_model_label_density_text(shortname):

    if 'r05b' in shortname:
        return ' dense, recombination'

    if 'r05k2' in shortname:
        return 'sparse, recombination'

    if 'r0k2' in shortname:
        return 'sparse, no recombination'

    return shortname

def get_model_label_text(shortname):

    if shortname == 'm4r05b':
        return get_model_label_text('m3r05b') + ' and target'

    if shortname == 'm3r05b':
        return 'selecting for stability'

    if shortname == 'm3r05p2' or shortname == 'm3r05p3':
        return 'selecting against stability'

    if shortname == 'm4r05p0':
        return 'neutral for stability, selecting for target'

    if shortname == 'm2r05b':
        return 'neutral'

    if 'u' in shortname:
        return get_model_label_mut_text(shortname)

    if 's' in shortname:
        return get_model_label_sel_text(shortname)

    if shortname == 'm4r0b':
        return get_model_label_rec_text(shortname)

def get_model_label_short_text(shortname):

    if shortname == 'm4r05b':
        return 'target'

    if shortname == 'm3r05b':
        return 'no target'

    if shortname == 'm3r05p2':
        return 'cycles of length 2'

    if shortname == 'm3r05p3':
        return 'cycles of length 3'

    if shortname == 'm4r05p0':
        return 'target'

    if shortname == 'm2r05b':
        return 'no target'

    if shortname == 'm3r05m':
        return 'off is 0'

    if shortname == 'm3r05a1':
        return 'real vectors (a=1)'

def get_stat_label(stat, diversity = None):

    if stat == 'p':
        return 'p'

    if stat == 's':
        return 'stability'

    if stat == 'equal':
        return 'robustness'#$r^=$'

    if stat == 'G.std(axis=0)':
        return 'diversity'

    if stat == 'period':
        return 'attractor length'#r'$l^{opt}$'

    if stat == 'diversity' and diversity:
        return get_diversity_label(diversity.keys()[0])

    if stat == 'conservation':
        return 'abs(average matrix).mean'

    if stat == 'a.matrix':
        return 'a.matrix.mean'

    if stat == '1-p':
        return 'q'

    if type(stat) == list:
        return str(map(get_stat_label, stat))

    return stat

def get_stat_label_text(stat):

    if stat == 'p':
        return 'diagonal'

    if stat == '1-p':
        return 'non-diagonal'

def get_diversity_label(stat):

    if stat == 'G.std(axis=0)':
        return 'within population,  between individuals'

    if stat == 'G.std(axis=1)':
        return 'within individuals, between rows'

    if stat == 'G.std(axis=2)':
        return 'within individuals, between columns'

    return stat

def get_y_label(stats, stat_label=get_stat_label):
    return ', '.join(map(stat_label, stats))
    #str(stats) if len(stats) > 1 else stat

def get_myrun_matrices(myrun):

    if myrun.binary:
        weights = 'binary'
    else:
        weights = 'real'

    if myrun.full_enum:
        stability = 'full enumeration of'
    elif myrun.stable:
        stability = 'stable'
    else:
        stability = 'random'

    return ' '.join((stability, weights))

def get_noise_label(noise, noise_function='random_noise', noise_time='before'):
    #noises = [1e-100, .01, .05, .1, .15, .2, .25]
    #labels = [
     #   'no noise', 'noise=0.01', 'noise=0.05', 'noise=0.10', 'noise=0.15',
      #  'noise=0.2', 'noise=0.25']

    if noise_function == 'random_noise' or noise_time == 'shmulevich':
        if noise == 1e-100:
            return 'no noise'
        elif noise:
            return 'noise=%.2f' %noise

    else:
        #flips_label = ' '.join(noise_function.split('_')[-2:])
        if noise == 1e-100:
            return 'no flip'
        elif noise == 1:
            return '1 flip'
        else:
            return ' '.join((str(noise), 'flips'))

def get_myrun_label(myrun):
    if myrun.noise:
        return get_noise_label(
            myrun.noise, myrun.noise_function, myrun.noise_time)
    else:
        return get_myrun_matrices(myrun)

def get_myrun_title(myrun, experiment='dims'):

    matrices = get_myrun_matrices(myrun) + ' matrices'

    if 'dims' in experiment:
        matrices += ' of size N = %d' %myrun.bits

    if 'overrepresentation' in experiment:
        matrices += ', k = %d' %myrun.dim

    if myrun.degree < myrun.bits:
        matrices += ' and density %.1f' %myrun.density

    if myrun.min == 0:
        matrices += ', min = 0'

    noise_func = ''
    if myrun.noise_function != 'random_noise':
        noise_func = ''.join(('\n',
                              ' '.join(myrun.noise_function.split('_')[-2:])))

    if myrun.noise_time == 'after':
        noise_func = '\nnoise after matrix multiplication'

    return matrices + noise_func
