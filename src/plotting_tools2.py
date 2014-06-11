# mac
import sys
if sys.platform == 'darwin':
    from pylab import *
    #from mayavi import mlab
    from mpl_toolkits.mplot3d import Axes3D
    import scikits.bootstrap as bootstrap
    import pandas as pd
    import ggplot
    import colorbrewer
import numpy.ma as ma
from string import uppercase
import classes3
reload(classes3)
from classes3 import *

import shortnames
reload(shortnames)
from shortnames import *

# load rcParams for specific journal
import src.plos as plos

#########################################################
## PLOTTING: HISTOGRAMS AND IMAGES AUXILIARY FUNCTIONS ##
#########################################################

# there's a precision issue here: using a workaround (1+1e-15)
# while i don't figure it out
#c = 1 + 1e-15

## Stability and path length parameters
# filename_exceptions
density_precisions = {174:3, 300:3, 418:3, 559:3, 747:3, 1000:3, 2000:3, 3000:4,
                      5000:3, 7000:3, 10000:4} # if N > 1000 then precision > 2
density_precisions = {300:3, 1000:3, 2000:3, 3000:4, 5000:3, 7000:3, 10000:4}
max_bits = {'stability_vs_N': {4:300, 6:60, 1.0:100, c:100},
            'path_leng_vs_N': {2:10000, 4:300, 1.0:100, c:100},
            'viability_vs_N': {3:100}}

# path_length
bitsk2 = [4, 10, 20] + range(50,1001,50)
bitsk2 += list(get_int_datapoints(log10(752), log10(2500), 10))
bitsk2.sort()
path_bits = {c: array(range(4,46) + [47, 50]),
             4: array([4, 5] + range(10,101,5)),
             2: array(bitsk2)}

# regression
fit_functions = {c: [array]*2,
                 4: [array, log],
                 2: [log]*2}

# stability
bitsc1  = range(4,100)
bitsk4  = bitsc1 + range(100, 200, 5)
bitsk4 += list(get_int_datapoints(log10(111), log10(201), 5))
bitsk4 += list(get_int_datapoints(log10(114), log10(204), 5))
bitsk4 += [200, 206, 220]
bitsk4  = unique(bitsk4)#.sort()
bitsk2  = list(bitsk4) + [250] + range(300, 700+1, 100)
bitsk2 += range(850, 1300+1, 150) + [1500, 1700]
bitsk2 += range(2000, 5000+1, 500)
bitsk2 += range(3001, 5001+1, 500) # mean
bitsk2 += range(3002, 5002+1, 500) # mean
bitsk2 += range(3003, 4503+1, 500) # mean
bitsk2 += [6000, 7000, 8000, 8500, 10000] # mean
bitsk2 += list(get_int_datapoints(log10(1499), log10(5000), 10)) # maximum
bitsk2 += list(get_int_datapoints(log10(4), log10(1000), 10))
bitsk2 += list(get_int_datapoints(log10(4), log10(1000), 20))
bitsk2 += list(get_int_datapoints(log10(1000), log10(3000), 4))
bitsk2 += list(get_int_datapoints(log10(3000), log10(10000), 4))
bitsk2 = unique(bitsk2)#.sort()
stability_bits = {c: array(bitsc1),
                  4: array(bitsk4),
                  2: array(bitsk2)}
delete_X = [[3002, 3003, 3346], [
    111, 114, 128, 131, 149, 152, 155, 165, 170,
    173, 176, 180], [37, 38, 39, 41, 42, 43, 44, 45, 47]]

min_N = 0
min_Samples = 0
ecoli = {'N': 1621, 'c': .002, 'p': 28./116}
v_lines = {'ecoli': {'N': 1621, 'c': .002, 'p': 28./116},
           'yeast': {'N': 4410+157}}

## Viability parameters
bitsc1  = range(4, 30)
bitsk4  = range(4, 20) + range(20, 150, 5)
bitsk3  = list(get_int_datapoints(log10(5), log10(2000), 20))[:-5]
bitsk3 += list(get_int_datapoints(log10(5), log10(100),  12))
bitsk3 += [53, 73, 100]
bitsk3.sort()
bitsk2  = range(5, 600, 25)
bitsk2 += list(get_int_datapoints(log10(600.1), log10(2000), 10))
bitsk2 += list(get_int_datapoints(log10(5), log10(2000), 20))
bitsk2.sort()
bitsk1  =      get_int_datapoints(log10(5), log10(2000), 20)
viability_bits = {c: array( bitsc1),
                  4: array( bitsk4),
                  3: unique(bitsk3),
                  2: unique(bitsk2)}#,
                  #1: array( bitsk1)}

#mutations = [[1]*3 + [0], [False, True, False, False],
 #            [False, False, True, False]]
mutations = {'uuall' : (1, False, False), 'u0all' : (1, True, False),
             'u-all' : (1, False, True), 'uu100' : (0, False, False)}

# pop.evolution
delete_x_index = [42, 43, 45, 46, 48, 49, 52, 58, 59, 61, 62, 64, 65, 67, 74,
                  75, 77, 78, 79, 81, 83]
delete_x = [  1102,   1270,   1687,   1944,   2582,   2976,   4556,  10672,
              12299,  16335,  18825,  25002,  28813,  38267, 103309, 119057,
              158121, 182224, 210001, 278903, 370413]

pop_stats = ['p', '1-p', 's'] + ['l', 'f'] + ['diversity'] + ['conservation']
pop_stats += ['density']# + ['epistasis']
pop_diversity = {'n(G)':0 , 'S(G)':1, 'G.std(axis=0)':2, 'G.std(axis=1)':3,
                 'G.std(axis=2)':4, 'n(P)':5, 'S(P)':6, 'P.std':7}
pop_normalization = {'np': N, 'nq': N*(N-1), 'q': N*(N-1),
                     'n(G)': P, 'n(P)': P,
                     #'conservation': P,
                     'S(G)': -log(P), 'S(P)': -log(P),
                     'stable': n_perturb, 'unique': n_perturb,
                     'equal': n_perturb, 'equal_fix': n_perturb,
                     'intersect': n_perturb}
pop_errors = {}
error_types = ('within', 'between', 'mean')
pop_generations = get_pop_datapoints('m3r05b', 1e7)[1:]

## Perturb and analyse data parameters
gr_dir = data_dir + 'scale_free_graph/'
#stats = ['equal cycle', 'entropy cycle', 'equal', 'entropy', 'gain',
 #        'hamming', 'stable', 'equal fix', 'entropy fix']
extra_stat = None#'entropy fix normed']
G_string   = 1000

## Figure parameters
formats      = {'plos': 'pdf', None: 'pdf'} #eps
dpis         = {'plos': 300, 'bmc': 600, None: 300}
fontsizes    = {'plos': 12,  'bmc': 12,  None: 10} #plos: 8-12 point
fontnames    = {'plos': 'Times',#'Arial 'Symbol'],
               'bmc': '',  None: 'Bitstream Vera Sans'}
#fontweight = 'bold'
panel_labels = {'plos': uppercase, None: uppercase}
fig_prefixes = {'plos': '', None: ''}
fig_widths   = {'plos': [3.27, 6.83, 6.83], None: [8, 12, 16]}
fig_height   = 6.#9.19
# fig.get_size_inches(); fig.set_size_inches(w,h)
screen_size  = {'macbook': (16., 8.775), 'imac': (24., 13.775), None:None}
bbox         = 'tight'
pad_inches   = 0.1
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio

#rcParams['axes.color_cycle']
#color = ''
colors = {None: ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'darkorange', 'sienna',
                 'lime', 'magenta', 'indigo', 'goldenrod', 'gray', 'tomato'],
          'plos': ['k', 'r', 'b', 'g', 'darkorange']}
markers    = ['o', 's', 'v', 'd', '*', '^', 'x', '+', 'p', 'h']
linestyles = ['-', '--', '-.', ':']
span_linestyles = ['solid', 'dashed', 'dashdot', 'dotted']
linewidth = 1
marker_size = 8

# Customizing matplotlib: Dynamic rc settings
def set_custom_rcParams(journal, experiment, n_cols=1):

    # figure size in inches (depends on number of columns)
    figwidth = journal.figure.pop('figwidth', (8,12))[n_cols-1]
    aspct_ratio = journal.figure.pop('aspct_ratio', 3./4)
    figheight = journal.figure.pop('figheight', figwidth*aspct_ratio)
    journal.figure['figsize'] = (figwidth, figheight)

    # default directory in savefig dialog box
    journal.savefig['directory'] = get_figs_save_dir('', experiment)

    matplotlib.rc('lines',   **journal.lines)
    matplotlib.rc('patch',   **journal.patch)
    matplotlib.rc('font',    **journal.font)
    matplotlib.rc('axes',    **journal.axes)
    matplotlib.rc('grid',    **journal.grid)
    matplotlib.rc('legend',  **journal.legend)
    matplotlib.rc('figure',  **journal.figure)
    matplotlib.rc('savefig', **journal.savefig)

def reset_default_rcParams():
    matplotlib.rcdefaults()


'''def set_fig_params(experiment = None, journal = None, prefix = None, n_col = 0):
    global fig_prefix, fontsize, dpi, fig_width, save_dir, format, color
    fig_prefix = prefix
    fontsize  = fontsizes[journal]
    dpi        = dpis[journal]
    fig_width  = fig_widths[journal][n_col]

    if journal:
        save_dir   = report_dir + experiment + '/figs/'
        format     = 'eps'
        color      = 'k'

    # Messing around, miscellaneous figures
    else:
        save_dir   = figs_dir
        format     = 'pdf'
        color      = ''

    ## Poster
    if prefix == 'poster':
        linewidth = 2
        rc("axes",  lw=linewidth)
        rc("lines", lw=linewidth, mew=linewidth, markersize=marker_size)
        rc("font",  size=fontsize)

    #else:
    #    rcdefaults()

set_fig_params(experiment, journal, prefix, n_col)
'''

########################
## AUXILIARY FUNCTIONS #
########################

def get_norm(stat, shortname=None):

    if stat in pop_normalization.keys():
        return pop_normalization[stat]

    #elif stat == 'p' and 'k2' in shortname: return 2
    else:
        return 1

def delete_data_points(fig, delete_x):
    axes = fig.get_axes()[0]
    for i, points in enumerate(delete_x):
        line  = axes.get_children()[2 + i]
        index = [where(line.get_xdata() == x)[0]
                 for x in points if x in line.get_xdata()]
        new   = delete(line.get_xydata(), index, axis=0)
        line.set_data(new.T)
    fig.canvas.draw()

#def find_data_points():
#    for point in points:
#        for xy in data:

# function is str
def get_ax_functions(function, ax = None):
    non_set_functions = ('axis', 'legend', 'plot', 'errorbar')
    if function not in non_set_functions:
        set_function = '_'.join(('set', function))
    else:
        set_function = function
    if ax:
        return getattr(ax, set_function)
    else:
        return globals()[function]

def set_plot_kwargs(
        x_ticks, x_label, y_label, title_='', axis_=None,
        x_scale='linear', y_scale='linear',
        labels=None, n_hist=2, n_bins=0, sub=0, journal=None,
        y_ticks=None, ax=None, x_ticklabels=None, y_ticklabels=None,
        loc='best', fancybox=True, markerscale=None, labelspacing=None,
        multi_image=False,
        panel_label='', panel_label_loc=2, panel_label_offset=.04,
        _grid=False, elev=None, azim=None, z_label=None,
        x_lim=None, y_lim=None, z_lim=None):

    # labels
    if x_label:
        get_ax_functions('xlabel', ax)(str(x_label))
    if y_label:
        get_ax_functions('ylabel', ax)(str(y_label))
    if z_label:
        get_ax_functions('zlabel', ax)(str(z_label))

    # title
    if title_:
        get_ax_functions('title',  ax)(str(title_))

    # x and y scales
    get_ax_functions('xscale', ax)(x_scale)
    get_ax_functions('yscale', ax)(y_scale)

    # legend
    if not n_bins:
        if labels:
            leg = get_ax_functions('legend', ax)(labels)
        leg = get_ax_functions('legend', ax)(
            loc=loc, markerscale=markerscale, labelspacing=labelspacing,
            fancybox=fancybox)
        if leg:
            for t in leg.get_texts():
                t.set_fontsize(fontsizes[journal]) # the legend text fontsize
                t.set_fontname(fontnames[journal])
    else:
        label_hist(
            gcf(), labels, n_hist, 0, n_bins, sub, fontsizes[journal], loc)

    # axis
    if axis_:
        get_ax_functions('axis', ax)(axis_)

    # ticks and ticklabels
    if x_ticks is not None:
        get_ax_functions('xticks', ax)(x_ticks)
    if y_ticks is not None:
        get_ax_functions('yticks', ax)(y_ticks)
    if x_ticklabels is not None:
        ax.set_xticklabels(x_ticklabels)
    if y_ticklabels is not None:
        ax.set_yticklabels(y_ticklabels)

    if _grid:
        grid(True)

    # multi-panel figures
    if panel_label:
        set_panel_labels(ax, panel_label, panel_label_loc, panel_label_offset)

    # 3d plots
    if elev and azim:
        ax.view_init(elev, azim)
    if any(x_lim):
        ax.set_xlim(x_lim)
    if any(y_lim):
        ax.set_ylim(y_lim)
    if any(z_lim):
        ax.set_zlim(z_lim)

    if not multi_image:
        draw()


# multi-panel figures
def set_panel_labels(ax, panel_label, panel_label_loc=2, panel_label_offset=.04):
    # upper right or upper left
    if panel_label_loc == 1 or panel_label_loc == 2:
        y = 1 - panel_label_offset
        if panel_label_loc == 1:
            x = 1 - panel_label_offset
        elif panel_label_loc == 2:
            x = 0 + panel_label_offset
    # upper left outside the box
    elif panel_label_loc == -1:
        x, y = -0.15, 1.05
    ax.text(x, y, panel_label, fontweight='bold', transform=ax.transAxes)

# from http://blog.olgabotvinnik.com/
# Remove top and right axes lines ("spines")
def remove_spines(ax, spines_to_remove=['top', 'right']):
    for spine in spines_to_remove:
        ax.spines[spine].set_visible(False)

# Get rid of ticks. The position of the numbers is informative enough of
# the position of the value.
def remove_ticks(ax):
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')

def plot_or_errorbar(
        x, y, error=False, yerr=None,
        color=None, marker=markers[0], linestyle=linestyles[0],
        label='', i=None, ax=None, xerr=None, filter_=None):

    #i = nonzero(y)[0]
    # cast to avoid TypeError: list indices must be integers, not NoneType
    x = array(x)
    y = array(y)

    if any(yerr):
        yerr = array(yerr)

    if x.size > y.size:
        x = delete(x, filter_)

    if error:
        if color:
            get_ax_functions('errorbar', ax)(
                x[i].flatten(), y[i].flatten(), yerr[i].flatten(), xerr,
                fmt=marker, mfc=color, ecolor=color, label=label)
        else:
            get_ax_functions('errorbar', ax)(
                x[i].flatten(), y[i].flatten(), yerr[i].flatten(), xerr,
                fmt=marker, label=label)
    else:
        if color:
            get_ax_functions('plot', ax)(
                x[i].flatten(), y[i].flatten(), c=color, ls=linestyle,
                label=label)
        else:
            get_ax_functions('plot', ax)(
                x[i].flatten(), y[i].flatten(), ls=linestyle, label=label)



###############
## HISTOGRAMS #
###############

## lable multiple histrogams: workaround (example for 2 hist i.e. hist([x0,x1]))
def label_hist(
        fig, labels, n_hist, start=0, n_bins=10, sub=0,
        fontsize=fontsizes[None], loc=0):

    axes = fig.get_axes()[sub]

    for i in range(start, n_hist):
        rec = axes.get_children()[2 + n_bins*i]
        rec.set_label(labels[i])
    leg = legend(loc=loc)

    if leg:
        # the legend text fontsize
        for t in leg.get_texts():
            t.set_fontsize(fontsize)

    fig.canvas.draw()

def rnorm_hist(fig, n_hist, start=0, n_bins=17, sub=0):
    axes = fig.get_axes()[sub]
    for i in range(start, n_hist):
        for j in range(2 + n_bins*i, 2 + n_bins*(i+1)):
            rec = axes.get_children()[j]
            rec.set_height(rec.get_height()/n_bins)
    fig.canvas.draw()

def color_hatch_hist(
        fig, facecolors, edgecolors, hatches, n_hist, start=0, n_bins=10,
        sub = 0):
    axes = fig.get_axes()[sub]
    for i in range(start, n_hist):
        for j in range(2 + n_bins*i, 2 + n_bins*(i+1)):
            rec = axes.get_children()[j]
            rec.set_facecolor(facecolors[i])
            rec.set_edgecolor(edgecolors[i])
            rec.set_hatch(hatches[i])
    fig.canvas.draw()

def hist_count_label(
        stats_, labels, colors, stat='stable', normed=False, code=None,
        zero=0):
    lixo = 10
    stats_.append(lixo)
    fig = figure()
    n, bins, patches = hist(
        stats_, bins = arange(0,1.1,0.1), histtype='bar', normed=normed)
    [setp(patch, facecolor=color)
     for patch, color in zip(patches[:-1], colors)]

    [setp(patch[0], label = labels[i]) for i, patch in enumerate(patches[:-1])]

    #y = [stat.compress((stat == 0) | (stat == 1))
     #    for stat in array(stats_[:-1])]
    y = []
    for stat in stats_[:-1]:
        stat = array(stat)
        n = stat.compress((stat == zero) | (stat == 1))
        # to use normed = True, array n as to be of same size as stats_
        #n = append(n, ones(stat.size - n.size)*10)
        y.append(n)
    y.append(lixo)
    n, bins, patches = hist(y, bins = bins, histtype='bar', normed = normed)

    if stat == 'stable':
        labels2 = ['none stable', 'all stable']
    else:
        labels2 = ['none equal', 'all equal']

    for patch, label in zip(patches[:-1], [labels2, ['',''], ['','']]):
        for i, color in zip([0,-1], ['r', 'g']):
            setp(patch[i], facecolor=color, label=label[i])

    #xlabel('% stable runs')
    #ylabel('# genomes')
    legend(loc = code)
    return fig

# data is Transposed
def bar_stacked(
        data, width=.35, journal=None, x_label='p',
        x_ticklabels=arange(0,1.1,.1), y_label=''):

    fig = figure()
    x = range(data[0].size)
    # normalize
    if data.sum(axis=0)[0] != 1:
        data = data/float(data.sum(axis=0)[0])
    [bar(
        x, data[i], width, bottom=data[:i].sum(axis=0),
        color=colors[journal][i])
     for i in range(len(data))]
    axis_ = (x[0], x[-1] + width, 0, 1)
    title_ = zip(colors[journal], range(len(data)))
    title_ = '\n'.join(
        (str(title_[:len(title_)/2]), str(title_[len(title_)/2:])))
    set_plot_kwargs(
        x, x_label, y_label, title_, axis_=axis_, journal=journal, ax=gca(),
        x_ticklabels=x_ticklabels)

    return fig



###################
## MULTIPLE PLOTS #
###################

def get_robustness_percentile(
        stability=None, robustness=None, step=.1, stable=True, norm=100.):

    p = ps_all
    r_bins = append(arange(0, 1.1, step), [1.00001])
    i = where(r_bins > 1)[0][0]

    if not any(stability):
        stability =  non_evolved_autoregulation('s')

    if not any(robustness):
        robustness = non_evolved_autoregulation('equal')

    stability  =  stability.reshape((p.size, -1))
    robustness = robustness.reshape((p.size, -1))/norm

    r_stable = array([y[where(x == stable)]
                      for x, y in zip(stability, robustness)])
    # shape(r_stable) = (p.size,): robustness of stable nets for each p
    r_i = array([digitize(x, r_bins[:i+1]) for x in r_stable])
    # shape(r_i) = (p.size,):
    # discretizes/digitizes robustness of stable nets for each p
    r_n = array([bincount(x, minlength=r_bins.size + 1) for x in r_i])
    # shape(r_n) = (p.size, r_bins.size + 1): p vs robustness
    p_n = r_n.T[1:-1]
    # shape(p_n) = (r_bins.size, p.size): robustness vs. p

    data = []
    for x in p_n:
        r = []
        for i, y in enumerate(x):
            r.extend(y*[i/10.])
        data.append(r)
    data[-2].extend(data[-1])
    return data[:-1], p_n

def get_boxplot_data(stat='s', data=None, norm=10., xboxes=[0, 1]):
    if not any(data):
        data = non_evolved_autoregulation(stat)#.reshape((ps_all.size, -1))
    boxplot_data = []
    for x in xboxes:
        n = []
        for i, s in enumerate(data):
            #n.extend([i/norm]*where(s == x)[0].size)
            n.extend([i/norm]*ma.nonzero(s == x)[0].size) #same thing
        boxplot_data.append(array(n))
    return data, array(boxplot_data)

# y = '1-p'
def get_boxplot_data2(x='s', y='p', data=None, xboxes=[0, 1]):
    if not any(data):
        # data is 2-dimensional
        data = non_evolved_autoregulation([x, y]).T.reshape((2, -1))
    return [data[1, where(data[0] == x)[0]] for x in xboxes]

def multiple_plot(
        shortnames=six_shortnames, i_ax=[0, 0, 1, 1, 2, 2], stat='p',
        colors=['k', 'r']*3, n_rows=1, n_cols=3, journal='plos',
        shortname_label=get_model_label_short_text,
        threshold=.67, n_col=2, min_n=min_N,
        axis_=(1,1e7,.2,1), x_scale='log', x_label='generations',
        fig_name='p_evolution-3panel-same_scale', titles=None,
        #marker, linestyle,
        experiment='', errors=False, save_=True):

    if experiment == 'autoregulation':
        axis_ = (-.05,1.05,0,1.2)
        x_scale = 'linear'

    # figsize = (width, height)  in inches
    fig_width = fig_widths[journal][n_col]
    fig       = figure(figsize=(fig_width*n_cols, fig_width*n_rows))
    ax        = [fig.add_subplot(n_rows, n_cols, 1)]
    ax.extend(  [fig.add_subplot(n_rows, n_cols, i, sharex=ax[0])
                 for i in range(2, n_rows*n_cols+1)])

    for shortname, i, color in zip(shortnames, i_ax, colors):

        if experiment == 'autoregulation':
            x, y, yerr = pop_autoregulation(
                shortname, min_n=min_n, plotting_=False)[0][-3]
        else:
            x, y, yerr, k = all_models_all_stats(
                [stat], [shortname], threshold=threshold, plotting_=False)[1:]

        plot_or_errorbar(
            x, y, errors, yerr, color, label=shortname_label(shortname), i=k,
            ax=ax[i])

        #ax[i].plot(x[:y.size], y, color, label = shortname_label(shortname))

        if titles is not None:
            title_ = titles[i]
        else:
            title_ = get_model_label_text(shortname)

        set_plot_kwargs(
            None, '', '', title_, axis_, x_scale, journal=journal, ax=ax[i],
            multi_image=True, panel_label=panel_labels[journal][i])

    ax[0].set_xlabel(x_label)
    ax[0].set_ylabel(get_stat_label(stat))
    show()

    if save_:
        save_fig(
            fig_name, get_figs_save_dir(journal, 'autoregulation'),
            fig_prefixes[journal], dpis[journal], formats[journal])

    return fig

# plot all_models_all_stats
def one_stat_per_panel(
        shortnames, stats_, n_rows, n_cols, fig=None, ax=[],
        journal=None, n_col=2, stat_label=get_stat_label,
        shortname_label=str, threshold=.49,
        diversity={'G.std(axis=0)':2}, vlines_=None,
        axis_=None, x_scale='log', x_label='', y_ticks=None,
        _grid=False, save_=False, fig_name=''):

    if not fig:
        fig_width = fig_widths[journal][n_col]
        fig = figure(figsize=(fig_width*n_rows, fig_width*n_cols))

    if not ax:
        ax = []

    for i, stat in enumerate(stats_):
        a = fig.add_subplot(n_rows, n_cols, i+1)

        for j, shortname in enumerate(shortnames):

            x, y, yerr, k = all_models_all_stats(
                [stat], [shortname], diversity, threshold=threshold,
                plotting_=False)[1:]

            if (shortname == 'b05p05q05' and shortnames == mutbias_shortnames):
                linestyle = linestyles[1]
            else:
                linestyle = linestyles[0]

            a.plot(
                x[k], y[k], c=colors[journal][j], ls=linestyle,
                label=shortname_label(shortname))
#        if not axis_: axis_ = axis()
        if any(vlines_):
            a.vlines(
                vlines_, axis_[2], axis_[3], color='k',
                linestyle=linestyles[-1])

        set_plot_kwargs(
            None, x_label, '', stat_label(stat), axis_, x_scale,
            journal=journal, y_ticks=y_ticks, ax=a, y_ticklabels=y_ticks,
            multi_image=True, panel_label=panel_labels[journal][i],
            _grid=_grid)

        ax.append(a)

    if save_:
        if not fig_name:
            fig_name = get_figname(stats_, shortnames)
        save_fig(fig_name)#, get_figs_save_dir(journal, experiment),
        #         fig_prefixes[journal], dpis[journal], formats[journal])

    return fig, ax

# plot all_models_all_stats
def one_shortname_per_panel(
        shortnames, stats_, n_rows, n_cols, fig=None, ax=[],
        journal=None, stat_label=get_stat_label, shortname_label=str,
        threshold=.49, diversities=[{'G.std(axis=0)':2}], vlines_=None,
        axis_=None, x_scale='log', x_label=''):

    if not fig:
        fig_width = fig_widths[journal][2]
        fig = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    if not ax:
        ax = []

    for i, shortname in enumerate(shortnames):
        a = fig.add_subplot(n_rows, n_cols, i+1)
        colors_ = iter(colors[journal])
        for stat, diversity in itertools.product(stats_, diversities):

            x, y, yerr, k = all_models_all_stats(
                [stat], [shortname], diversity, threshold=threshold,
                plotting_=False)[1:]

            a.plot(
                x[k], y[k], colors_.next(), label=stat_label(stat, diversity))

        set_plot_kwargs(
            None, x_label, '', shortname_label(shortname), axis_, x_scale,
            journal=journal, ax=a, multi_image=True,
            panel_label=panel_labels[journal][i])

        ax.append(a)

    return fig, ax

def multiple_plot2(
        shortnames, stats_, n_rows, n_cols, journal=None,
        stat_label=get_stat_label, shortname_label=str,
        threshold=.49, diversity={'G.std(axis=0)':2}, vlines_=None,
        axis_=None, x_scale='log', x_label='',
        panel_function=one_stat_per_panel,
        n_col=2, fig_name='', experiment='autoregulation', save_=True):

    if 'm4r05u1' in shortnames:
        vlines_ = (3e2, 3e3, 3e4)
        axis_   = (1,1e6,.1,1)
    if 'm4r05s1' in shortnames:
        axis_ = (1,1e6,.08,1)
    if 'm4r0b'   in shortnames:
        axis_ = (1,1e6,.1,1)

    fig_width = fig_widths[journal][n_col]
    fig = figure(figsize = (fig_width*n_cols, fig_width*n_rows))

    panel_function(
        shortnames, stats_, n_rows, n_cols, fig, [],
        journal, stat_label, shortname_label,
        threshold, diversity, vlines_, axis_, x_scale, x_label)

    show()

    if save_:
        if not fig_name:
            fig_name = get_figname(stats_, shortnames, diversity)
        save_fig(fig_name, get_figs_save_dir(journal, experiment),
                 fig_prefixes[journal], dpis[journal], formats[journal])

    return fig

#
def double_plot(
        fnames, labels, xlabels, ylabels,
        fig_width=fig_widths[None][0], fig_height=fig_height,
        color='', markers=markers, linestyles=['', ''],
        titles=['(a)', '(b)'], axis_=None,
        x_scale='linear', y_scale='linear', experiment='',
        dpi=dpis[None], format=format, bbox_inches=bbox,
        pad_inches=pad_inches, fig_prefix=fig_prefixes[None], save=True):

    fig = figure(figsize = (fig_width, fig_height))

    for i, fname in enumerate(fnames):
        x, y, z = loadtxt(paper_dir + fname + '.csv', delimiter=',')
        subplot(1,2,i+1)
        plot(x, y, color + markers[0] + linestyles[0], label=labels[i][0])
        plot(x, z, color + markers[1] + linestyles[1], label=labels[i][1])
        set_plot_kwargs(
            x, xlabels[i], ylabels[i], titles[i], axis_, x_scale, y_scale)

    if save:
        # report_dir
        fname = '.'.join((experiment, format))
        fname = '-'.join((fig_prefix, fname))
        savefig(
            report_dir + fname, dpi=dpi, format=format, bbox_inches=bbox,
            pad_inches=pad_inches)

    return fig

#
def double_hist(
        fnames, labels, xlabels, ylabels, L, bins, dx, y_max, normed=True,
        align='left', fig_width=fig_widths[None][0], fig_height=fig_height,
        color='', markers=markers, linestyles=['', ''],
        titles=['(a)', '(b)'], axis_=None,
        x_scale='linear', y_scale='linear', experiment='',
        dpi=dpis[None], format=format, bbox_inches=bbox,
        pad_inches=pad_inches, fig_prefix=fig_prefixes[None],
        save_dir=get_figs_save_dir(None), save=True):

    fig = figure(figsize = (fig_width, fig_height))
    n_hist = 2

    for i, fname in enumerate(fnames):
        data = loadtxt(paper_dir + fname + '.csv', delimiter=',')
        subplot(1,2,i+1)
        y = data[0]
        z = insert(data[1], 0, -1) if data[0].size == data[1].size else data[1]
        if any(L):
            y /= L[i]
            z /= L[i]
        n, bins[i], patches = hist(
            [y, z], bins[i], normed = normed, align = align)
        n_bins = bins[i].size - 1
        color_hatch_hist(
            fig, ['k', 'w'], ['k', 'k'], ['', 'x'], n_hist, 0, n_bins, i)
        rnorm_hist(fig, n_hist, 0, n_bins, i)
        if any(L):
            axis_ = [-dx[i]/L[i], (L[i] + dx[i])/L[i], 0, y_max[i]]

        set_plot_kwargs(
            None, xlabels[i], ylabels[i], titles[i], axis_, x_scale, y_scale,
            labels[i], n_hist, n_bins, i)

    if save:
        save_fig(
            experiment, save_dir, fig_prefix, dpi, format, bbox, pad_inches)

    return fig



########################
## PHENOTYPE FUNCTIONS #
########################

def get_gene_count_variation(gene, min = min_state):
    gene = list(gene)
    return array([gene.count(min), gene.count(1)]).min() / float(len(gene))

def get_all_phenotype_variation(
        phenotypes, min = min_state, var = get_gene_count_variation):

    return array([[var(gene) for gene in cycle.T] for cycle in phenotypes])

def get_all_vector_phenotypes(
        phenotypes, bits, basename = def_basename, base = def_base):

    return (get_vector_phenotypes(phenotype, bits, basename, base)
            for phenotype in phenotypes)

def get_vector_phenotypes(
        phenotypes, bits, basename = def_basename, base = def_base):

    return array([get_vector(decimal, basename, base, bits)
                  for decimal in phenotypes])

def get_fraction_of_fixed_points_count(var):
    return 1 - nonzero(var.mean(axis=0))[0].size / float(shape(var)[-1])

def get_fraction_of_fixed_points_mean(var):
    return 1 - var.mean()

def get_cycling_genes_vs_p(individuals, phenotypes, myrun,
                           var = get_gene_count_variation,
                           fp = get_fraction_of_fixed_points_count):
    # get vector phenotypes
    phenotypes = get_all_vector_phenotypes(
        phenotypes, myrun.bits, myrun.basename, myrun.base)
    all_var =  get_all_phenotype_variation(
        phenotypes, myrun.min, var).reshape(
            (myrun.samples, 2**myrun.bits, myrun.bits))
    all_var /= all_var.max()

    # load functions
    a = get_activating_fraction
    d = diag

    return array([(a(d(genotype)), fp(var))
                  for genotype, var in zip(individuals, all_var)])

def get_fixed_points_from_tuple_phenotypes(phenotypes):
    return array([phenotype for phenotype in phenotypes if len(phenotype) == 1])

def get_stability_from_phenotypes(filename):
    phenotypes = cPickle.load(open(filename))
    fixed = get_fixed_points_from_tuple_phenotypes(phenotypes)
    return fixed.size

def plot_phenotype_distribution(bit, filename, all_bits = False):
    phenotypes = array(cPickle.load(open(data_dir + filename)))
    fixed = get_fixed_points_from_tuple_phenotypes(phenotypes).flatten()
    fig = figure()
    if all_bits:
        hist(fixed, range(2**bit))
    else:
        hist(fixed)
    fixed = get_vector_phenotypes(fixed, bit).flatten()
    fig = figure()
    hist(fixed)



###################
## FILE FUNCTIONS #
###################

def load_stability0(
        bit, density, min, dim, noise, suffix, prefix=data_dir + 'stability',
        extension='.dat'):
    filename = get_filename(bit, density, min, dim, noise, precision)
    return cPickle.load(open(prefix + filename + suffix + extension))

def load_stability(
        bit, samples, density, min, dim, noise, suffix, min_n=min_N,
        min_samples = min_Samples, prefix=data_dir + 'stability',
        extension='.dat', precision=2, verbose=False):
    filename = get_filename(bit, density, min, dim, noise, precision)
    try:
        # .dat file
        file_dat = prefix + filename + suffix + extension
        data = cPickle.load(open(file_dat))
    except IOError:
        if verbose:
            print 'cannot open', file_dat
        try:
            # .tsv file
            file_tsv = prefix + filename + suffix + '.tsv'
            samples, n = loadtxt(file_tsv, int)[-1]
        except IOError:
            if verbose:
                print 'cannot open', file_tsv
            try:
                # phenotypes
                file_phe = (data_dir + 'phenotypes' + filename + suffix +
                            extension)
                n = get_stability_from_phenotypes(file_phe)
            except IOError:
                if verbose:
                    print 'cannot open', file_phe
                return 0
            else:
                if verbose:
                    print 'opening', file_phe
                if n > min_n:
                    return n / samples
        else:
            if verbose:
                print 'opening', file_tsv
            if n > min_n and samples > min_samples:
                return n / float(samples)
    else:
        if verbose:
            print 'opening', file_dat
        n =  data['%d' %min]['%d' %dim]['%d' %bit]['%g' %density]
        if n > min_n:
            return n / samples
    return 0

def load_stability_cSigmoid(
        steepness, converge_time, converge_threshold, devo_time, filter,
        filename, prefix=data_dir +' stability', extension='.dat'):
    try:
        suffix  = '_a%sT%dE%.0e' %(steepness, converge_time, converge_threshold)
        if devo_time:
            suffix += 'M%s' %devo_time
        suffix += '-filter%g' %filter
        file = open(prefix + filename + suffix + extension)
    except: #IOError ?
        suffix = ('_a%.2fT%dE%.0e' %
                  (steepness, converge_time, converge_threshold))
        if devo_time:
            suffix += 'M%s' %devo_time
        suffix += '-filter%g' %filter
        file = open(prefix + filename + suffix + extension)
    return cPickle.load(file)

def load_viability(
        suffix, bit, density, min = -1., dim = 2, noise = 0, min_n = min_N,
        min_samples = min_Samples, prefix = data_dir + 'viability',
        extension = '.tsv', precision = 2, verbose = False):

    # filename_exceptions
    precision = density_precisions[bit] if bit in density_precisions else 2
    filename = get_filename(bit, density, min, dim, noise, precision)
    if verbose:
        print filename + suffix
    try:
        samples, n = loadtxt(prefix + filename + suffix + extension, int)[-1]
    except: #IOError?
        if verbose:
            print 'Cant open file'
        return 0
    else:
        if n > min_n and samples > min_samples:
            return n / float(samples)
        else:
            return 0

def load_samples(
        bit, samples, density, min, dim, noise,
        suffix, prefix = data_dir + 'stability', extension = '.dat',
        precision = 2, verbose = False):

    n_samples = 0
    filename  = get_filename(bit, density, min, dim, noise, precision)
    if verbose:
        print filename
    try:
        # .dat file
        data = cPickle.load(open(prefix + filename + suffix + extension))
    except: #IOError ?
        # .tsv file
        try:
            samples, n = loadtxt(prefix + filename + suffix + '.tsv', int)[-1]
        except: #IOError ?
            pass
        else:
            n_samples = max(n_samples, samples)
    else:
        n_samples = max(n_samples, samples)
    return n_samples

# old
def save_plot(
        filename, experiment='stability', save_dir=figs_dir, dpi=None,
        format='pdf', separator='-'):
    fname = '.'.join((filename, format))
    if experiment:
        fname = separator.join((experiment, fname))
    savefig(save_dir + fname, dpi=dpi, format=format)

# new; report_dir
def save_fig(
        fig_name, save_dir=get_figs_save_dir(None),
        fig_prefix=fig_prefixes[None],
        dpi=dpis[None], format=formats[None],
        bbox_inches=bbox, pad_inches=pad_inches, separator='-',
        fig_number=None):

    fname = '.'.join((fig_name, format))

    if fig_number:
        fname = separator.join(('Fig_%s'%fig_number, fname))

    if fig_prefix:
        fname = separator.join((fig_prefix, fname))

    savefig(
        save_dir + fname, dpi=dpi, format=format, bbox_inches=bbox,
        pad_inches=pad_inches)



################
## REGRESSIONS #
################

def get_regression_label(x_fit, a, b, x='N', y='y', label=None):
    if label:
        return label
    if x_fit == array:
        return '$%s = %.3g e^{%.2g%s}$' %(y, exp(b), a, x)
    if x_fit == log:
        return '$%s = %.2g %s^{%.2g}$'  %(y, exp(b), x, a)

# linear regression in lin-log (exponential fit) or
# log-log (power-law: '$y \sim N^{a}$') scales
def polyfit_regression(
        x, y, w=None, x_min=None, x_max=None,
        x_fits1=[array], x_fits2=[], y_fit=log, x_step=1,
        x_label='N', y_label='y', label=None,
        color='b', linestyle=linestyles[1], plotting_=True):

    # default behaviors
    x = array(x)
    y = array(y)
    if not x_min:
        x_min = x[0]
    if not x_max:
        x_max = x[-1]
    if not x_fits2:
        x_fits2 = x_fits1

    # regression range
    xx = arange(x_min, x_max+1*x_step, x_step)

    # removes x or y zeros because of log()
    #i = slice(None)
    if log in x_fits1 or (log in x_fits2 and any(w)):
        i = nonzero(x)[0]
        x = x[i]
        y = y[i]
    if y_fit == log:
        i = nonzero(y)[0]
        x = x[i]
        y = y[i]

    # the regression
    for x_fit1, x_fit2 in zip(x_fits1, x_fits2):

        # mean
        ya, yb = polyfit(x_fit1(x), y_fit(y), 1) # with data
        yy = ya*x_fit1(xx) + yb                  # full
        if y_fit == log:
            yy = exp(yy)

        if plotting_:
            label = get_regression_label(
                x_fit1, ya, yb, x_label, y_label, label)
            plot(xx, yy, color + linestyle, label=label)

        # max
        if any(w):
            wa, wb = polyfit(x_fit2(x), y_fit(w), 1) # with data
            wy = wa*x_fit2(xx) + wb                  # full
            if y_fit == log:
                wy = exp(wy)

            if plotting_:
                label = get_regression_label(x_fit2, wa, wb)
                plot(xx, wy, 'k' + linestyle, label=label)

            # max: powers of two
            powers_wy = 2**ceil(log2(wy)).astype(long) # full

            if plotting_:
                label = get_regression_label(x_fit2, wa, wb) + ' $2^n$'
                plot(xx, powers_wy, 'g' + linestyle, label=label)

    if any(w):
        return xx, yy, wy, powers_wy
    else:
        return xx, yy, ya, yb



##############
## STABILITY #
##############

## STABILITY VS. C (OLD)

def stability_vs_c_vs_k(
        dims, min, bit = 4, noise=0, samples=1e6,
        full_enum=False, symm=False, stable=False, binary=False,
        filter=None, experiment='stability', extension ='.dat',
        stability={}):

    densities = arange(1, bit+1.) / bit
    G = array([float(sizes[str(bit)][density]['G'])
               for density in map(str, densities)])
    G_red = array([float(sizes[str(bit)][density]['G_red'])
                   for density in map(str, densities)])

    prefix = data_dir + experiment
    suffix = get_suffix(samples, full_enum, symm, stable, binary, filter)

    if '%d' %min not in stability.keys():
        stability['%d' %min] = {}

    fig = figure()

    for dim in dims:
        y = [load_stability0(
            bit, density, min, dim, noise, suffix, prefix, extension)
             ['%d' %min]['%d' %dim]['%d' %bit]['%g' %density]
             for density in densities]

        if full_enum:
            L = float(sizes[str(bit)][str(dim)]['L'])
            if symm:
                normalize = L * G_red
            else:
                normalize = L * G
        else:
            normalize = samples

        y  = array(y, dtype = float)
        y /= normalize
        plot(densities, y, 'o--', label = 'K = %d' %dim)

        stability['%d' %min]['%d' %dim] = y

    set_plot_kwargs(
        densities, 'c', experiment, 'N = %s, min = %d, ' %(bit, min) + suffix)

    prefix = figs_dir + experiment
    filename = '_vs_c_vs_k-N%s_min_%s-noise%s-' %(bit, [min], noise)
    extension = format
    savefig(prefix + filename + suffix + extension)

    return stability

def stability_vs_c_vs_N(
        bits, min, dim=2, noise=0, samples=1e6,
        stable=False, binary=False, filter=None,
        experiment= 'stability', extension='.dat',
        stability={}, axis_=None, save=True):

    normalize = samples
    basename = get_base(dim, min=min)[1]

    prefix = data_dir + experiment
    suffix = get_suffix(
        samples, stable=stable, binary=binary, filter=filter)

    if '%d' %min not in stability.keys():
        stability['%d' %min] = {}

    fig = figure()

    for bit in bits:
        #densities = arange(c, bit+1.) / bit
        n_zeros   = arange(bit, dtype = float)
        densities = (c-n_zeros/bit)[::-1]
        y = [load_stability0(
            bit, density, min, dim, noise, suffix, prefix, extension)
             ['%d' %min]['%d' %dim]['%d' %bit]['%g' %density]
             for density in densities]

        y  = array(y, dtype = float)
        y /= normalize
        plot(densities, y, 'o--', label = 'N = %d' %bit)

        stability['%d' %min]['%d' %bit] = y

    set_plot_kwargs(
        densities, 'c', experiment, suffix + ', %s' %basename, axis_)

    if save:
        filename = '_vs_c_vs_N_%s_%s-noise%s-' %(bits, basename, noise)
        save_plot(filename + suffix, experiment)

    return stability


## STABILITY VS. N (OLD)

def stability_vs_N(
        bits, min, dim=2, degrees=[4], noise=0, samples=1e6, stable=False,
        binary=False, filter=None, experiment='stability', extension ='.dat',
        stability={}, axis_=None, save=True, linestyle='o--'):

    degrees = array(degrees, dtype=float)
    normalize = samples
    basename = get_base(dim, min=min)[1]
    prefix = data_dir + experiment
    suffix = get_suffix(samples, stable=stable, binary=binary, filter=filter)

    fig = figure()

    for i, degree in enumerate(degrees):
        degree += c-1
        y = [load_stability0(bit, degree/bit, min, dim, noise, suffix, prefix, extension)['%d' %min]['%d' %dim]['%d' %bit]['%g' %(degree/bit)] for bit in bits[i:]]

        y  = array(y, dtype = float)
        y /= normalize
        plot(bits[i:], y, linestyle, label = 'cN = %d' %degree)

    set_plot_kwargs(bits, 'N', experiment, suffix + ', %s' %basename, axis_)

    if save:
        filename = '_vs_N_vs_degree_%s-noise%s-' %(basename, noise)
        save_plot(filename + suffix, experiment)

    return y


## STABILITY with cSIGMOID2 (OLD)

def stability_vs_a(
        a, converge_time=T, converge_threshold=converge, devo_time=None,
        samples=1e6, filter=.01, linestyle='o--'):

    dim = 2
    min = -1.
    bit = N
    density = c
    noise = 1e-100
    experiment = 'stability'
    extension  = '.dat'

    prefix   = data_dir + experiment
    filename = ('-N_%s_c_%s_min_%s_dim_%s-noise%s-random-n%.2g-cSigmoid2' %
                ([bit], get_basename(array([density])), [min], [dim], noise,
                 samples))

    y = [load_stability_cSigmoid(
        steepness, converge_time, converge_threshold, devo_time, filter,
        filename, prefix, extension)['%d' %min]['%d' %dim]['%d' %bit][
            '%g' %density]
         for steepness in a]

    y  = array(y, dtype = float)
    y /= samples
    fig = figure()
    plot(a, y, linestyle)
    xlabel('a')
    ylabel(experiment)
    suffix = '_T%dE%.0e-filter%g' %(converge_time, converge_threshold, filter)
    title(filename + suffix)

    return y

def stability_vs_T(
        T, steepness=a, converge_threshold=converge, devo_time=None,
        samples=1e6, filter=.01, linestyle='o--'):

    dim = 2
    min = -1.
    bit = N
    density = c
    noise = 1e-100
    experiment = 'stability'
    extension  = '.dat'
    basename  = get_base(dim, min=min)[1]

    prefix   = data_dir + experiment
    filename = ('-N_%s_c_%s_min_%s_dim_%s-noise%s-random-n%.2g-cSigmoid2' %
                ([bit], get_basename(array([density])), [min], [dim], noise, samples))

    y = [load_stability_cSigmoid(steepness, converge_time, converge_threshold, devo_time, filter, filename, prefix, extension)['%d' %min]['%d' %dim]['%d' %bit]['%g' %density] for converge_time in T]

    y  = array(y, dtype = float)
    y /= samples
    fig = figure()
    plot(T, y, linestyle)
    xlabel('T')
    ylabel(experiment)
    suffix = '_a%.2fE%.0e-filter%g' %(steepness, converge_threshold, filter)
    title(filename + suffix)

    return y

def stability_vs_E(converge, steepness = a, converge_time = T, devo_time = None, samples = 1e6, filter = .01, linestyle = 'o--'):

    dim = 2
    min = -1.
    bit = N
    density = c
    noise = 1e-100
    experiment = 'stability'
    extension  = '.dat'
    basename  = get_base(dim, min=min)[1]

    prefix   = data_dir + experiment
    filename = ('-N_%s_c_%s_min_%s_dim_%s-noise%s-random-n%.2g-cSigmoid2' %
                ([bit], get_basename(array([density])), [min], [dim], noise, samples))

    y = [load_stability_cSigmoid(steepness, converge_time, converge_threshold, devo_time, filter, filename, prefix, extension)['%d' %min]['%d' %dim]['%d' %bit]['%g' %density] for converge_threshold in converge]

    y  = array(y, dtype = float)
    y /= samples
    fig = figure()
    plot(converge, y, linestyle)
    xscale('log')
    xlabel('$\epsilon$')
    ylabel(experiment)
    suffix = '_a%.2fT%d-filter%g' %(steepness, converge_time, filter)
    title(filename + suffix)

    return y


## STABILITY VS. N (NEW)

def stability_vs_N2(
        bits, density=c, degree=None, min=-1., dim=2, noise=0, min_n=min_N,
        min_samples=min_Samples, sample_sizes=[1e6], devo_times=[inf],
        alt_filenames=['path_leng_vs', 'power_vs_mu'], alt_bits=N_dict,
        binary=False, graph=False, scale_free=False, fix_gpmap=False,
        gpmap=None, n_zeros=None, samples2=1e6,
        experiment='stability', extension ='.dat', verbose=False):

    if experiment == 'stability_vs_N':
        experiment = 'stability'
    prefix = data_dir + experiment
    extra = 0
    stability = []
    for i, bit in enumerate(bits):

        if degree:
            density  = None
            samples2 = 1e8
            alt_bits = 25
        #if verbose: print density, n_zeros, degree, bit

        # filename_exceptions
        if bit in density_precisions:
            precision = density_precisions[bit]
        else:
            precision = 2

        # extra measurements available: paper_dir .csv
        if (density == c and inf in devo_times and min == -1. and not binary and
            not fix_gpmap and not graph):
            if bit < 11:
                file_csv = 'figs/paper/figure_1-stability_vs_N.csv'
                if verbose:
                    print 'opening', file_csv
                stability = append(
                    stability, loadtxt(file_csv, delimiter=',')[-1][i])
            else:
                stability = append(stability, 0)
            extra += 1

        # extra measurements available: path_leng_vs_N and power_vs_mu_N
        if inf in devo_times and min == -1.:
            for filename in alt_filenames:
                if bit > alt_bits:
                    filename += '_N%d' %bit
                    if degree:
                        filename += 'k%d' %degree
                    try:
                        file_alt = data_dir + filename + '.txt'
                        data = loadtxt(file_alt)
                    except: #IOError ?
                        if verbose:
                            print 'cannot open', file_alt
                        n = 0
                    else:
                        if verbose:
                            print 'opening', file_alt
                        n = len(data)/samples2
                else:   n = 0
                stability = append(stability, n)
            extra += len(alt_filenames)

        # extra measurements available: noise 1e-100
        if density == c and not binary:
            samples   = 1e7
            stability = append(
                stability,
                load_stability(
                    bit, samples, density, min, dim, 1e-100,
                    get_suffix(samples), min_n, min_samples, prefix, extension,
                    precision, verbose))
            extra += 1

        # main files for N > 10: stability-N
        density = get_density_degree_nzeros(density, n_zeros, degree, bit)[0]
        #if verbose: print density, n_zeros, degree, bit
        for samples in sample_sizes:
            for devo_time in devo_times:
                suffix = get_suffix(
                    samples, binary=binary, fix_gpmap=fix_gpmap,
                    gpmap=gpmap, devo_time=devo_time,
                    graph=graph, scale_free=scale_free)
                #if verbose: print suffix
                stability = append(
                    stability,
                    load_stability(
                        bit, samples, density, min, dim, noise, suffix, min_n,
                        min_samples, prefix, extension, precision, verbose))

    extra /= len(bits)
    shape = [len(bits), len(sample_sizes)*len(devo_times) + extra]

    stability = stability.reshape(shape).T.max(0)
    index = nonzero(stability)
    x = array(bits)[index]
    y = stability[index]

    return x, y

#
def plot_stability_vs_N(
        bits=arange(4,101), density=c, degree=None, min=-1., dim=2,
        min_n=min_N, min_samples=min_Samples,
        sample_sizes=[1e8, 1e7, 1e6], devo_times=[inf, 'maximum', 'mean'],
        alt_bits=N_dict, noise=0, binary=False, graph=False,
        scale_free=False, fix_gpmap=False, gpmap=None,
        alt_filenames=['path_leng_vs', 'power_vs_mu'], regression=None,
        color='', linestyles = markers, x_scale='log', y_scale='log',
        experiment='stability_vs_N', extension='.csv', plotting_=True,
        save_=True, verbose=False):

    max_bit = bits[-1]

    if plotting_:
        fig = figure()
    if plotting_ or save_:
        for devo_time, linestyle in zip(devo_times, linestyles):
            x, y = stability_vs_N2(
                bits, density, degree, min, dim, noise, min_n, min_samples,
                sample_sizes, [devo_time], alt_filenames, alt_bits, binary,
                graph, scale_free, fix_gpmap, gpmap, experiment=experiment,
                verbose=verbose)

            if plotting_:
                plot (x, y, color + linestyle, label=str(devo_time))

            if save_:
                prefix   = data_dir + experiment + print_min_max(x, x)
                filename = '-c%g' %density if density else '-k%d' %degree
                suffix   = '-devo_time_%s' %str(devo_time)
                if graph:
                    suffix += '-graph'
                    if scale_free:
                        suffix += '-scale_free'
                    else:
                        suffix += '-exp_pow'
                if fix_gpmap:
                    suffix += '-' + gpmap.__name__
                savetxt(
                    prefix + filename + suffix + extension, (x, y),
                    delimiter=',')

        max_bit = xticks()[0][-1]

    # regression == 0 is True
    if regression is not None:
        # max of all data: devo_times = [inf, 'maximum', 'mean']
        x, y = stability_vs_N2(
            bits, density, degree, min, dim, noise, min_n, min_samples,
            sample_sizes, devo_times, alt_filenames, alt_bits, binary, graph,
            scale_free, fix_gpmap, gpmap, experiment=experiment,
            verbose=verbose)

        # fit the data
        if not density:
            x_fit = log
            label = '$y \sim N^{-a}$'
        else:
            x_fit = array
            label = '$y \sim exp^{-aN}$'
        ya, yb = polyfit(x_fit(x[regression:]), log(y[regression:]), 1)
        # range of fitted values
        max_bit = max(max_bit, x[-1])
        xx = arange(bits[0], max_bit + 1)
        yy = exp(ya*x_fit(xx) + yb) # full
        if plotting_:
            plot(xx, yy, 'k--', label = label)

    if plotting_:
        title_ = '%s' %gpmap.__name__ if fix_gpmap else None
        set_plot_kwargs(None, 'N', 'stability', title_, None, x_scale, y_scale)
        if save_:
            save_plot(
                print_min_max(bits, [max_bit]) + filename + suffix, experiment)

    #returns max of all devo_times if regression is not None and
    # plotting_ = save_ = False (see plot_stability_vs_N_vs_k)
    return x, y, fig

#
def plot_stability_vs_N_vs_k(
        bits=None, densities=[None, None, c], degrees=[2, 4, None],
        sample_sizes=[[1e8], [1e8], [1e8, 1e7, 1e6]],
        devo_times=[inf, 'maximum', 'mean'], alt_bits=N_dict, min=-1.,
        dim=2, min_n = min_N, min_samples=min_Samples, regression=None,
        vline=None, color='', linestyle=linestyles[0], markers=markers,
        axis_=[3, 1.1e4, 1e-6, 1], x_scale='log', y_scale='log',
        experiment='stability_vs_N', extension = '.csv',
        format=format, bbox_inches=bbox, pad_inches=pad_inches,
        fig_prefix=fig_prefixes[None], save_dir=get_figs_save_dir(None),
        filter_=None, labels=False, save_=True, verbose=False):

    fig = figure()

    if not any(bits):
        bits = stability_bits[2]
    for density, degree, sample_size, marker in zip(
            densities, degrees, sample_sizes, markers):
        key = degree if degree else density
        if key != 2 and key !=1:
            bit = bits[:where(bits == max_bits[experiment][key])[0] + 1]
        else:
            bit = bits

        # regression not None garantees the max (x, y) of devo_times
        x, y = plot_stability_vs_N(
            bit, density, degree, min, dim, min_n, min_samples, sample_size,
            devo_times, alt_bits, regression=0, plotting_=False,
            save=False, verbose=verbose)

        if save_:
            prefix   = data_dir + experiment + print_min_max(x, x)
            filename = '-c%g' %density if density else '-k%d' %degree
            savetxt(prefix + filename + extension, (x, y), delimiter=',')

        if not degree:
            degree = 'N'
        plot(x, y, color + marker, label = 'K = %s' %degree)

        if regression is not None:
            # fit the data
            if not density:
                x_fit = log
                label = '$y \sim N^{-a}$' if labels else ''
            else:
                x_fit = array
                label = '$y \sim exp^{-aN}$' if labels else ''
            ya, yb = polyfit(x_fit(x[regression:]), log(y[regression:]), 1)
            yy = exp(ya*x_fit(x) + yb) # full
            plot(x, yy, 'k--', label = label)

    #if delete_x: delete_data_points(fig, delete_x)
    if filter_:
        delete_data_points(fig, filter_)

    if vline:
        axvline(
            vline, color = 'k', linestyle = linestyle,
            label = 'e.coli\nN = %d' %vline)

    set_plot_kwargs(None, 'N', 'stability', '', axis_, x_scale, y_scale)

    if save_:
        # figs_dir
        #save_plot(print_min_max(bits, bits) + '_vs_k', experiment)
        # report_dir
        #fname = '.'.join((experiment, format))
        #fname = '-'.join((fig_prefix, fname))
        #savefig(
         #   report_dir + fname, dpi = dpi, format = format,
          #  bbox_inches = bbox, pad_inches = pad_inches)
        save_fig(
            experiment, save_dir, fig_prefix, format=format,
            bbox_inches=bbox_inches, pad_inches=pad_inches)

    return fig


## STABILITY VS. C (NEW)

def stability_vs_c_vs_N2(
        bits, min=-1., dim=2, noise=0, samples=1e8,
        devo_time=inf, axis_=None, save_=True, verbose=False):

    suffix = get_suffix(samples, devo_time = devo_time)
    fig = figure()

    for bit in bits:

        densities = array([get_density_degree_nzeros(None, None, degree, bit)[0]
                           for degree in range(1, bit+1)])

        stability = [
            load_stability(
                bit, samples, density, min, dim, noise, suffix, verbose=verbose)
            for density in densities]

        plot(densities, stability, 'o--', label = 'N = %d' %bit)

        if save_:
            savetxt(
                paper_dir + 'figure_6a-stability_vs_c_N%d.csv' %bit,
                (densities, stability), delimiter=',')

    set_plot_kwargs(None, 'c', 'stability', '', axis_ = axis_)

    #if save_: save_plot( '_vs_c-N%s' %(bits))

    return stability

def plot_stability_diff_vs_N(
        alternative=False, fname='add_file3-stability_diff',
        format=format, save_=True):
    x1, y1 = loadtxt(data_dir + 'stability_vs_N4_60-c1.csv',   delimiter = ',')
    x2, y2 = loadtxt(data_dir + 'stability_vs_N4_3000-k2.csv', delimiter = ',')
    x4, y4 = loadtxt(data_dir + 'stability_vs_N4_175-k4.csv',  delimiter = ',')
    x6, y6 = loadtxt(data_dir + 'stability_vs_N6_60-k6.csv',   delimiter = ',')
    x12    = [x for x in x1 if x in x2]
    x14    = [x for x in x1 if x in x4]
    x16    = [x for x in x1 if x in x6]
    yy12   = array([y for i, y in enumerate(y1) if x1[i] in x12])
    yy21   = array([y for i, y in enumerate(y2) if x2[i] in x12])
    yy14   = array([y for i, y in enumerate(y1) if x1[i] in x14])
    yy41   = array([y for i, y in enumerate(y4) if x4[i] in x14])
    yy16   = array([y for i, y in enumerate(y1) if x1[i] in x16])
    yy61   = array([y for i, y in enumerate(y6) if x6[i] in x16])
    z21    = yy21 - yy12
    z41    = yy41 - yy14
    z61    = yy61 - yy16
    fig = figure()
    plot(x12, z21, 'ko', label = "K' = 2")
    plot(x14, z41, 'ks', label = "K' = 4")
    plot(x16, z61, 'kv', label = "K' = 6")
    set_plot_kwargs(None, 'N', "stability(K=K') - stability(K=N)", '')
    axis([2, 61, -0.005, 0.07])
    if alternative:
        x42  = [x for x in x4 if x in x2]
        yy42 = array([y for i, y in enumerate(y4) if x4[i] in x42])
        yy24 = array([y for i, y in enumerate(y2) if x2[i] in x42])
        z24  = yy24 - yy42
        plot(x42, z24, 'b*', label = "K'=2 - K=4")
        x62  = [x for x in x6 if x in x2]
        yy62 = array([y for i, y in enumerate(y6) if x6[i] in x62])
        yy26 = array([y for i, y in enumerate(y2) if x2[i] in x62])
        z26  = yy26 - yy62
        plot(x62, z26, 'g*', label = "K'=2 - K=6")
        save_ = False

    fname = '.'.join((fname, format))
    if save_:
        savefig(report_dir + fname, dpi=dpi, format=format,
                bbox_inches=bbox, pad_inches=pad_inches)
    return x12[where(z21 == z21.max())[0]], x14[where(z41 == z41.max())[0]], x16[where(z61 == z61.max())[0]]


## STABILITY VS. P

# bins and min_samples are for compatibility with get_bins_mean_std()
def get_unique_mean_std(
        data_x, data_y, bins=None, min_samples=min_Samples,
        plotting_=False, x_label='', y_label='', axis_=None):

    x    = unique(data_x)
    y    = [data_y[where(data_x == t)[0]].mean() for t in x]
    yerr = [data_y[where(data_x == t)[0]].std()  for t in x]

    if plotting_:
        fig = figure()
        errorbar(x, y, yerr, fmt = 'o')
        set_plot_kwargs(None, x_label, y_label, '', axis_)

    return x, y, yerr

# NOTICE: data_y and _x are swapped!
def get_bins_mean_std(data_y, data_x, bins, min_samples=min_Samples):
    x    = [bins[0]]
    y    = []
    yerr = []
    for bin in bins[1:]:
        i = where(data_x < bin)[0]
        if i.size > min_samples:
            y.append(data_y[i].mean())
            yerr.append(data_y[i].std())
        else:
            x.pop()
        x.append(bin)
        data_x = delete(data_x, i)
        data_y = delete(data_y, i)
    if i.size > min_samples:
        y.append(data_y.mean())
        yerr.append(data_y.std())
    else:
        x.pop()
    return x, y, yerr

def plot_stability_vs_p(
        filenames, get_x_y_yerr=get_unique_mean_std, bins=None,
        min_samples=min_Samples, bits=(4, 10),
        x_label='p', y_label='stability', title_='',
        axis_=(-.05, 1.05, -.2, 1.2), x_ticks=True, x_scale='linear',
        color='', linestyle=linestyles[0], markers=markers,
        vline=None, experiment='stability_vs_p', extension='.txt',
        save_dir=get_figs_save_dir(None), error=True, save_=False):

    fig = figure()

    for filename, bit, marker in zip(filenames, bits, markers):
        ind_trace, p, stability = loadtxt(data_dir + filename + extension).T
        x, y, yerr = get_x_y_yerr(
            p, stability/stability.max(), bins, min_samples)

        if error:
            errorbar(
                x, y, yerr, fmt=marker, mfc=color, ecolor=color,
                label='N = %d' %bit)
        else:
            plot(x, y, color + marker, label='N = %d' %bit)

    if any(bins):
        x_label, y_label = y_label, x_label
        axis_ = (1e-3, 1, .2, 1)
        experiment = 'p_vs_stability'

    if vline:
        axvline(
            vline, color='k', linestyle=linestyle,
            label='e.coli: p = %.2f' %vline)

    if any(bins):
        x_ticks = None
    if x_ticks:
        x_ticks = x
    set_plot_kwargs(x_ticks, x_label, y_label, title_, axis_, x_scale)

    if save_:
        save_fig(experiment, save_dir)

    return x, y, yerr, fig


## STABILITY VS. GPMAP

'''def get_stability_from_phenotypes(gpmap, bit, density = None, degree = 2, noise = 0, devo_time = inf, samples = 1e4, min = -1.):
    density, n_zeros, degree = get_density_degree_nzeros(density, None, degree, bit)
    experiment = 'phenotypes'
    prefix   = data_dir + experiment
    filename = get_filename(bit, density, min, 2, noise)
    suffix   = get_suffix(samples, binary = True, fix_gpmap = True, gpmap = gpmap, devo_time = devo_time)
    phenotypes = cPickle.load(open(prefix + filename + suffix + '.dat'))
    return nonzero(phenotypes)[0].size / float(len(phenotypes))
'''

##################
## STABLE STATES #
##################

def stable_states_vs_c_vs_k(
        dims, min, bit=4, noise=0, samples=1e6, full_enum=True, symm=False,
        stable=False, binary=False, normalize=False, experiment='counting',
        extension='.dat', counting={}, y_scale='linear', save_=True,
        linestyle='o--', potential=False):

    densities = arange(1, bit+1.) / bit
    G = array([float(sizes[str(bit)][density]['G'])
               for density in map(str, densities)])
    G_red = array([float(sizes[str(bit)][density]['G_red'])
                   for density in map(str, densities)])

    prefix = data_dir + experiment
    if full_enum:
        suffix = 'binary-enum'
        if symm:
            suffix += '-symm'
        else:
            suffix += '-full'
    elif stable:
        suffix = 'stable-n%.2g' %(samples)
    elif binary:
        suffix = 'binary-n%.2g' %(samples)
    else:
        suffix = 'random-n%.2g' %(samples)

    if '%d' %min not in counting.keys():
        counting['%d' %min] = {}

    fig = figure()

    for dim in dims:

        # state space size, ie, total number of possible phenotypes
        L = 'L_red' if symm else 'L'
        L = float(sizes[str(bit)][str(dim)][L])

        y = []
        for density in densities:
            filename = get_filename(bit, density, min, dim, noise)
            try:
                data = cPickle.load(
                    open(prefix + filename + suffix + extension))
            except: #IOError ?
                n = loadtxt(prefix + filename + suffix + '.tsv', int)[-1][-1]
            else:
                n = data['%d' %min]['%d' %dim]['%d' %bit]
            y.append(n)

        if normalize:
            y  = array(y, dtype=float)
            # state space size, ie, total number of possible phenotypes
            L = 'L_red' if symm else 'L'
            L = float(sizes[str(bit)][str(dim)][L])
            if normalize == 'L':
                y /= L
            if normalize == 'G':
                if symm:
                    y /= G_red
                else:
                    y /= G
            if normalize == 'LG':
                if symm:
                    y /= L * G_red
                else:
                    y /= L * G

        plot(densities, y, linestyle, label='K = %d' %dim)
        if potential:
            plot(densities, [L]*densities.size, 'k--')

        counting['%d' %min]['%d' %dim] = y

    xticks(densities)
    xlabel('c')
    experiment = 'stable_states'
    ylabel(experiment)
    if normalize:
        suffix += '-%s' %normalize
    title('N = %s, min = %d, ' %(bit, min) + suffix)
    leg = legend(loc=0)
    yscale(y_scale)

    if save_:
        prefix = figs_dir + experiment
        filename = '_vs_c_vs_k-N%s_min_%s-noise%s-' %(bit, [min], noise)
        extension = format
        savefig(prefix + filename + suffix + extension)

    return counting

def stable_states_vs_c_vs_N(
        bits, min, dim=2, noise=0, samples=1e7, stable=False, binary=False,
        normalize=False, experiment='counting', extension='.dat', counting={},
        y_scale='linear', save_=True):

    prefix = data_dir + experiment
    if stable:
        suffix = 'stable-n%.2g' %(samples)
    elif binary:
        suffix = 'binary-n%.2g' %(samples)
    else:
        suffix = 'random-n%.2g' %(samples)

    if '%d' %min not in counting.keys():
        counting['%d' %min] = {}

    fig = figure()

    for bit in bits:

        #densities = arange(1, bit+1.) / bit
        n_zeros = arange(bit, dtype = float)
        densities = (c-n_zeros/bit)[::-1]
        G = array([float(sizes[str(bit)][density]['G'])
                   for density in map(str, densities)])
        G_red = array([float(sizes[str(bit)][density]['G_red'])
                       for density in map(str, densities)])

        y = []
        for density in densities:
            filename = get_filename(bit, density, min, dim, noise)
            try:
                data = cPickle.load(
                    open(prefix + filename + suffix + extension))
            except: #IOError ?
                n = loadtxt(prefix + filename + suffix + '.tsv', int)[-1][-1]
            else:
                n = data['%d' %min]['%d' %dim]['%d' %bit]
            y.append(n)

        if normalize:
            y  = array(y, dtype=float)
            # state space size, ie, total number of possible phenotypes
            L = float(sizes[str(bit)][str(dim)]['L'])
            if normalize == 'L':
                y /= L
            if normalize == 'G':
                y /= G
            if normalize == 'LG':
                y /= L * G

        plot(densities, y, 'o--', label='N = %d' %bit)

        counting['%d' %min]['%d' %bit] = y

    xticks(densities)
    xlabel('c')
    experiment = 'stable_states'
    ylabel(experiment)
    if normalize:
        suffix += '-%s' %normalize
    title('K = %s, min = %d, ' %(dim, min) + suffix)
    leg = legend(loc=0)
    yscale(y_scale)

    if save_:
        prefix = figs_dir + experiment
        filename = '_vs_c_vs_N-k%s_min_%s-noise%s-' %(dim, [min], noise)
        extension = format
        savefig(prefix + filename + suffix + extension)

    return counting

def stable_states_vs_k(
        dims, min, bits=[4], n_zeros=0, noise=0, samples=1e6,
        full_enum=True, symm=False, stable=False, binary=False,
        experiment='counting', extension='.dat', potential=False,
        save_=False, linestyle='o--'):

    densities = unique(c-n_zeros/array(bits, dtype=float))
    density = densities[0]
    bit = bits[0]

    prefix = data_dir + experiment
    if full_enum:
        suffix = 'binary-enum'
        if symm:
            suffix += '-symm'
        else:
            suffix += '-full'
    elif stable:
        suffix = 'stable-n%.2g' %(samples)
    elif binary:
        suffix = 'binary-n%.2g' %(samples)
    else:
        suffix = 'random-n%.2g' %(samples)

    y = []
    for dim in dims:

        filename = get_filename(bit, density, min, dim, noise)
        try:
            data = cPickle.load(open(prefix + filename + suffix + extension))
        except: #IOError ?
            n = loadtxt(prefix + filename + suffix + '.tsv', int)[-1][-1]
        else:
            n = data['%d' %min]['%d' %dim]['%d' %bit]
        y.append(n)

    fig = figure()
    if potential:
        plot(dims, array(dims)**bit, 'k--', label='potential')
        linestyle = 'ko'

    plot(dims, y, linestyle, label = suffix)
    xlabel('k')
    experiment = 'stable_states'
    ylabel(experiment)
    title('N = %s, min = %d, ' %(bit, min) + suffix)
    leg = legend(loc=0)

    if save_:
        prefix = figs_dir + experiment
        filename = '_vs_k-N%s_min_%s-noise%s-' %(bit, [min], noise)
        extension = format
        savefig(prefix + filename + suffix + extension)

    return y


def new_stable_states_vs_k(
        myrun_fname, dims=None, i=0, width=0,
        fig=None, color='b', normed=False, marginal=False,
        experiment='dims', bits=None, noises=None, potential=False,
        n=None, samples=[1e8, 1e6], devo_times=[16, 32, 100, 200],
        plotting_=True, ls='o', mec=None, jitter=0,
        x_label='', title_='', label='', noise_label='noise',
        y_scale='linear', fig_name='', verbose=False):

    potential_label = 'potential'

    # load myrun file with all parameters
    myrun = cPickle.load(open(logs_dir + myrun_fname + '.run'))

    if experiment == 'dims':
        bits = int(myrun_fname[1])
        if not any(dims):
            dims = arange(2, 10)
        x = dims
        # get all filename variations
        fnames = [get_replaced_filename(myrun=myrun, dim=dim)
                  for dim in dims]
        if not x_label:
            x_label = '# gene expression states, k'
            # '# intermediary gene expression states (2 = on or off)'

    elif experiment == 'bits':
        dims = 2
        if not any(bits):
            bits = arange(4, 20)
        x = bits
        # get all filename variations
        fnames = [get_replaced_filename(myrun=myrun, bits=bit)
                  for bit in bits]
        if not x_label:
            x_label = 'network size N'

    elif experiment == 'noise':
        dims = 2
        bits = 4
        if not any(noises):
            noises = [1e-100, .01, .05, .1, .15, .2, .25]
        x = noises
        # get all filename variations
        fnames = [get_replaced_filename(myrun=myrun, noise=noise)
                  for noise in noises]
        if not x_label:
            x_label = 'noise'
        potential = False

    if not label:
        label = get_myrun_label(myrun, noise_label)

    if myrun.experiment == 'phenotypes':
        func_load = load_stable_states_from_txt_file
        y_label = '# phenotypes' # '# stable states'

    elif myrun.experiment == 'stability':
        func_load = load_stability_from_txt_file
        y_label = 'stability'

    # load data
    y = [func_load(fname + '.tsv', n=n, samples=samples, devo_times=devo_times,
                   verbose=verbose)
         for fname in fnames]

    if normed:
        y /= (dims**bits).astype(float)
        y_label = '% accessible phenotypes'

    elif marginal:
        x = x[1:]
        y = array(y)
        y = y[1:]/y[:-1]
        y_label = 'marginal contribution'

    if plotting_:
        if not fig:
            fig = figure()

        if potential:
            plot(x, dims**bits, 'k--', label=potential_label)

        if width:
            bar(x + i*width, y, width, color=color, label=label)
        else:
            plot(fixed_jitter(x, jitter), y, ls, mec=mec, color=color, label=label)

        if not title_ and not i:
            title_ = get_myrun_title(myrun, experiment)

        set_plot_kwargs(None, x_label, y_label, title_, y_scale=y_scale)

        if fig_name:
            save_fig(fig_name)

    return fig, array(y).astype(int)


def stable_states_vs_k_vs_all(
        fnames, n=None, devo_times=[16, 32, 100, 200], samples=[1e8, 1e6],
        colormap=None, dims=None, potential=True, normed=False, marginal=False,
        experiment='dims', bits=None, labels=None, noise_label='noise',
        ls='o', y_scale='linear', title_='', journal='plos',
        x_label='', y_label='', pot_lbl='potential', fig=None, fig_name='',
        verbose=False, save_=True):

    if not labels:
        labels = ['']*len(fnames)

    if colormap:
        if type(colormap) is dict:
            colormap = array(colormap[len(fnames)])/255.#colors[journal]
    else:
        colormap = array(colorbrewer.OrRd[len(fnames)])/255.#colors[journal]

    # if width == 0: scatter plot else bar plot in new_stable_states_vs_k()
    if normed or marginal:
        width = 0
    elif len(fnames) < 5:
        width = 0
    elif len(fnames) == 5:
        width = .15
    elif len(fnames) == 7:
        width = .11

    bar_xoffset = len(fnames)/2.
    hline_xoffset = .4
    text_yoffset = 7

    if experiment == 'dims':
        bits = int(fnames[0][1])
        if bits == 1:
            bits = 10.
        if not any(dims):
            dims = arange(2, 10)
        x = dims
        if not x_label:
            x_label = '# gene expression states, k'
            # '# intermediary gene expression states (2 = on or off)'
    else:
        dims = 2
        if not any(bits):
            bits = arange(4, 20)
        x = bits
        if not x_label:
            x_label = 'network size, N'

    if isinstance(bits, int):
        if bits == 2:
            text_yoffset = .1
        elif bits == 3:
            text_yoffset = 1

    if not fig:
        fig = figure()

    data = array([new_stable_states_vs_k(
        fname, dims, i - bar_xoffset, width, fig, color, normed, marginal,
        experiment, bits, n=n, samples=samples, devo_times=devo_times,
        ls=ls, label=label, noise_label=noise_label, verbose=verbose)
                  for i, (fname, color, label)
                  in enumerate(zip(fnames, colormap, labels))])

    fig = data[0][0]
    data = data.T[1]
    heights = array([array(y) for y in data])

    axis_ = list(axis())
    if isinstance(bits, int) and bits == 2:
        axis_[-1] += 7

    if potential:
        y = dims**bits
        if width:
            hlines(
                y, x - hline_xoffset, x + hline_xoffset,
                color='k', linestyles='--', label=pot_lbl)
            [text(w, h + text_yoffset, '%d'%h, ha='center', va='bottom')
             for w, h in zip(x, y) if h < axis_[-1]]
        else:
            plot(x, y, 'k--', label=pot_lbl)

    # set_plot_kwargs
    myrun = cPickle.load(open(logs_dir + fnames[0] + '.run'))

    if not title_:
        title_ = get_myrun_title(myrun)

    if not y_label:

        if myrun.experiment == 'phenotypes':
            if normed:
                y_label = '% visible phenotypes'
            elif marginal:
                y_label = 'marginal contribution'
            else:
                y_label = '# phenotypes' # '# stable states'
        else:
            y_label = 'stability'

    set_plot_kwargs(x, x_label, y_label, title_, axis_, y_scale=y_scale)

    # save
    if save_:
        if not fig_name:
            fig_name = get_noise_figname(fnames[0])
            if normed:
                fig_name = '-'.join((fig_name, 'visible'))
            elif marginal:
                fig_name = '-'.join((fig_name, 'marginal'))
        save_fig(fig_name)

    return fig, heights


def stable_states_vs_k_vs_meta(meta_myruns, color=None, save_=True):
    return [stable_states_vs_k_vs_all(fnames, colormap=color, save_=save_)
            for fnames in meta_myruns]

def stable_states_vs_N(
        bits=arange(4,20), fname='N4bun8inf', n=None,
        samples=[1e8, 1e6], devo_times=[16, 32, 100, 200], fig=None,
        color='b', potential=True, label='visible', noise_label='noise',
        ls='o', mec=None, jitter=0, fig_name='stable_states_vs_N',
        verbose=False, y_scale='log'):

    return new_stable_states_vs_k(
        fname, fig=fig, color=color, experiment='bits', bits=bits,
        potential=potential, n=n,
        samples=samples, devo_times=devo_times,
        ls=ls, mec=mec, jitter=jitter,
        label=label, noise_label=noise_label,
        y_scale=y_scale, fig_name=fig_name, verbose=verbose)

# scatter plot
def stable_states_vs_N_vs_noise(
        fnames=myrun_noiseN4ux, n=None, samples=[], devo_times=[],
        fig=None, ls='o-', mec=None, jitter=0.3,
        pot_lbl='potential',  noise_lbl='noise', samp_lbl='samples',
        experiment='transition', fig_name='', save_=False, verbose=False):

    bits = arange(4,20)

    myrun = cPickle.load(open(logs_dir + fnames[-1] + '.run'))

    if myrun.noise and len(fnames) > 2:
        colormap = array(colorbrewer.OrRd[len(fnames)])/255.
    else:
        colormap = colors[None]

    if not fig:
        fig = figure()

    # potential
    plot(bits, 2**bits, 'k--', label=pot_lbl)

    # sample size
    if n:
        xmin, xmax, _, _ = axis()
        # true sample size is only approximately n
        x = get_int_datapoints(0, log10(1e8), 100)
        y = find_nearest(x, n)
        hlines(y, xmin, xmax, 'k', '-.', label=samp_lbl)

    for fname, color in zip(fnames, colormap):
        fig, y = stable_states_vs_N(
            bits, fname, n, samples, devo_times,
            fig, color, False, '', noise_lbl, ls, mec, jitter, '', verbose)

    # set title
    if len(devo_times) == 1:
        myrun.devo_time = devo_times[0]
    title(get_myrun_title(myrun, experiment, n, False))

    if save_:
        fig_name = get_transition_fig_name(myrun, n, fig_name, experiment)
        save_fig(fig_name, get_figs_save_dir('', experiment))

    return fig

# scatter plot
def stable_states_vs_k_vs_noise(
        fnames, n=None, samples=[], devo_times=[],
        fig=None, ls='o-', mec=None, jitter=0.3, y_scale='log', hline=False,
        pot_lbl='potential',  noise_lbl='noise', samp_lbl='samples',
        experiment='transition', fig_name='', save_=False, verbose=False):

    dims = arange(2,10)

    myrun = cPickle.load(open(logs_dir + fnames[-1] + '.run'))

    if myrun.noise and len(fnames) > 2:
        colormap = array(colorbrewer.OrRd[len(fnames)])/255.
    else:
        colormap = colors[None]

    if not fig:
        fig = figure()

    # potential
    plot(dims, dims**myrun.bits, 'k--', label=pot_lbl)

    # sample size
    if n:
        n = 10**n
        if hline:
            xmin, xmax, _, _ = axis()
            # true sample size is only approximately n
            x = get_int_datapoints(0, log10(1e8), 100)
            y = find_nearest(x, n)
            hlines(y, xmin, xmax, 'k', '-.', label=samp_lbl)

    for fname, color in zip(fnames, colormap):
        fig, y = new_stable_states_vs_k(
            fname, fig=fig, color=color, n=n,
            samples=samples, devo_times=devo_times,
            ls=ls, mec=mec, jitter=jitter, noise_label=noise_lbl,
            y_scale=y_scale, verbose=verbose)

    # set title
    if len(devo_times) == 1:
        myrun.devo_time = devo_times[0]
    title(get_myrun_title(myrun, experiment, n, False))

    if save_:
        fig_name = get_transition_fig_name(myrun, n, fig_name, experiment)
        save_fig(fig_name, get_figs_save_dir('', experiment))

    return fig

# bar chart
def stable_states_vs_N_vs_noise_bar(
        bits=arange(4,11), fig_name='stable_states_vs_N_vs_noise'):
    return stable_states_vs_k_vs_all(
        myrun_noiseN4s, colormap=colorbrewer.OrRd,
        experiment='bits', bits=bits, fig_name=fig_name)

def stability_vs_k_vs_noise():
    return stable_states_vs_k_vs_all(
        myrun_stabilityN4, colormap=colorbrewer.OrRd, potential=False)

def percent_visible_states(
        fnames, color=None, x_label='', labels=[''], title_='', ls='--o',
        fig_name='percent_visible_states', save_=False):
    return stable_states_vs_k_vs_all(
        fnames, colormap=color, potential=False, normed=True, labels=labels,
        title_=title_, ls=ls, x_label=x_label, fig_name=fig_name, save_=save_)

def marginal_visible_states(fnames, color=None, save_=False):
    return stable_states_vs_k_vs_all(
        fnames, colormap=color, potential=False, marginal=True, save_=save_)

def plot_access_vs_samples_vs_noise(
        fnames=myrun_noiseN4s, k=4, devo_times=None, hline=True, fig=None,
        verbose=False):

    if not fig:
        fig = figure()

    colormap = array(colorbrewer.OrRd[len(fnames)])/255.

    for fname, color in zip(fnames, colormap):
        # load myrun file with filename and parameters
        myrun = cPickle.load(open(logs_dir + fname + '.run'))
        # replace myrun.dim with k in filename
        filename = get_replaced_filename(myrun=fname, dim=k)
        # get all filenames with parameter substitution
        filenames = get_replaced_suffix(filename, [], devo_times)
        # load all runs
        #return [loadtxt(filename + '.tsv', int) for filename in filenames]
        data = []
        for filename in filenames:
            _data = load_file(filename + '.tsv', verbose, lambda x: loadtxt(x, int))
            if any(_data):
                data.append(_data)
        # concatenate and mask invalid (there's actually no NaNs, but see below)
        data = concatenate(data)

        # mask the first point of runs #2 and #3 to avoid connecting straight lines
        x = data.T[0]
        y = data.T[1]
        x = ma.masked_where(x==1, x)
        x.mask[0] = False
        y = ma.masked_where(x.mask, y)

        label = get_myrun_label(myrun, '')
        plot(x, y, c=color, lw=1.5, label=label)

        # print sample size
        if verbose:
            print data.T[0][-1]

    # set_plot_kwargs
    xscale('log')
    xlabel(r'$\lambda$')
    ylabel(r'$\mathcal{U}$')
    title('stable' if myrun.stable else 'random')
    legend(loc='best', fancybox=True, title=r'$\sigma$')
    gca().set_axis_bgcolor("#E5E5E5")

    # dashed lines with parameters/limits (e.g. sample size)
    if hline:
        ## hlines and vlines
        # max visible states (state space size)
        #myrun.dim = k
        umax = k**myrun.bits
        xmin, xmax, _, _ = axis()
        hlines(umax, xmin, xmax, 'k', '-.')
        label = r'$\mathcal{U}_{max} = k^N$'
        annotate(label, (xmax,umax), xytext=(-65, 2),
                 textcoords='offset points', va='bottom')

        # max number of matrices (genotype size)
        wmax = 2**myrun.bits**2
        _, _, ymin, ymax = axis()
        vlines(wmax, ymin, umax, 'k', ':')
        label = r'$2^{N^2}$'
        annotate(label, (wmax,umax), xytext=(2, -30),
                 textcoords='offset points', ha='left')

        # max discovery rate (u==lambda)
        x = y = [1, umax]
        plot(x, y, 'k--')
        label = r'$\mathcal{U} = \lambda$'
        annotate(label, (umax,umax), xytext=(-12, -30),
                 textcoords='offset points', ha='right')

    return fig

def stable_states_overrepresentation(
        dim=3, fnames=('N4bso1e-100n8200'),
        noises=[1e-100, .01, .05, .1, .15, .2, .25],
        fig=None, save_=True):

    experiment = 'overrepresentation'

    myrun = cPickle.load(open(logs_dir + fnames[0] + '.run'))
    bits = myrun.bits

    L = sizes[str(bits)][str(dim)]['L']

    base, basename = get_base(dim, bits)[:2]
    delta = base[-1] - base[-2]
    bins = append(base, 1 + delta)

    phenotypes = generate_full_phenotype_space(noises[0], basename, base, bits)
    phenotypes = array(list(phenotypes))

    filename = get_replaced_filename(myrun=myrun, dim=dim)
    if len(fnames) == 1:
        filenames = [get_replaced_filename(filename, myrun, noise=noise)
                     for noise in noises]
    else:
        filenames = [get_replaced_filename(myrun=fname, dim=dim)
                     for fname in fnames]
    data = (loadtxt(filename + '.txt', int) for filename in filenames)

    # frequency
    #vectors = array([get_vector(x, basename, base, bits) for x in data])
    #fig = figure();hist(vectors.flatten(), bins)

    n = [histogram(x, arange(L+1))[0] for x in data]
    masks = [x == 0 for x in n]

    # count the number of occurrences of each state in each phenotype vector
    n = [[unique([count_nonzero(v == s)
                  for v in phenotypes[-mask]]) for s in base] for mask in masks]

    y = [histogram(phenotypes[-mask].flatten(), bins)[0] for mask in masks]
    potential = histogram(phenotypes.flatten(), bins)[0]

    n_bars = len(noises) if len(noises) > 1 else len(fnames)
    if n_bars == 1:
        width = .4
    elif n_bars == 5:
        width = .15
    elif n_bars == 7:
        width = .1*3/dim
    bar_xoffset = n_bars/2.
    colors = array(colorbrewer.OrRd[n_bars])/255.#colors[journal]
    labels = ['']*n_bars#[get_noise_label(noise) for noise in noises]

    x = base
    heights = y/potential.astype(float)

    if not fig:
        fig = figure()

    [bar(x + (i-bar_xoffset)*width, y, width, color=color, label=label)
     for i, (y, color, label) in enumerate(zip(heights, colors, labels))]

    x_label = r'gene expression state, $s_i$'
    y_label = '% visible'
    myrun.dim = dim
    title_ = get_myrun_title(myrun, '-'.join(('dims', experiment)))
    set_plot_kwargs(x, x_label, y_label, title_)

    if save_:
        fig_name = [get_noise_figname(fname) for fname in fnames]
        fig_name = '-'.join(fig_name)
        fig_name = ('-k%d-'%dim).join((fig_name, experiment))
        save_fig(fig_name)

    return fig, heights, masks, n

def plot_noise_functions(
        noise_function, dims, samples=1e2, bits=4,
        distance_function=int_vector_hamming_distance,
        plotting_=True, save_=True):

    y = []
    for dim in dims:
        base, basename, jumps = get_base(dim, bits)
        phenotypes = [generate_vector_phenotype(base, bits)
                      for i in range(int(samples))]

        y.append([
            mean([
                [distance_function(phenotype, noise_function(
                    phenotype.copy(), noise_size, base=base))
                 for i in range(int(samples))]
                for phenotype in phenotypes])
            for noise_size in range(bits+1)])

    if plotting_:
        fig = figure()
        lines = plot(dims, y)
        set_plot_kwargs(
            None, 'k', distance_function.__name__, noise_function.__name__,
            labels=range(bits+1), _grid=True)

    if save_:
        fig_name = '-'.join((noise_function.__name__,
                             distance_function.__name__))
        save_fig(fig_name)

    return fig, y


############################################
## PERIOD or CYCLE SIZE or ATRACTOR LENGTH #
############################################

#
def plot_period_dist(
        dim=2, y_scale='log', color='', experiment='period_dist',
        save_dir=get_figs_save_dir(None), save=True):

    fig = figure()

    #for min in mins:
    min = -1.
    filename = 'period-N10c1_%s-n1e+06-uniform.dat' %get_base(dim, min=min)[1]
    period   = array(cPickle.load(open(data_dir + filename)))

    # even period
    i = 1
    bins = arange(period.min() + i, period.max() + 2 - i, 2)
    even = array([p for p in period if not p%2])
    n, x = histogram(even, bins, normed=True)
    marker = markers[0]
    label = 'off $= %d$, even $l$' %min
    plot(x[:-1], n, color + marker, label=label)

    #odd
    i = 0
    bins = arange(period.min() + i, period.max() + 2 - i, 2)
    odd = array([p for p in period if p%2])
    n, x = histogram(odd, bins, normed=True)
    marker = markers[1]
    label = 'off $= %d$,  odd $l$' %min
    plot(x[:-1], n, color + marker, label=label)

    min = 0.
    filename = 'period-N10c1_%s-n1e+06-uniform.dat' %get_base(dim, min=min)[1]
    period   = array(cPickle.load(open(data_dir + filename)))
    bins     = arange(period.min(), period.max() + 2)
    n, x     = histogram(period, bins, normed=True)
    marker   = markers[2]
    label    = 'off $= \,\, %d$,    all $l$' %min
    plot(x[:-1], n, color + marker, label = label)

    polyfit_regression(x[:-1], n, x_label='l', color=color)

    set_plot_kwargs(
        None, 'attractor length $l$', 'probability', y_scale=y_scale)

    if save:
        save_fig(experiment, save_dir)

    return fig



################
## PATH LENGTH #
################

#
def plot_path_length_vs_N(
        density=c, degree=None, max_bit=None, min_bit=4, dim=2, regression=True,
        x_scale='linear', y_scale='log', experiment='path_leng_vs_N',
        extension='.txt', plotting_=True, save_=True):

    if save_:
        regression = True

    # get x
    key = degree if degree else density
    bits = path_bits[key]

    #load data
    if not degree:
        lengths = [loadtxt(data_dir + experiment + '%d'%bit + extension)
                   for bit in bits]
        for i, bit in enumerate(bits):
            filename = data_dir + experiment + '%d-n1e+04'%bit + extension
            try:
                data = loadtxt(filename)
            except: #IOError ?
                pass
            else:
                lengths[i] = append(lengths[i], data)
    else:
        density = None
        lengths = [loadtxt(
            data_dir + experiment + '%dk%d' %(bit, degree) + extension, int)
                   for bit in bits]

    # stats
    y = array([data.mean() for data in lengths])
    z = array([data.std()  for data in lengths])
    w = array([data.max()  for data in lengths])

    # powers of two
    powers_w = 2**ceil(log2(w)).astype(int) # with data

    # plot
    if plotting_:
        fig = figure()
        errorbar(bits, y, z, c='b', label = 'mean + std')
        plot(bits, w, 'ko', label = 'max')
        plot(bits, powers_w, 'r--', label = '$2^n$') # with data
        axis_ = [bits[0], bits[-1], y[0] - z[0], powers_w[-1]]

    # linear regression in log scale (exponential fit)
    if regression:
        if not max_bit: max_bit = max_bits[experiment][key]
        xx, yy, wy, powers_wy = polyfit_regression(
            bits, y, w, min_bit, max_bit, [fit_functions[key][0]],
            [fit_functions[key][1]], plotting_=plotting_)
        axis_ = [min_bit, max_bit, yy[0], powers_wy[-1]]

    if plotting_:
        set_plot_kwargs(
            None, 'N', 'path length to equilibrium', '', axis_, x_scale,
            y_scale)

    # save path length, plot and dict
    if save_:

        # save path_length_vs_N
        prefix = data_dir + experiment + print_min_max(bits, bits)
        filename = '-c%g' %density if density else '-k%d' %degree
        savetxt(prefix + filename + '.csv', (bits, y, z), delimiter=',')

        # save plot
        if plotting_:
            filename  = '%d_%d' %(bits[0], bits[-1])
            filename += '-c%g' %density if density else '-k%s' %degree
            save_plot(filename, experiment)

        # load dict
        devo_times = load_file(maps_dir + 'devo_times.dict')

        # initialize empty dict
        if not devo_times:
            devo_times = dict(mean = {'%d'%dim: {}}, maximum = {'%d'%dim: {}})

        if key not in devo_times['mean']['%d' %dim]:
            devo_times['mean']['%d' %dim]['%s' %key] = {}
        if key not in devo_times['maximum']['%d' %dim]:
            devo_times['maximum']['%d' %dim]['%s' %key] = {}

        # mean
        yy = ceil(yy).astype(int)

        # conservative limit = max(powers_w, powers_wy)
        powers = zeros_like(powers_wy) # objects have to have same shape
        powers[:powers_w.size] = powers_w
        powers = maximum(powers, powers_wy)

        # save dict
        for i, x in enumerate(xx):
            devo_times['mean']   ['%d' %dim]['%s' %key]['%d' %x] = yy[i]
            devo_times['maximum']['%d' %dim]['%s' %key]['%d' %x] = powers[i]
        save_file(maps_dir + 'devo_times.dict', devo_times)

    return y, z, w, powers_w

#
def plot_path_length_vs_N_vs_k(
        densities=[None, None, c], degrees=[2, 4, None],
        markers=markers, linestyles=linestyles, color='',
        axis_=[3, 1500, 0, 2e4], x_scale='log', y_scale='log',
        experiment='path_leng_vs_N', extension='.csv', format=format,
        bbox_inches=bbox, pad_inches=pad_inches, fig_prefix=fig_prefixes[None],
        plotting_=True, error=False, save_=True):

    if plotting_:
        fig = figure()

    loop = zip(densities, degrees, markers, linestyles)
    for density, degree, symbol, linestyle in loop:

        # get x
        key = degree if degree else density
        bits = path_bits[key]

        # filename
        prefix = data_dir + experiment + print_min_max(bits, bits)
        filename = '-c%g' %density if density else '-k%d' %degree

        # get y (mean) and z (std)
        try:
            bits, y, z = loadtxt(prefix + filename + extension, delimiter=',')
        except: #IOError ?
            y, z = plot_path_length_vs_N(
                bits, density, degree, regression=False, plotting_=False,
                save_=False)[:2]
            if save_:
                savetxt(
                    prefix + filename + extension, (bits, y, z), delimiter=',')

        # regression
        x_fit = fit_functions[key][0]
        x, yy, ya, yb = polyfit_regression(
            bits, y, None, bits[0], bits[-1], [x_fit], plotting_ = False)

        if plotting_:
            if not degree:
                degree = 'N'
            label = 'K = %s' %degree
            # mean +- std vs N
            if error:
                errorbar(
                    bits, y, z, fmt = color + symbol, label=label)
                    #mfc = color)
            # just mean  vs N
            else:
                plot(bits, y, color + symbol, label=label)
            # regression vs N
            plot(
                x, yy, color + linestyle,
                label=get_regression_label(x_fit, ya, yb))

    if plotting_:
        set_plot_kwargs(
            None, 'N', 'path length to equilibrium', '', axis_, x_scale,
            y_scale)
        if save_:
            fname = '.'.join((experiment, format))
            fname = '-'.join((fig_prefix, fname))
            savefig(
                report_dir + fname, dpi=dpi, format=format,
                bbox_inches=bbox, pad_inches=pad_inches)

    return bits, y, z, yy

#
def plot_power_vs_mu(bits, filename='power_vs_mu_N', extension='.txt'):
    #bits = range(4,13) + [13, 15, 20, 25, 30, 35, 40]
    powers = [loadtxt(data_dir + filename + '%d'%bit + extension)
              for bit in bits]
    x = [data.T[1] for data in powers]
    y = [data.T[0] for data in powers]
    fig = figure()
    '''x = list(flatten(x))
    y = list(flatten(y))
    plot(x, y, 'k.')
    '''
    for xx, yy in zip(x,y):
        plot(xx, yy, 'k.')

    #plot(x, 2*array(x), 'k--')

    xlabel('path length to equilibrium (mu in brent)')
    ylabel('power of two (in brent)')


###############################
## SAMPLES AND DISCOVERY TIME #
###############################

#
def discovery_time_vs_c_vs_N(
        bits, min, dim=2, noise=0, samples=1e7, stable=False, binary=False,
        normalize=False, experiment='counting', extension='.tsv', sampled={},
        y_scale='linear', axiss=None, save_=True):

    prefix = data_dir + experiment
    if   stable:
        suffix = 'stable'
    elif binary:
        suffix = 'binary'
    else:
        suffix = 'random'

    if '%d' %min not in sampled.keys():
        sampled['%d' %min] = {}

    fig = figure()

    for bit in bits:

        #densities = arange(1, bit+1.) / bit
        n_zeros = arange(bit, dtype = float)
        densities = (c-n_zeros/bit)[::-1]
        G = array([float(sizes[str(bit)][density]['G'])
                   for density in map(str, densities)])
        G_red = array([float(sizes[str(bit)][density]['G_red'])
                       for density in map(str, densities)])

        y = []
        for density in densities:
            filename = get_filename(bit, density, min, dim, noise)
            try:
                n = loadtxt(prefix + filename + suffix + '-n%.2g' %(samples)
                            + extension, int)[-1][0]
            except: #IOError ?
                n = loadtxt(prefix + filename + suffix + '-n%.2g' %(samples/10)
                            + extension, int)[-1][0]
            y.append(n)

        if normalize:
            y  = array(y, dtype = float)
            # state space size, ie, total number of possible phenotypes
            L = float(sizes[str(bit)][str(dim)]['L'])
            if normalize == 'L':
                y /= L
            if normalize == 'G':
                y /= G
            if normalize == 'LG':
                y /= L * G

        plot(densities, y, 'o--', label='N = %d' %bit)

        sampled['%d' %min]['%d' %bit] = y

    xticks(densities)
    xlabel('c')
    experiment = 'discovery_time'
    ylabel(experiment)
    if normalize: suffix += '-%s' %normalize
    title('K = %s, min = %d, ' %(dim, min) + suffix)
    leg = legend(loc=0)
    for t in leg.get_texts():
        t.set_fontsize('small') # the legend text fontsize
    yscale(y_scale)
    if axiss:
        axis(axiss)

    if save_:
        prefix = figs_dir + experiment
        filename = '_vs_c_vs_N-k%s_min_%s-noise%s-' %(dim, [min], noise)
        extension = format
        savefig(prefix + filename + suffix + extension)

    return sampled

#
def samples_vs_N(
        bits, density=c, degree=None, n_zeros=None, min=-1., dim=2,
        noise=0, min_n=min_N, sample_sizes=[1e6], devo_times=[inf],
        experiment='stability', extension='.dat', verbose=False):

    samples_vs_N = []
    for i, bit in enumerate(bits):

        n_samples = 0
        # filename_exceptions
        precision = density_precisions[bit] if bit in density_precisions else 2

        # degree and density
        if degree:
            density  = None
        if verbose:
            print density, n_zeros, degree, bit
        density = get_density_degree_nzeros(density, n_zeros, degree, bit)[0]
        if verbose:
            print density, n_zeros, degree, bit

        # load data
        prefix  = data_dir + experiment
        for samples in sample_sizes:
            for devo_time in devo_times:
                suffix = get_suffix(samples, devo_time=devo_time)
                if verbose: print suffix
                n_samples = max(
                    n_samples, load_samples(
                        bit, samples, density, min, dim, noise, suffix, prefix,
                        extension, precision, verbose))

        samples_vs_N.append(n_samples)

    return samples_vs_N

#
def plot_samples_vs_N_vs_k(
        bits=None, densities=[None, None, c], degrees=[2, 4, None],
        sample_sizes=[[1e8], [1e8], [1e8, 1e7, 1e6]],
        devo_times=[inf, 'maximum', 'mean'],
        filter_=False, color='', linestyles=markers,
        axis_=[3, 1.1e4, 1, 2e8], x_scale='log', y_scale='log',
        experiment='samples_vs_N', format=format,
        bbox_inches=bbox, pad_inches = pad_inches,
        fig_prefix=fig_prefixes[None],
        save_dir=get_figs_save_dir(None), save_=True):

    fig = figure()

    if not bits:
        bits = stability_bits[2]

    loop = zip(densities, degrees, sample_sizes, linestyles)
    for density, degree, sample_size, linestyle in loop:
        key = degree if degree else density
        bit = bits[:where(bits == max_bits['path_leng_vs_N'][key])[0] + 1]

        y = samples_vs_N(
            bit, density, degree, sample_sizes=sample_size,
            devo_times=devo_times)

        if filter_:
            yy = unique(y)[::-1]
            x = [bit[where(y == z)[0].max()] for z in yy]
            xy = [(x1, y1) for x1, x2, y1 in zip(x[1:], x, yy[1:]) if x1 > x2]
            if xy[1][0] < xy[0][0]:
                xy.pop(1)
            x0 = array([x for i, x in enumerate(bit) if y[i] == yy[0]])
            xx = append(x0, transpose(xy)[0])
            yy = append([yy[0]]*(x0.size), transpose(xy)[1])
        else:
            xx = bit
            yy = y

        if not degree:
            degree = 'N'
        plot(xx, yy, color + linestyle, label = 'K = %s' %degree)

    set_plot_kwargs(None, 'N', '# samples', '', axis_, x_scale, y_scale)

    if save_:
        save_fig(experiment, save_dir, fig_prefix, format=format,
                 bbox_inches=bbox_inches, pad_inches=pad_inches)

    return fig


##############
## VIABILITY #
##############

def viability_vs_N(
        bits=None, density=c, degree=None, samples=1e5, stable=True,
        binary=False, devo_time=inf, n_mutants=100, all_mutants=1,
        deletion=False, change_sign=False, min_n=min_N, min_samples=min_Samples,
        color='', marker=markers[0], linestyle=linestyles[0], axis_=None,
        x_scale='linear', label='', plotting_=False, verbose=False):

    if change_sign:
        binary = True
    suffix = get_suffix(
        samples, stable=stable, binary=binary, devo_time=devo_time,
        n_mutants=n_mutants, all_mutants=all_mutants, deletion=deletion,
        change_sign=change_sign)
    key = degree if degree else density
    if not any(bits):
        bits = viability_bits[key]
    if key == 3 and (deletion or change_sign):
        bits = bits[:where(bits == max_bits['viability_vs_N'][key])[0] + 1]

    viability = []
    key = 'degree' if degree else 'density'
    for bit in bits:
        if key == 'degree':
            density = None
        if verbose:
            print density, degree, bit
        density, degree = get_density_degree_nzeros(
            density, None, degree, bit)[::2]
        n  = load_viability(
            suffix, bit, density, min_n=min_n, min_samples=min_samples,
            verbose=verbose)
        n /= degree*bit if all_mutants == 1 else n_mutants
        viability.append(n)

    # removes y zeros
    index = nonzero(viability)
    x = array(bits)[index]
    y = array(viability)[index]

    if plotting_:
        fig = figure()
        plot(x, y, color + marker + linestyle, label=label)
        set_plot_kwargs(None, 'N', 'viability', '', axis_, x_scale)

    return x, y

def plot_viability_vs_N(
        bits=None, density=c, degree=None, mutations=mutations,
        n_mutants=100, min_n=min_N, min_samples=min_Samples, stable=True,
        binary=False, color='', markers=markers, linestyle=linestyles[0],
        axis_=None, x_scale='linear', experiment='viability_vs_N',
        extension='.csv', save_dir=get_figs_save_dir(None),
        fig_prefix=fig_prefixes[None], plotting_=True, save_=True,
        verbose=False):

    key = degree if degree else density
    if not any(bits):
        bits = viability_bits[key]

    if plotting_:
        fig = figure()
    if save_:
        prefix = data_dir + experiment
        filename = '-c%g' %density if density else '-k%d' %degree
        if stable:
            filename += '-stable'

    if plotting_ or save_:
        for key, value, marker in zip(
                mutations.keys(), mutations.values(), markers):
            x, y = viability_vs_N(
                bits, density, degree, stable=stable, binary=binary,
                n_mutants=n_mutants, all_mutants=value[0],
                deletion=value[1], change_sign=value[2], min_n=min_n,
                min_samples=min_samples, verbose=verbose)
            if any(y):
                if plotting_:
                    plot(x, y, color + marker + linestyle, label = key)
                if save_:
                    savetxt(
                        prefix + '-'.join((filename, key)) + extension, (x, y),
                        delimiter=',')

    if plotting_:
        set_plot_kwargs(None, 'N', 'viability', '', axis_, x_scale)
        if save_:
            save_fig(experiment + filename, save_dir, fig_prefix)

    return fig


##################################
## DATA from PERTURB_AND_ANALYSE #
##################################

def unite(data):
    new = {}
    for g in data[0].keys():
        new[g] = {}
        for stat in data[0][g].keys():
            if stat != 'gain':
                new[g][stat] = []
                for data_i in data:
                    new[g][stat].extend(data_i[g][stat])
            else:
                new[g][stat] = zeros(shape(data[0][g][stat]))
                for data_i in data:
                    new[g][stat] = new[g][stat] + data_i[g][stat]
    del data
    return new

def unite_and_get_data(
        model, generations = G, period = 1, n_neighbours = 3,
        n_runs = range(n_runs), initials = True, mutations = False,
        extension  = '.dat'):

    filename  = get_evolution_filename(
        model, generations, period, n_neighbours, initials, mutations)

    #data = [cPickle.load(open(data_dir + filename + str(run) + extension))
     #       for run in n_runs]

    n = 0
    data = []
    for run in n_runs:
        try:
            file = open(data_dir + filename + print_run(run) + extension)
        except: #IOError ?
            n +=1
        else:
            data.append(cPickle.load(file))

    print '%s samples' %(len(n_runs) - n)
    return unite(data)

def plot_perturbation_stats(
        data, experiment, stats_=pa_stats, extra_stat=extra_stat, color='',
        linestyle=linestyles[0], axis_=None, x_scale='linear',
        fig_width= fig_widths[None][2], fig_height=fig_height*1.2,
        save_dir= get_figs_save_dir(None), fig_prefix=fig_prefixes[None],
        error=False, plotting_=True, save_=True, logging=True,
        verbose=False):

    if plotting_:
        fig = figure(figsize=(fig_width, fig_height))

    x = data.keys()
    x.sort()
    for i, stat in enumerate(stats_):
        if stat != 'gain' and stat != 'hamming':
            y = [mean(data[g][stat]) for g in x]
            z = [std (data[g][stat]) for g in x]

            if plotting_:
                if error:
                    errorbar(x, y, z, label = stat)
                else:
                    plot(x, y, color + markers[i] + linestyle, label = stat)
            if save_:
                savetxt(
                    data_dir + experiment + '-' + stat + '.csv', (x, y, z),
                    delimiter=',')
            if logging:
                savetxt(
                    data_dir + experiment + '-' + stat + '.log' ,
                    map(str,
                        [between_and_within_variation(data[g][stat])
                         for g in x]),
                    fmt = '%s', delimiter = '\n')

    if extra_stat:
        y = []; z = []
        for g in x:
            normed = - array(data[g]['entropy fix']) / log(array(data[g]['stable']))
            normed = nan_to_num(normed)
            y.append(normed.mean())
            z.append(normed.std() )

        if plotting_:
            errorbar(x, y, z, label = extra_stat)
        if save_:
            savetxt(
                data_dir + experiment + '-' + extra_stat + '.csv', (x, y, z),
                delimiter=',')

    if plotting_:
        set_plot_kwargs(
            None, 'generation', '', experiment, axis_=axis_, x_scale=x_scale)
        if save_:
            save_fig(experiment, save_dir, fig_prefix)

    del data
    return fig

def plot_single_stat(
        stat, m3_periods=[1,2,3], generations=1e4, n_neighbours=3,
        initials=True, mutations=False,
        axis_=None, x_scale='log', loc=0, fig_width=fig_widths[None][2],
        save_dir=get_figs_save_dir(None), fig_prefix=fig_prefixes[None],
        error=False, save_=True):

    fig = figure(figsize=(fig_width, fig_height))

    models  = [4] + [3]*len(m3_periods)
    periods = [1] + m3_periods
    for model, period in zip(models, periods):

        filename = get_evolution_filename(
            model, generations, period, n_neighbours, initials, mutations)
        suffix   = '-' + stat
        suffix  += ' fix' if period == 1 else ' cycle'

        x, y, z  = loadtxt(
            data_dir + filename + suffix + '.csv', delimiter = ',')

        label = 'm%sp%s' %(model, period) + suffix

        if error:
            errorbar(x, y, z, label=label)
        else:
            plot(x, y, label=label)

    set_plot_kwargs(
        None, 'generation', stat, '', axis_=axis_, x_scale=x_scale, loc=loc)

    if save_:
        filename = get_evolution_filename(
            model, generations, None, n_neighbours, initials, mutations)
        save_fig(filename + suffix, save_dir, fig_prefix)

    return fig

def plot_perturbation_stats_all_models(
        m3_periods=range(1,6), generations=1e6, n_neighbours=3, n_runs=50,
        stats_=['entropy', 'unique', 'equal'], x_scale='log'):

    models  = [4] + [3]*len(m3_periods)
    periods = [1] + m3_periods

    for initials, mutations in zip([True, False], [False, True]):

        for model, period in zip(models, periods):

            data = unite_and_get_data(
                model, generations, period, n_neighbours, range(n_runs),
                initials, mutations)

            filename = get_evolution_filename(
                model, generations, period, n_neighbours, initials, mutations)

            fig = plot_perturbation_stats(data, filename, x_scale=x_scale)

        for stat in stats_:

            # legend localization
            loc = 0
            if initials:
                if stat == 'equal':
                    loc = 8
                if stat == 'entropy':
                    loc = 4

            fig = plot_single_stat(
                stat, m3_periods, generations, n_neighbours, initials, mutations,
                x_scale=x_scale, loc=loc)


########################################
## STATISTICS: mean, std, median, mode #
########################################

def confidence_interval(
        data, statfunction=ma.mean, alpha=0.05, n_samples=10000,
        output='errorbar'):#'lowhigh'):
    return bootstrap.ci(data, statfunction, alpha, n_samples, output=output)

def compare_std(shortnames, stat):
    all_data = (get_data_corrected_interval(shortname, stat)
                for shortname in shortnames)
    yerr = [between_and_within_variation(data)[-5::2] for data in all_data]
    return array(yerr)

def compare_variation(shortnames, stat):
    all_data = (get_data_corrected_interval(shortname, stat, False, False)
                for shortname in shortnames)
    yerr = [(between_and_within_variation_pop(data, 'mean').mean(),
             between_and_within_variation_pop(data, 'between').mean(),
             between_and_within_variation_pop(data, 'within').mean())
            for data in all_data]
    return array(yerr)

# 2D data
def between_and_within_variation(data):
    return ('mean(axis=0).std', data.mean(axis=0).std(),
            'mean(axis=1).std', data.mean(axis=1).std(),
            'std(axis=0).mean', data.std(axis=0).mean(),
            'std(axis=1).mean', data.std(axis=1).mean(),
            'std', data.std())

# data is a masked array
def between_and_within_variation_pop(data, error='within'):

    ## Average population
    # deviations from the mean, averaged over individuals? shape:
    # (n_pops, time_points, n_individuals) ->
    # (time_points, n_individuals) ->
    # (time_points)
    if error == 'mean':
        return data.mean(axis=0).std(axis=1).data

    # between populations std, averaged over individuals? shape:
    # (n_pops, time_points, n_individuals) ->
    # (time_points, n_individuals) ->
    # (time_points)
    elif error == 'between':
        return data.std(axis=0).mean(axis=1).data

    ## Average individual?
    # within population std, averaged over runs. shape:
    # (n_pops, time_points, n_individuals) ->
    # (n_pops, time_points) ->
    # (time_points)
    elif error == 'within':
        return array([run.std(axis=1).data for run in data]).mean(axis=0)

# mean, median and mode
def get_y_stats(data, statistic='mean'):

    if data.ndim > 2:
        data = ma.masked_invalid(data.swapaxes(0,1).reshape(shape(data)[1], -1))

    # mode
    if statistic == 'mode':
        return stats.mstats.mode(data, axis=1)[0]

    # mean, median
    return getattr(ma, statistic)(data, axis=1)

# mean and std
def get_y_yerr(
        data, stat='', statistic='mean', errors={}, shortname='', d_stat='',
        k=-1):

    if not stat:
        return array([(mean(x), std(x)) for x in data]).T

    # diversity's shape: (runs.size, datapoints.size, -1)
    if stat == 'diversity':
        norm = get_norm(stat, shortname)
        y    = data.T[k].mean(axis=1) / norm
        yerr = data.T[k].std (axis=1) / norm if stat in errors else None

    # conservation, robustness and survivability shape:
    # (runs.size, datapoints.size, bits, bits)
    # other stats' shape:
    # (runs.size, datapoints.size, pop_size)
    else:
        if 'diff' in stat:

            if data.ndim == 3:
                data = abs(data)

            elif data.ndim == 4:
                data = [[abs(difference_diag_nondiag(x, mean))
                         for x in run]
                        for run in data]
                data = ma.masked_invalid(data).reshape(shape(data)[:2] + (1,))

        # 'conservation', 'a.conservation' and 'a.matrix'
        if ('conservation' in stat or
            'matrix' in stat or
            stat == 'survivability' or
            stat == 'robustness'):
            data = data.reshape(shape(data)[:2] + (-1,))
            if stat == 'conservation':
                data = abs(data)

        if stat == 'p' and shortname == 'm4r05b' and shape(data)[1] == 105:
            data = data[:,1:,:] # delete generation 0 (saved by mistake)

        norm = get_norm(stat, shortname)
        #y = data.mean(axis=0).mean(axis=1)
        y = get_y_stats(data, statistic)
        y/= norm
        # y = y.compressed
        if stat in errors:
            yerr = between_and_within_variation_pop(data, errors[stat]) / norm
        else:
            yerr = None

    return y, yerr, data



##############################
## EVOLUTION: AUTOREGULATION #
##############################

def get_pop_function(stats_, experiment='', shortname='', n=None):
    if 'equal'      in stats_:
        # loads analysis = get_pop()
        return pop_robustness
    if 'qp' in stats_ and 'm3r05' in shortname and n == 'all1000':
        # loads get_ind_stat(get_pop())
        return get_qp_pop_matrix
    if 'qp.pop'     in experiment:
        # loads get_non_evolved_pop()
        return get_qp_pop_robustness_and_survivability
    if 'stable.pop' in experiment:
        # loads get_non_evolved_pop()
        return get_qp_pop_robustness_and_survivability
    if 'phenotype'  in stats_:
        return get_stat
    return pop_evolution # loads pop = get_pop()

def get_stat_function(
        stat, data, pop_size=P, bits=N, experiment='', pop_stat_func=None,
        normalize=None, p=None):

    if 'qp.pop' in experiment or 'stable.pop' in experiment:
        data = globals()[pop_stat_func](
            data, stat, pop_size, bits=bits, experiment=experiment, p=p)
        if pop_size > 1e3:
            data = [ma.masked_invalid(data).mean()]
        return data

    if 'robustness' in stat or 'survivability' in stat:
        pop   = data
        size_ = bits**2 if 'qp' in experiment else pop_size
        axis_ =       0 if 'qp' in experiment else 1
        if not (pop and pop.size) or not hasattr(pop, stat):
            return array([NaN]*size_)
        else:
            return ma.masked_invalid([getattr(individual, stat) for individual
                                      in pop.individuals]).mean(axis=axis_)

    # data.__class__ is Population
    if data.__class__ == Population:
        if stat == 'density':
            return [evaluate_genotype_density(data)]
        if stat == 'diversity':
            return evaluate_genotype_diversity(data, normalize),
            evaluate_phenotype_diversity(data, normalize)
        if stat == 'conservation':
            return get_average_matrix(data, bits)
        if stat == 'a.matrix':
            return get_a_matrix(data, bits)
        if 'clusters' in stat:
            return list(
                itertools.chain.from_iterable(
                    get_ind_stat(data, stat, pop_size)))
        if 'matrix'   in stat:
            return get_ind_stat(data, stat, pop_size).mean(axis=0)
        # redundant:
        else:
            return get_ind_stat(data, stat, pop_size)

    # data.__class__ is numpy.ndarray
    if type(data) == ndarray:
        if stat == 'survivability':
            if any(data):
                return data.mean(axis=0) #mean() because of MemoryError
            else:
                return [NaN]*bits**2#*pop_size

    # data.__class__ is EnsembleAnalysis
    if data.__class__ == EnsembleAnalysis:
        if stat == 'robustness':
            if data and hasattr(data, stat):
                return getattr(data, stat)
            else:
                return [NaN]*pop_size

    return get_ind_stat(data, stat, pop_size)

# reshape data
def get_stat_shape(
        stat, runs, datapoints, pop_size=P, bits=N, ps=1, qs=1, experiment=''):
    if 'clusters' in stat:
        return -1
    if 'stable.pop' in experiment:
        if 'k' in stat:
            return (datapoints, runs*pop_size, bits)
        else:
            return (datapoints, runs*pop_size)
    if 'qp.pop' in experiment and datapoints == 1:
        if 'matrix' in stat:
            return (ps, qs, -1, bits, bits)
        return (ps, qs, runs*pop_size)
    if stat == 'density':
        return (runs, datapoints)
    if stat == 'diversity':
        return (runs, datapoints, -1)
    if 'conservation' in stat:
        return (runs, datapoints, bits, bits)
    if 'matrix' in stat:
        return (runs, datapoints, bits, bits)#, pop_size
    if (stat == 'survivability' or stat == 'robustness') and 'qp' in experiment:
        return (runs, datapoints, bits, bits)#, pop_size
    if stat == 'qp':
        return (runs, datapoints, pop_size, 2)
    else:
        return (runs, datapoints, pop_size)

def get_stat(
        shortname, stat='phenotype', pop_stat_func=get_ind_stat, save_=True,
        **args):

    if type(pop_stat_func) is str:
        pop_stat_func = globals()[pop_stat_func]
    if type(stat) is list:
        stat = stat[0]

    generations = get_params(shortname, stat)[2]
    generations = get_pop_datapoints(shortname, generations)[1:]
    runs = get_runs(shortname)
    data = array([[pop_stat_func(get_pop(shortname, run, g)[0], stat)
                   for g in generations]
                  for run in runs])
    # save data
    if save_:
        filename = atleast_1d(pop_filenames[shortname])[-1]
        filename += print_run(run) + print_generations(g)
        save(data_dir + filename.replace('pop.evolution', stat), data)

    return data

# getattr(pop, stat)
def non_evolved_autoregulation(
        stats_, stable = False, binary = True, n_runs = 10, pop_size = int(1e4),
        p = None, q = None, period = 1, ps = ps_all,
        genotype_func = 'diagonal_p_genotype',
        experiment = 'autoregulation', extension = '.dat', verbose = True):

    pops = ((get_non_evolved_pop(
        experiment, stable, run, x, p, q, period,
        pop_size, binary, genotype_func, extension,
        verbose)[0]
             for run in range(1, n_runs+1))
            for x in ps)

    data = ma.masked_invalid([[[get_pop_stat(pop, stat, pop_size, experiment)
                                for stat in stats_]
                               for pop in runs]
                              for runs in pops])

    return data.swapaxes(-2, -1).reshape(ps.size, -1, len(stats_))

def get_qp_pop_robustness_and_survivability(
        stats_, n_runs=10, period=1,
        act_fraction=p, p='ps_all', q='qs_all',
        pop_size=1000, pop_stat_func='get_ind_stat',
        stable=True, experiment='qp.pop.stats',
        genotype_func='qp_genotype',
        reshape_=False, **args):

    # population's mean robustness.matrix to save space and memory
    if 'matrix' in stats_[0]:
        f = lambda x: x.mean(axis=0)
    else:
        f = None

    act_fractions = [None] if 'qp.pop' in experiment else [act_fraction]
    ps = ps_all if p == 'ps_all' else [p]
    qs = qs_all if p == 'qs_all' else [q]

    #ps = [p] if ('clusters' in stats_[0] or not p or p >= .5) else (ps_all if p == 'ps_all' else arange(0, 1 + p, p))
    #qs = [q] if ('clusters' in stats_[0] or not q or q == .5) else (qs_all if q == 'qs_all' else arange(0, 1 + q, q))

    return get_stable_pop_all_data(
        stats_, n_runs, act_fractions, ps, qs, period, pop_size, pop_stat_func,
        f, stable, True, experiment.rstrip('.stats'), genotype_func)

def get_stable_pop_all_data(
        stats_, n_runs=200, act_fractions=arange(0, 1.1, .1),
        ps=[None], qs=[None], period=1, pop_size=1000,
        pop_stat_func='get_pop_stat', function = None, stable=False,
        binary=True, experiment='stable.pop',
        genotype_func='diagonal_p_genotype', bits=N, extension='.dat',
        reshape_=False, verbose=True, save_=True):

    # init
    data = {}
    for stat in stats_:
        data[stat] = []
    for act_fraction, p, q in itertools.product(act_fractions, ps, qs):

        for run in range(1, n_runs + 1):
            pop = get_non_evolved_pop(
                experiment, stable, run, act_fraction, p, q, period, pop_size,
                binary, genotype_func, extension, verbose)[0]

            for stat in stats_:
                _data = get_stat_function(
                    stat, pop, pop_size, bits, experiment, pop_stat_func, p = p)
                if any(function):
                    _data = function(_data)
                if pop and 'clusters' in stat:
                    _data = list(itertools.chain.from_iterable(_data))
                if type(_data) == float:
                    data[stat].append(_data)
                else:
                    data[stat].extend(_data)
    if save_:
        act_fraction = act_fractions[0] if len(act_fractions) == 1 else None
        p = ps[0] if len(ps) == 1 else ps[1]-ps[0] # p_step
        q = qs[0] if len(qs) == 1 else qs[1]-qs[0] # q_step
        suffix = get_suffix(
            stable = stable, binary = binary,
            act_fraction = act_fraction, p = p, q = q,
            period = period, genotype_func = genotype_func,
            pop_size = pop_size, run = n_runs, experiment = experiment)

        for stat in stats_:
            if reshape_:
                new_shape = get_stat_shape(
                    stat, n_runs, len(act_fractions), pop_size, bits,
                    len(ps), len(qs), experiment)
                data[stat] = reshape(data[stat], new_shape)

            filename = '-'.join((stat, experiment, suffix))
            save(data_dir + filename + '.npy', data[stat])

    return data

# stat = 'robustness' or 'survivability'
def get_stable_pop_robustness_or_survivability(
        stat, n_runs=200, pop_size=1000, genotype_func= 'diagonal_p_genotype',
        experiment='stable.pop', extension='.dat', verbose=True, save_=True):

    boxplot_data = []
    heatmap_data = []

    for p in arange(0, 1.1, .1):
        matrix = []
        for run in range(1, n_runs + 1):
            pop = get_non_evolved_pop(
                experiment, p, run, pop_size, genotype_func, extension,
                verbose)[0]
            boxplot_data.append(([p]*pop_size, getattr(pop, stat)))
            matrix.append(
                mean([getattr(ind, stat) for ind in pop.individuals], axis = 0))
        heatmap_data.append(mean(matrix, axis = 0))
    return boxplot_data
    heatmap_data = reshape(heatmap_data, (-1, ind.size, ind.size))
    boxplot_data = swapaxes(boxplot_data, 0, 1).reshape(2, -1)
    if save_:
        fname = ('p_vs_%s-stable.pop-diagonalp_P%s-n%d.npy' %
                 (stat, print_int_or_float(pop_size), n_runs))
        save(data_dir + fname, boxplot_data)
        fname = ('%s-matrix-stable.pop-diagonalp_P%s-n%d.npy' %
                 (stat, print_int_or_float(pop_size), n_runs))
        save(data_dir + fname, heatmap_data)
    return boxplot_data, heatmap_data

# data=load(data_dir + 'p_vs_robustness-stable.pop-diagonalp_P1000-n200.npy') or
#     =load(data_dir + 'p_vs_survivability-stable.pop-diagonalp_P1000-n200.npy')
def get_percentiles(data, bins=None, step=.1, discrete=False):

    if not any(bins):
        if step == .1 and discrete:
            bins = ps_all
        else:
            bins = arange(0, 1.1, step)

    # data[0] = p
    if not discrete:
        bins[-1] += 1e-10
        return array([data[0][where((data[1]>=a) & (data[1]<b))[0]]
                      for a, b, in zip(bins, bins[1:])])

    # data[1] = robustness or survivability
    else:
        return array([data[1][where(data[0] == p)[0]]
                      for p in bins])

def bar_stack_states(data, plotting_=False):

    data = array([[get_vector(decimal, def_basename, def_base)
                   for decimal in decimals]
                  for decimals in data])

    a = array([[get_number_of_positive(v)
                for v in vectors]
               for vectors in data])

    n = array([bincount(x) for x in a]).T

    if plotting_:
        return bar_stacked(n)
    else:
        return data, a, n

# data = load(data_dir + 'robustness.matrix-qp.pop-stable_P1000-n10.npy') or
#      = load(data_dir + 'survivability.matrix-qp.pop-stable_P1000-n10.npy')
# shape(data) = (11, 11, 10000, 10, 10) = (p, q, runs*pop_size, bits, bits)
# ma.masked_invalid(data) ?
def difference_diag_nondiag_robustness_or_survivability(data, plotting_=True):
    if not plotting_:
        return array([[ma.mean([difference_diag_nondiag(x, ma.mean)
                                for x in q])
                       for q in p]
                      for p in data])

#def plot_robustness_and_survivability_heatmaps():


##################
## QPG, QPS, QPR #
##################

## Get data

# in_data  shape: (-1, 2)
# out_data shape: (11, 91)
def compute_density_matrix(data):
    qp = zeros((ps_all.size, qs_all.size))
    for x in data:
        if not all(isnan(x)):
            qp[x[0], x[1]] += 1
    return qp

# ATTENTION: run3.py passes myrun.run as g and myrun.samples as samples
def get_qp_pop_matrix(
        shortname, g, samples=None, stat='qp', save_=True, **args):

    data = compute_density_matrix(reshape([get_ind_stat(get_pop(shortname, run, g)[0], stat) for run in get_runs(shortname, samples)], (-1, 2)))
    if save_:
        filename = '-'.join(('qpg.matrix', shortname)) + print_run(run) + print_generations(g)
        save(get_experiment_dir(stat, shortname) + filename, data)
    return data

#
def get_qpg_pop_matrix(
        shortname, n_runs=None, samples=None, save_=True, **args):
    return array([get_qp_pop_matrix(shortname, g, n_runs)
                  for g in get_generations(shortname, samples)])

#
'''def get_qpg_matrix(shortname, norm_ = [10,90], save_ = True):
    data = get_stats_vs_g(shortname, ['p', '1-p']) # shape = (2, 104, 150000)
    data = data.transpose(1,2,0) # shape = (104, 150000, 2)
    data = data*reshape(norm_, (1,1,2))
    #(datapoints.size, stats.size, runs.size*pop.size)
    data = array([compute_density_matrix(g) for g in data])
    if save_: save(data_dir + '.matrix-'.join(('qpg', shortname)), data)
    return data
'''
#
def get_qpg_matrix(shortname, save_=True, **args):
    in_data = get_pop(shortname, extension = '.npy', experiment = 'qp')[0]
    if any(in_data):
        out_data = [compute_density_matrix(g) for g in in_data]
        if save_: save(data_dir + '.matrix-'.join(('qpg', shortname)), out_data)
        return out_data
    return 0

def get_stat_vs_g(shortname, stat, swapaxes_=True, reshape_=True):

    # load
    data = get_pop(shortname, None, None, '', '.npy', True, stat, data_dir)[0]

    # swap generations (datapoints) to first (0) dimension
    if swapaxes_:
        data = data.swapaxes(0, 1)

    # reshape data
    if reshape_:
        #generations = pop_model(shortname)[2]
        #datapoints = get_pop_datapoints(shortname, generations)[1:]
        #n_datapoints = len(datapoints)
        n_datapoints = data.shape[0]
        new_shape = (n_datapoints, -1)
        data = data.reshape(new_shape)

    # mask NaNs
    data = ma.masked_invalid(data)

    return data

# returns shape = (3, 104, 150000)
def get_stats_vs_g(shortname, stats_ = ['p', '1-p', 'robustness']):
    n_datapoints = len(get_pop_datapoints(shortname, pop_model(shortname)[2])[1:])
    return ma.masked_invalid([
        swapaxes(
            get_pop(
                shortname, None, None, '', '.npy',
                True, stat, data_dir)[0], 0, 1).reshape(n_datapoints, -1)
        for stat in stats_])

#
def get_mean_stats_at_g(g, shortname='m3r05b', stats_=['p', 'q'], f=ma.mean):
    return f(get_stats_vs_g(shortname, stats_)[:,g], axis=1) / [
        get_norm(stat) for stat in stats_]

#
def get_qp_from_qp_matrix(shortname):
    p = []; q = []
    for g in get_generations(shortname):
        fname = '-'.join(('qpg.matrix', shortname))
        fname += print_run(pop_model(shortname)[9])
        fname += print_generations(g)
        data  = load(get_experiment_dir('qp', shortname) + fname + '.npy')
        p.append(dot(ps_all, data.sum(axis=1)) / data.sum())
        q.append(dot(qs_all, data.sum(axis=0)) / data.sum())
    return p,q

#
# P=1e3 * n=10 = samples=1e4 for each p
def plot_evolution_oversampling(
        q, g, data=None, min_n=min_N, samples=100, qs_stable=arange(0,1,.05),
        fig=None, shortname='m3r05b', loc=2, journal='plos', size=20,
        marker='o', multi_image=False, panel=0, save_=True):

    if not fig:
        fig = figure()
    colors_ = iter(colors[journal])
    axis_ = [.5, 11.5, -.02, 1.02]

    # stable networks
    fname = ('robustness-qp.pop-stable-p0.1-q%.2f_P1000-n10' %
             find_nearest(qs_stable, q))
    stable = ma.masked_invalid(load(data_dir + fname + '.npy')).reshape(11, -1)
    y, yerr = get_y_yerr(stable)

    color = colors_.next()
    label = 'stable'

    # line plot
    #plot(ps_all, y, color, label=label)

    # shaded region
    #fill_between(ps_all, y+yerr, y-yerr, facecolor=color, alpha=0.5)

    # scatter plot
    #pa = [[scatter(x, y, size, color, marker)
     #      for y in p[:samples]]
      #    for x, p in zip(ps_all, stable)]

    # box plot
    bp = boxplot(list(stable))
    # unicolor all boxplot objects
    [setp(x, color=color) for x in bp.itervalues()]
    # set one label in random line
    x[0].set_label(label)
    # set x for next plot
    x = range(1,12)

    del stable, bp#, pa


    # evolved networks
    if not any(data):
        data = get_stats_vs_g(shortname) # shape = (3, 104, 150000)
        data = swapaxes(data, 0, 1)
    i = find_nearest_i(pop_generations, g)
    robustness = [get_all_x_at_qpg(data[i], p, find_nearest(qs_all, q))
                  for p in ps_all]

    # mean
    # cannot broadcast robustness to masked array, but every x is masked array
    y = ma.masked_invalid([p.mean() for p in robustness])
    yerr = ma.masked_invalid([p.std() for p in robustness])
    y.mask = yerr.mask = [p.size < min_n for p in robustness]

    color = colors_.next()
    label = 'generation %s' %print_generations(g, True)[2:]

    # line plot
    plot(x, y, color, label=label)

    # shaded region
    #fill_between(ps_all, y+yerr, y-yerr, facecolor=color, alpha=0.5)

    # scatterplot: add random noise in x to separate marks
    pa = [[scatter(fixed_jitter(x), y, size, color, marker, alpha=0.5)
           for y in p[:samples] if y]
          for x, p in zip(x, robustness)]
    del pa

    set_plot_kwargs(
        range(1,12), 'p', 'robustness', 'q=%.2f'%q, axis_, ax=fig.gca(),
        x_ticklabels=ps_all, loc=loc,
        multi_image=multi_image, panel_label=panel_labels[journal][panel])

    if save_:
        fname = 'robustness_stable_vs_evolved-q%.2f'%q
        fname += print_generations(g, True)
        save_fig(fname)

    return fig, data, robustness, r

#
def get_all_x_at_qpg(data, p, q, mask=False):
    qp = where((data[0] == p) & (data[1] == q))[0]
    if mask:
        return qp
    else:
        return data[2][qp]

# x could be robustness and stability
def get_qpx_at_g(
        data, g, min_n=min_N, stats_=['p', '1-p', 'robustness'],
        shortname=None, save_=False):

    qs = qs_all if stats_[1] == '1-p' else arange(91)

    if not any(data):
        data = get_stats_vs_g(shortname, stats_) # shape = (3, 104, 150000)
    data = swapaxes(data, 0, 1)

    qpx = ma.masked_invalid(
        [[get_all_x_at_qpg(data[g], p, q).mean()
          for q in qs]
         for p in ps_all])

    qpx.mask = array(
        [[get_all_x_at_qpg(data[g], p, q, True).size < min_n
          for q in qs]
         for p in ps_all])

    #if save_:
    return qpx

#
def get_qpr_vs_g(
        shortname, stats_=['p', '1-p', 'robustness'], min_n=min_N,
        fname='qpr', save_=True):

    n_datapoints = len(
        get_pop_datapoints(
            shortname, pop_model(shortname)[2])[1:])
    in_data      = get_stats_vs_g(shortname, stats_) # shape = (3, 104, 150000)
    out_data     = ma.masked_invalid([
        get_qpx_at_g(in_data, i, min_n, stats_)
        for i in xrange(n_datapoints)])
    if save_:
        fname = '.matrix-'.join((fname, shortname))
        save(data_dir + fname,               out_data.data)
        save(data_dir + fname + '.mask.npy', out_data.mask)
    return out_data

#
def get_qps_vs_g(
        shortname, stats_=['p', '1-p', 's'], min_n=min_N, fname='qps',
        save_=True):

    return get_qpr_vs_g(shortname, stats_, min_n, fname, save_)

#
def get_qpgsr_matrices(
        shortname, fnames=['qpg', 'qps', 'qpr'], q='1-p', save_=True):

    f = dict(
        zip(
            ['qpg', 'qps', 'qpr'],
            [get_qpg_matrix, get_qps_vs_g, get_qpr_vs_g]))

    s = dict(
        zip(
            ['qpg', 'qps', 'qpr'],
            [None, ['p', q, 's'], ['p', q, 'robustness']]))

    data = [f[fname](shortname, stats_=s[fname], save_=save_)
            for fname in fnames]

    return data


## Plot data

# ffmpeg -qscale 5 -r 5 -b 9600 -i qp*.matrix-m3r05b-%04d.png qp*.mp4
# ffmpeg -qscale 10 -r 6 -b 9600 -i qpgsr.matrix-m3r05b-%04d.png -s 1024x768 qpgsr.avi
# ffmpeg -r 6 -i qpgsr.matrix-m3r05b-%04d.png -b 2000k -minrate 2000k -maxrate 2000k -bufsize 1835k -s 1280x700 qpgsr.mp4
def plot_qpgsr(
        shortname='m3r05b', start=0, end=95,
        n_rows=1, n_cols=2, experiment='qpgsr.matrix',
        screen='macbook', titles=('stability', 'robustness'),
        elev=30, azim=45,
        y_ticks=None, plot_points=True, contour_projection=True,
        cmap1=None, cmap2=None, alpha1=.8, alpha2=.8, alpha3=.5,
        extension='png', plotting_=True):

    journal = 'plos'; n_col = 2
    qpg, qps, qpr = load_qpgsr_matrices(shortname)
    X, Y = meshgrid(qs_all, ps_all)
    zlim = qpg.max()

    if not cmap1:
        cmap1 = cm.jet
    if not cmap2:
        cmap2 = cm.binary#coolwarm

    fig_size = screen_size[screen]
    if not fig_size:
        fig_width = fig_widths[journal][n_col]
        # (width x height) in inches
        fig_size  = (fig_width*n_cols, fig_width*n_rows)

    for i in range(start, end):
        print i
        Z = qpg[i]
        S = qps[i]
        R = qpr[i]
        fig = figure(figsize=fig_size)

        # 2x 4d plot_surface
        if len(experiment.split('.')[0]) == 5:
            for j, N in enumerate((S, R)):
                a = fig.add_subplot(n_rows, n_cols, j+1, projection='3d')
                m = a.plot_surface(
                    X, Y, Z, rstride=1, cstride=1, alpha=alpha1,
                    facecolors=cmap1(N),
                    linewidth=0, antialiased=False, shade=False)
                if plot_points:
                    #k = where(Z == Z.max())
                    #x, y = X[k], Y[k]
                    x = sum(
                        [[z*k/90. for k,z in enumerate(p)]
                         for p in Z])/Z.sum() # q.mean
                    y = sum(
                        [[z*k/10. for z in p]
                         for k,p in enumerate(Z)])/Z.sum() # p.mean
                    c = a.plot(
                        x,y, zs=[0], zdir='z', marker='o', color = 'k',
                        alpha=alpha2)
                    #c = a.plot(1,y, zs=[0], zdir='z', marker='o', color = 'k')
                    c = a.plot([x,x], [0,1], 'k', alpha=alpha3)
                    c = a.plot([0,1], [y,y], 'k', alpha=alpha3)
                if contour_projection:
                    #c = a.contourf(X, Y, Z, zdir='z', offset=0, cmap=cmap2)
                    #c = a.contourf(X, Y, Z, zdir='x', offset=0, cmap=cmap2)
                    c = a.contourf(
                        X, Y, Z, zdir='y', offset=0, cmap=cmap2, alpha=alpha2)
                title  = 't=%s' %print_generations(pop_generations[i])[2:]
                title += '\nmean ' + titles[j] + ' = %.2f'%N.mean()

                set_plot_kwargs(
                    y_ticks, 'q', 'p', title, y_ticks=y_ticks, ax=a,
                    x_ticklabels=y_ticks, y_ticklabels=y_ticks,
                    elev=elev, azim=azim, z_label='#', x_lim=(0,1),
                    y_lim=(0,1), z_lim=(0,zlim))

                del N

        # 4d plot_surface
        elif len(experiment.split('.')[0]) == 4:
            N = R if 'qpgr' in experiment else S
            a = fig.gca(projection='3d')
            m = a.plot_surface(
                X, Y, Z, rstride=1, cstride=1, facecolors=cm.jet(N),
                linewidth=0, antialiased=False, shade=False)
            c = a.view_init(elev, azim)
            a.set_zlim(0, zlim)
            a.set_xlabel('q')
            a.set_ylabel('p')
            a.set_zlabel('#')
            title(print_generations(pop_generations[i])[2:])
            del N

        # 3d imshow
        elif len(experiment.split('.')[0]) == 3:
            if 'qpg' in experiment:
                data = Z
                cmap = cm.hot
                vmax = zlim
            else:
                data = R if 'qpr' in experiment else S
                cmap = cm.jet
                vmax = 1
            m = imshow(
                data, cmap, Normalize(vmin=0, vmax=vmax), 'auto', 'nearest')
            c = colorbar()
            a = fig.gca()
            title  = 't=%s' %print_generations(pop_generations[i])[2:]
            title += '\nmean ' + titles[1] + ' = %.2f'%data.mean()
            set_plot_kwargs(
                arange(0,qs_all.size,9), 'q', 'p', title,
                y_ticks = arange(ps_all.size), ax = a,
                x_ticklabels = qs_all[arange(ps_all.size)*9],
                y_ticklabels = ps_all)
            del data

        #show()
        if end-start == 1:
            return fig, X, Y, Z, S, R, a, m, c
        fig_name = '-'.join((experiment, shortname))
        save_fig(
            fig_name + '-%04d'%(i),
            get_figs_save_dir(None, experiment, shortname), format=extension)
        close()
        del Z, S, R, fig, a, m, c
    del qpg, qpr
    return 1

# ffmpeg -qscale 5 -r 5 -b 9600 -i qpg*%04d.png movie.mp4
# mlab.surf or imshow
def plot_qp_matrix(
        shortname, end=None, start=0, extension='.png', stat='qp', save_=True,
        _3d=True):

    generations = list(get_generations(shortname))[start:end]
    _close = False if len(generations) == 1 else True

    for i, g in enumerate(generations):
        fname = '-'.join(('qpg.matrix', shortname)) + print_run(pop_model(shortname)[9])
        data  = load(get_experiment_dir(stat, shortname) + fname + print_generations(g) + '.npy')

        if _3d:
            fig = mlab.figure(size=(800, 700))
            s   = mlab.surf(data, extent=(1,100,10,10,1,100))
            if save_:
                mlab.savefig(get_figs_save_dir(None, stat, shortname) + fname + '-%04d'%(i+start) + extension)
            if _close:
                mlab.close()
        else:
            fig = figure()
            Im  = imshow(data, cm.hot, None, 'auto', 'nearest')
            colorbar()
            set_plot_kwargs(
                None, '# positive in non-diagonal', '# positive in diagonal',
                shortname)
            if save_:
                save_fig(fname)
    return fig, data

# boxplot
def plot_robustness_vs_p_boxplot(
        qpr, g, min_n=min_N, plotting_=True, save_=True):
    data = array([qpr[2,g][where(qpr[0,g] == p)[0]]
                  for p in ps_all])
    mask = reshape([where(qpr[0,g] == p)[0].size > min_n
                    for p in ps_all], (ps_all.size, ))
    if plotting_:
        fig = figure()
        boxplot(data[mask])
        title = ('robustness: m3r05b at generation %s (min_n=%d)' %
                 (print_generations(pop_generations[g]), min_n))
        set_plot_kwargs(
            arange(data[mask].size)+1, 'p', 'robustness', title, ax=fig.gca(),
            x_ticklabels = ps_all[mask])
        if save_:
            save_fig(
                'robustness_vs_p-m3r05b' + print_generations(
                    pop_generations[g]) + '-boxplots')
        return fig, data, mask
    return data, mask

# imshow
def plot_qpx_at_g(
        data, g, cmap=None, norm=None, aspect='auto', interpolation='nearest',
        save_=True):

    if not cmap:
        cmap = cm.jet
    if not norm:
        norm = Normalize(vmin=0, vmax=1)
    fig = figure(); imshow(data, cmap, norm, aspect, interpolation); colorbar()
    title = ('robustness: m3r05b at generation %s (min_n=%d)' %
             (print_generations(pop_generations[g])[2:], min_n))
    xticklabels = qs_all[::10]
    xticks = xticklabels*qs_all.size
    xticklabels = ['%.2f'%x for x in xticklabels]
    yticklabels = ps_all
    yticks = yticklabels*ps_all.size
    set_plot_kwargs(
        xticks, 'q', 'p', title, y_ticks=yticks, ax=fig.gca(),
        x_ticklabels=xticklabels, y_ticklabels=yticklabels)
    if save_:
        save_fig(
            ('robustness_vs_qp-m3r05b' + print_generations(pop_generations[g])
             + '-heatmap'))
    return fig

##
def pop_autoregulation(
        shortname, stats_=['p', 'equal'], min_n=20, g_samples=9,
        error=False, n_datapoints=100, stat_label=get_stat_label,
        plotting_=True, _set_plot_kwargs=True, save_=False,
        extension='.npy', _dir=data_dir, verbose=True):

    data = []
    for stat in stats_:

        # get parameters
        params = get_params(shortname, stat)
        generations = params[2]
        n_runs, precision, filter_, suffix = params[9:13]
        datapoints = get_pop_datapoints(
            shortname, generations, filter_, n_datapoints)[1:]

        # get data
        _data = get_pop(
            shortname, n_runs, generations, suffix, extension, verbose, stat,
            _dir)[0]
        _data = ma.masked_invalid(_data)
        data.append(_data)

    x_norm = get_norm(stats_[0], shortname)
    y_norm = get_norm(stats_[1], shortname)

    i = range(0, datapoints.size, g_samples)# + [-1]

    z = [array([(x/x_norm, data[1][:,j][where(data[0][:,j] == x)].mean()/y_norm, data[1][:,j][where(data[0][:,j] == x)].std()/y_norm)
                for x in unique(data[0][:,j]) if data[1][:,j][where(data[0][:,j] == x)].size > min_n]).T for j in i]
    n_samples = [[(x, data[1][:,j][where(data[0][:,j] == x)].size) for x in unique(data[0][:,j])] for j in i]

    if plotting_:
        fig = figure()
        g = datapoints[i]
        labels = (print_generations(x)[1:] for x in g)
        #[errorbar(x[0], x[1], x[2], color = c, label = labels.next()) for x, c in zip(z, colors[None])]
        [plot_or_errorbar(x[0], x[1], error, x[2], color, label = labels.next()) for x, color in zip(z, colors[None]) if any(x)]

        if _set_plot_kwargs: set_plot_kwargs(arange(0,1.1,.1), stat_label(stats_[0]), stat_label(stats_[1]), shortname)
        # save plot
        if save_:
            fig_name = '-'.join(('pop.autoregulation', shortname))
            fig_name = shortname.join(('_vs_'.join(map(stat_label, stats_[::-1])) + '-', '-min_n%s'%min_n))
            save_fig(fig_name)

        return fig, data, z, n_samples

    else: return    data, z, n_samples


##################
## POP EVOLUTION #
##################

def pop_evolution(
        shortname, stats_=pop_stats, start=1, end=None,
        samples=None, datapoints=None, n_datapoints=100,
        pop_size=P, normalize=None, axis_=None, x_scale='log',
        experiment='pop.evolution', pop_dir=pop_dir,
        save_dir=get_figs_save_dir(None),
        fig_width=fig_widths[None][2],
        verbose=False, plotting_=True, _set_plot_kwargs=True,
        save_=True, extension='.dat', reshape_=True, **args):

    # init
    data = {}
    for stat in stats_:
        data[stat] = []
    if normalize == 'pop.size':
        normalize = pop_size

    # get parameters
    params = get_params(shortname, experiment)
    generations = params[2]
    bits = params[6]
    n_runs, precision, filter_, suffix = params[9:13]
    if not end:
        end = n_runs
    runs = arange(start, end + 1)

    if not any(datapoints):
        datapoints = get_pop_datapoints(
            shortname, generations, filter_, n_datapoints)
    datapoints = datapoints[1:samples]

    if plotting_:
        fig = figure(figsize=(fig_width, fig_height))

    # load and evaluate
    for run in runs:
        for g in datapoints:
            if stats_ == ['survivability']:
                # ATTENTION: pop is actually numpy.ndarray
                pop = get_pop(
                    shortname, run, g, suffix, extension, verbose, stats_[0])[0]
            else:
                pop = get_pop(shortname, run, g, suffix, extension, verbose)[0]
            for stat in stats_:
                new_data = get_stat_function(
                    stat, pop, pop_size, bits, experiment, normalize=normalize)
                data[stat].extend(new_data)

    # reshape and save output
    if plotting_ or save_:
        filename = atleast_1d(pop_filenames[shortname])[-1]
        filename += print_run(runs.size) + print_generations(g, precision)
        filename += suffix + extension.split('.')[0]

        for stat in stats_:

            # reshape data
            if reshape_:
                new_shape = get_stat_shape(
                    stat, runs.size, datapoints.size, pop_size, bits)
                data[stat] = reshape(data[stat], new_shape)

            if stat == 'qp':
                data[stat] = swapaxes(data[stat], 0, 1).reshape(
                    datapoints.size, -1, 2)
            # save data
            if save_:
                fname = data_dir + filename.replace('pop.evolution', stat)
                save(fname, array(data[stat]))

            # plot
            if plotting_ and stat != 'diversity' and stat != 'density':
                # mask an array where invalid values occur (NaNs or infs)
                y = ma.masked_invalid(data[stat])
                if stat == 'conservation':
                    y = abs(reshape(y, shape(y)[:2] + (-1,)))
                norm = get_norm(stat, shortname)
                y = y.mean(axis=0).mean(axis=1) / norm
                plot(datapoints, y, label = stat)

        if plotting_:
            if _set_plot_kwargs:
                set_plot_kwargs(
                    None, 'generations', str(stats_), shortname, axis_, x_scale)
            # save plot
            if save_:
                fname  = str(stats_) if len(stats_) > 1 else stat
                save_fig('-'.join((fname, shortname)), save_dir)
            return fig, data

    return data

def all_models_all_stats(
        stats_, shortnames, diversity=pop_diversity, statistic='mean',
        errors=pop_errors, individual_runs=False, threshold=.67,
        vline='', vlines_=None, step=1, runs_axis = 0,
        stat_label=get_stat_label, shortname_label=get_model_label,
        n_datapoints=100, pop_size=P, x_ticks=None, y_ticks=None,
        markers=markers, linestyle=linestyles[0], _grid=False,
        vspan=None, y_label='', title_='', axis_=None, x_scale='log',
        journal=None, n_col=2, suffix='', experiment='autoregulation',
        fig_name=None, _dir=data_dir, plotting_=True, verbose=True,
        fig=None, save_=False, extension='.npy'):

    stats_ = atleast_1d(stats_)
    shortnames = atleast_1d(shortnames)

    if plotting_ and not fig:
        fig = figure(figsize=(fig_widths[journal][n_col], fig_height))

    if shortnames.size <= len(colors[journal]):
        colors_ = iter(colors[journal])
    else:
        colors_ = iter(colors[None])

    color = colors_.next()

    linestyles_ = iter(linestyles[1:])

    for shortname in shortnames:
        for stat, marker in zip(stats_, markers):
            if verbose:
                print shortname, stat

            # get parameters
            params = get_params(shortname, stat)
            generations = params[2]
            bits = params[6]
            n_runs, precision, filter_, suffix = params[9:13]
            x = get_pop_datapoints(
                shortname, generations, filter_, n_datapoints)[1:]

            if stat == 'a.conservation':
                x = x[step::step]
                suffix += '-step%d'%step

            # load data
            data = get_pop(
                shortname, n_runs, generations, suffix, extension, verbose,
                stat, _dir)[0]

            # mask an array where invalid values occur (NaNs or infs)
            data = ma.masked_invalid(data)

            #if not plotting_: return data
            if data.ndim == 1:
                new_shape = get_stat_shape(stat, n_runs, x.size, pop_size, bits)
                data = reshape(data, new_shape)

            # diversity's shape:    (runs.size, datapoints.size, -1)
            if stat == 'diversity':
                for d_stat, k in diversity.items():
                    norm = get_norm(d_stat, shortname)
                    # mean and std
                    y    = data.T[k].mean(axis=1) / norm
                    if stat in errors:
                        yerr = data.T[k].std (axis=1) / norm
                    else:
                        yerr = None
                    # correct for small samples
                    x, y, yerr = correct_sample_size(
                        x, y, yerr, data, runs_axis, threshold)

                    if not plotting_:
                        return data, x, y, yerr

                    # label
                    label = ''

                    if shortnames.size > 1:
                        label += shortname_label(shortname)
                    if stats_.size > 1 or len(diversity) > 1:
                        label += stat_label(d_stat)
                    label += suffix

                    # plot
                    plot_or_errorbar(
                        x, y, stat in errors, yerr, color, marker, linestyle,
                        label)

                    color = colors_.next()

            # conservation's shape: (runs.size, datapoints.size, bits, bits)
            # other stats' shape:   (runs.size, datapoints.size, pop_size))
            else:
                y, yerr, data = get_y_yerr(
                    data, stat, statistic, errors, shortname)

                # correct for small samples
                #xx, y, yerr = correct_sample_size(
                 #   x, y, yerr, data, runs_axis, threshold)
                i = correct_sample_size2(data, runs_axis, threshold)

                if not plotting_:
                    return data, x, y, yerr, i

                # label
                label = ''
                if shortnames.size > 1:
                    label += shortname_label(shortname)
                if stats_.size > 1:
                    label += stat_label(stat)
                #label                         += suffix

                # color and linestyle for random model
                if 'm2' in shortname:
                    color = 'k'
                    linestyle = linestyles_.next()

                # plot
                plot_or_errorbar(
                    x, y, stat in errors, yerr, color, marker, linestyle, label,
                    i)

                if stat == vline:
                    xline = x[where(y == y.max())[0]]
                    axvline(xline, color='k', linestyle='--',
                            label=('max(%s): g = %s' %
                                   (vline, print_generations(xline, True)[2:])))

                if individual_runs:
                    norm  = get_norm(stat)
                    [plot(x[:run.mean(axis=1).compressed().size],
                          run.mean(axis=1).compressed() / norm, color + '.')
                     for run in data]

                color = colors_.next()

    if not y_label:
        y_label = get_y_label(stats_, stat_label)

    if not title_ and not journal:
        pass#title_ = str(shortnames)

    if not axis_:
        axis_ = list(axis())
        #axis_ = (1, generations, 0, 1)
        if any(stats_ != ['l']) and errors.values():
            axis_ = (.9, 1.1*generations, -.2, 1.2)

    if any(vlines_):
        vlines(
            vlines_, axis_[2], axis_[3], color='k', linestyle=linestyles[-1])

    if any(vspan):
        p = [axvspan(
            z[0], z[1], facecolor='gray', alpha=0.5, linestyle='solid')
             for z in vspan]

    set_plot_kwargs(
        x_ticks, 'generations', y_label, title_, axis_, x_scale,
        y_ticks = y_ticks)

    grid(_grid)

    if save_:
        if not fig_name:
            fig_name = get_figname(
                stats_, shortnames, stat_label, shortname_label, diversity,
                errors, suffix + statistic)

        save_fig(
            fig_name, get_figs_save_dir(journal, experiment),
            fig_prefixes[journal], dpis[journal], formats[journal])

    del data
    return fig, x, y, yerr, i

def compare_models(
        stat, filenames, labels, x, color='', markers=markers,
        n_runs=50, generations=1e6, error = 'within', precision=True,
        suffix='', experiment='pop.evolution', y_label='',
        title_='', x_scale='log', save_dir=get_figs_save_dir(None),
        save=False):

    fig = figure()

    for filename, label, marker in zip(filenames, labels, markers):
        # load data
        filename = filename.replace(experiment, stat) + print_run(n_runs)
        filename += print_generations(generations, precision) + suffix
        data = load(data_dir + filename + '.npy')
        # mask an array where invalid values occur (NaNs or infs)
        data = ma.masked_invalid(data)
        y = data.mean(axis=0).mean(axis=1)
        if error:
            yerr = between_and_within_variation_pop(data, error)
            errorbar(x, y, yerr, fmt = marker, label = label)
            #errorbar(
             #   bits, y, z, fmt = color + symbol, label = label)#, mfc = color)
        else:
            plot(x, y, color + marker, label=label)

    set_plot_kwargs(None, 'generations', y_label, title_, x_scale=x_scale)

    return fig

def compare_errors(
        stats_, shortnames, errors=error_types, n_datapoints=100,
        label_function=str, threshold=.67, runs_axis=0,
        colors=colors[None], markers=markers, linestyles=linestyles,
        linewidth=linewidth, y_label='std', title_='', x_scale='log',
        fig_width=fig_widths[None][2], fig_height=fig_height*1.2,
        suffix='', experiment='pop.evolution',
        save_dir=get_figs_save_dir(None), save_=False):

    fig = figure(figsize=(fig_width, fig_height))
    colors = iter(colors)

    for stat, shortname in itertools.product(stats_, shortnames):
        color = colors.next()
        # params
        generations = pop_model(shortname)[2]
        n_runs, precision, filter_, suffix = pop_model(shortname)[9:13]
        x = get_pop_datapoints(
            shortname, generations, filter_, n_datapoints)[1:]
        # load data
        fname = atleast_1d(pop_filenames[shortname])[-1]
        fname = fname.replace(experiment, stat)
        fname += print_run(n_runs) + print_generations(generations, precision)
        fname += suffix
        data = load(data_dir + fname + '.npy')
        # mask an array where invalid values occur (NaNs or infs)
        data = ma.masked_invalid(data)
        for error, linestyle, marker in zip(errors, linestyles, markers):
            yerr = between_and_within_variation_pop(data, error)
            # correct for small samples
            x, y, yerr = correct_sample_size(
                x, yerr, yerr, data, runs_axis, threshold)
            # label
            label = ''
            if len(stats_) > 1:
                label += stat
            if len(shortnames) > 1:
                label += label_function(shortname)
            if len(errors) > 1:
                label += error
            label += suffix
            # plot
            plot(
                x, yerr, color + linestyle + marker, lw=linewidth, label=label)

    if not title_:
        title_  = str(shortnames)
    set_plot_kwargs(None, 'generations', y_label, title_, x_scale=x_scale)

    if save_:
        fname = '-'.join((y_label, str(errors), str(shortnames)))
        save_fig(fname, save_dir)

    return fig

def pop_selection(
        shortname, n_runs=50, samples=None, generations=1e6, start=1,
        n_datapoints=100, filter_=None, precision=False,
        experiment='-selection', pop_dir=pop_dir, extension='.dat'):

    datapoints = get_int_datapoints(
        -0.1, log10(generations), n_datapoints, filter_)[1:samples]

    # load and evaluate
    for run in range(start, n_runs + 1):
        for g in datapoints:
            pop, fname = get_pop(shortname, run, g)
            if pop:
                pop.selection()
                save_file(fname + experiment + extension, pop)

def plot_conservation_profile(
        shortnames, bits=N, stat='conservation', journal=None,
        colors=colors[None], save_=False):
    fig = figure()
    colors = iter(colors)
    norm = pop_normalization[stat] if stat in pop_normalization else 1
    for shortname in shortnames:
        # load data
        filename = atleast_1d(
            pop_filenames[shortname])[-1].replace('pop.evolution', stat)
        runs, precision, filter_, suffix = pop_model(shortname)[9:13]
        generations = pop_model(shortname)[2]
        filename += print_run(runs) + print_generations(
            generations, precision) + suffix
        data = load(data_dir + filename + '.npy')
        # mask an array where invalid values occur (NaNs or infs)
        data = ma.masked_invalid(data)
        # flatten genotypes:
        # conservation's shape = (runs.size, datapoints.size, bits, bits)
        data = reshape(data, shape(data)[:2] + (-1,))
        # we're interested in the normalized absolute value of:
        # -pop_size <= conservation <= +pop_size
        data = abs(data)/norm
        # average over runs and generations
        y = data.mean(axis=0).mean(axis=0)
        # x-axis is gene ##
        plot(y, color=colors.next(), label = shortname)
    x = range(bits**2)[0::bits+1]
    vlines(
        x, axis()[2], axis()[3],
        colors='k', linestyles='--', label='autoregulation')
    set_plot_kwargs(x, 'gene ##', stat)
    if save_:
        save_fig(
            '-'.join(('_'.join((stat, 'profile')), str(shortnames))),
            get_figs_save_dir(journal))
    return fig

def heatmap_conservation(
        shortname, abs=array, cmap=None, norm=None, interpolation='nearest',
        experiment='conservation', extension='.npy', journal='plos', n_col=1,
        verbose=True, save_=True):

    # makeshift init because of corn (pylab has problems there),
    # if sys.platform == 'darwin': from pylab import *
    if not cmap:
        cmap = cm.RdBu
    if not norm:
        norm = Normalize(vmin = -1, vmax = 1)

    fig = figure(figsize=(fig_widths[journal][n_col], fig_height))

    # get parameters
    params = get_params(shortname, experiment)
    generations = params[2]
    n_runs, precision, filter_, suffix = params[9:13]

    # load data
    data = get_pop(
        shortname, n_runs, generations, suffix, extension, verbose,
        experiment)[0]

    # mask an array where invalid values occur (NaNs or infs)
    data = ma.masked_invalid(data)

    # plot; conservation's shape = (runs.size, datapoints.size, bits, bits)
    x = abs(data.mean(axis=0).mean(axis=0)/pop_normalization[experiment])
    imshow(x, cmap, norm, interpolation = interpolation)
    colorbar()
    return fig, x

def multi_heatmap_conservation(
        shortnames, abs=array,
        cmap=None, aspect='auto', norm=None, interpolation='nearest', origin=None,
        stat='conservation', experiment='autoregulation',
        fig_name='average_matrix', journal=None, n_col=1,
        label_function=get_model_label_text,
        plotting_=False, verbose=True, save_=True, extension='.npy'):

    # makeshift init because of corn (pylab has problems there):
    # if sys.platform == 'darwin': from pylab import *
    if not cmap:
        cmap = cm.RdBu
    if not norm:
        norm = Normalize(vmin=-1, vmax=1)

    # even number of shortnames
    if len(shortnames) == 4:
        n_rows = n_cols = 2
    if len(shortnames) == 9:
        n_rows = n_cols = 3

    data = []
    for shortname in shortnames:

        # get parameters
        params = get_params(shortname, stat)
        generations = params[2]
        n_runs, precision, filter_, suffix = params[9:13]

        # load data
        _data = get_pop(
            shortname, n_runs, generations, suffix, extension, verbose, stat)[0]
        # mask an array where invalid values occur (NaNs or infs)
        _data = ma.masked_invalid(_data)
        # conservation's shape: (runs.size, datapoints.size, bits, bits)
        x = abs(_data.mean(axis=0).mean(axis=0)/pop_normalization[stat])
        data.append(x)

    fig = multi_image2(
        data, map(label_function, shortnames), n_rows, n_cols,
        cmap, norm, aspect, interpolation, fig_widths[journal][n_col])

    if save_:
        save_fig(
            fig_name, get_figs_save_dir(journal, experiment),
            fig_prefixes[journal], dpis[journal], formats[journal])

    return fig

def test_a_conservation(shortname):
    fig = figure()
    x = list(get_generations(shortname))[1:]
    for step in range(1,5):
        pop = get_pop(
            shortname, suffix='-step%d'%step, extension='.npy',
            experiment='a.conservation')[0]
        data = ma.masked_invalid(pop.swapaxes(0,1))
        y = [ma.median([ma.diag(run) for run in g]) -
             ma.median([flatnondiag(run) for run in g])
             for g in data]
        plot(x[step::step], y, label='step %d'%step)
    xscale('log'); leg = legend(loc=0)
    return fig

def get_a_conservation(
        shortname, step=1, conservation='a.conservation', stat='a.matrix',
        verbose=True, extension='.npy'):

    # get parameters
    params = get_params(shortname, stat)
    generations = params[2]
    n_runs = params[9]

    # load data
    data, fname = get_pop(
        shortname, n_runs, generations, extension=extension,
        verbose=verbose, experiment=stat)

    # a-matrix shape: (runs.size, datapoints.size, bits, bits)
    # difference between generations, with somewhat fancy slicing
    data = data[:,:-step:step,:] - data[:,step::step,:]

    # manuel
    if conservation:
        data = 1 - abs(data)
        fname = fname.replace(stat, conservation)

    save(fname + '-step%d' %step + extension, data)
    return data

def get_a_matrix(pop, bits=N):
    if pop and pop.size:
        genotypes = pop.get_genotypes(flat=True)
        # (pop.size, n_rows, n_cols) -> (n_rows, n_cols)
        return array([get_activating_fraction_binary(gene)
                      for gene in genotypes.T]).reshape((bits, bits))
    else:
        return reshape(array([NaN]*bits*bits), (bits, bits))

def get_average_matrix(pop, bits=N):
    if pop and pop.size:
        genotypes = pop.get_genotypes()
        # (pop.size, n_rows, n_cols) -> (n_rows, n_cols)
        return genotypes.mean(axis=0)
    else:
        return reshape(array([NaN]*bits*bits), (bits, bits))

def evaluate_genotype_density(pop):
    if pop and pop.size:
        genotypes = pop.get_genotypes()
        # (pop.size, n_rows, n_cols) -> ()
        return nonzero(genotypes)[0].size #numpy.count_nonzero()
    else:
        return NaN

def evaluate_genotype_diversity(pop, normalize=None):
    if pop and pop.size:
        genotypes = pop.get_genotypes(flat=True, _tuple=True)
        # unique
        n = len(set(genotypes))
        # entropy
        entropy = evaluate_entropy2(genotypes, normalize)
        # variance
        genotypes = reshape(genotypes, (pop.size, pop.bits, pop.bits))
        ## variance between genotypes (whole matrix)
        # minimized by all_equal and maximized by all different
        # (pop.size, n_rows, n_cols) -> (n_rows, n_cols) -> ()
        matrix  = genotypes.std(axis=0).mean()
        ## variance between rows
        # (pop.size, n_rows, n_cols) -> (pop.size, n_cols) -> ()
        rows    = genotypes.std(axis=1).mean()
        ## variance between columns
        # (pop.size, n_rows, n_cols) -> (pop.size, n_rows) -> ()
        columns = genotypes.std(axis=2).mean()
        return n, entropy, matrix, rows, columns
    else: return array([NaN]*5)

def evaluate_phenotype_diversity(pop, normalize = None):
    if pop and pop.size:
        phenotypes = [ind.phe_decimal
                      for ind in pop.individuals if ind.period == 1]
        # unique
        n = len(set(phenotypes))
        # entropy
        entropy = evaluate_entropy2(phenotypes, normalize)
        # variance
        phenotypes = get_vector_phenotypes(
            phenotypes, pop.bits, pop.basename, pop.base)
        # distance between phenotypes.
        # minimized by all_equal and maximized by all different.
        # (phenotypes.size, phenotype.size) -> (phenotype.size) -> ()
        distance = phenotypes.std(axis=0).mean()
        return n, entropy, distance
    else: return array([NaN]*3)

def pop_perturb_and_analyse(
        shortname, stats_ = pa_stats, n_runs=50,
        samples=None, generations=1e6, start=1,
        n_datapoints=100, filter_=None, precision=False,
        pop_size=P, normalize=n_perturb,
        colors=colors[None], axis_=None, x_scale='log',
        fig_width=fig_widths[None][2],
        suffix='', experiment='perturb_and_analyse',
        save_dir=get_figs_save_dir(None), plotting_=True,
        save_=True, extension='.dat'):

    # init
    data = {}
    for stat in stats_:
        data[stat] = []
    filename = pop_filenames[shortname].replace('pop.evolution', experiment)
    datapoints = get_int_datapoints(
        -0.1, log10(generations), n_datapoints, filter_)[1:samples]

    if plotting_:
        fig = figure(figsize = (fig_width, fig_height))

    # load and evaluate
    for run in range(start, n_runs + 1):
        _dir = get_pop_dir(shortname, run, pa_dir)

        for g in datapoints:
            fname = filename + print_run(run) + print_generations(g, precision)
            pa_data = load_file(_dir + fname + suffix + extension)

            for stat in data.keys():
                if pa_data:
                    data[stat].extend(pa_data[stat])
                else:
                    data[stat].extend([NaN]*pop_size)

    # reshape and save output
    filename += print_run(run) + print_generations(g, precision) + suffix
    for stat, color in zip(stats_, colors):
        new_shape = (n_runs, datapoints.size, pop_size)
        data[stat] = reshape(data[stat], new_shape)

        if save_:
            save(data_dir + filename.replace(experiment, stat), data[stat])

        if plotting_:
            norm = float(normalize) if 'entropy' not in stat else -log(normalize)
            # mask an array where invalid values occur (NaNs or infs)
            y = ma.masked_invalid(data[stat])
            plot(
                datapoints, y.mean(axis=0).mean(axis=1)/norm, c = color,
                label = stat)

    if plotting_:
        if not axis_:
            axis_ = (1, generations, 0, 1)
        set_plot_kwargs(None, 'generations', '', filename, axis_, x_scale)

        if save_:
            save_fig('-'.join((experiment, shortname)), save_dir)

    if plotting_:
        return fig, data

    return data

def pop_epistatis(
        shortname, n_runs=50, samples=None, generations=1e6, start=1,
        n_datapoints=100, filter_=None, precision=False, verbose=False,
        axis_=None, x_scale='log', fig_width=fig_widths[None][2],
        errors=pop_errors, suffix='', experiment='perturb_and_analyse',
        save_dir=get_figs_save_dir(None), plotting_=True, save_=True,
        extension='.dat'):

    # init
    y = []
    yerr = []
    stat = 'epistasis'
    filename = pop_filenames[shortname].replace('pop.evolution', experiment)
    datapoints = get_int_datapoints(
        -0.1, log10(generations), n_datapoints, filter_)[1:samples]

    if plotting_:
        fig = figure(figsize = (fig_width, fig_height))

    # load and evaluate
    for run in range(start, n_runs + 1):
        _dir = get_pop_dir(shortname, run, pa_dir)
        for g in datapoints:
            pa = load_file(
                _dir + filename + print_run(run) +
                print_generations(g, precision) + suffix + extension, verbose)
            y.   append(pa[stat].mean)
            yerr.append(pa[stat].std)

    # reshape and save output
    filename += print_run(run) + print_generations(g, precision) + suffix
    y    = reshape(y,    (n_runs, datapoints.size)).mean(axis=0)
    yerr = reshape(yerr, (n_runs, datapoints.size)).mean(axis=0)

    if save_:
        save(data_dir + filename.replace(experiment, stat), (y, yerr))

    if plotting_:
        plot_or_errorbar(datapoints, y, stat in errors, yerr, label=shortname)

    if plotting_:
        #if not axis_: axis_ = (1, generations, -1, 1)
        set_plot_kwargs(None, 'generations', '', filename, axis_, x_scale)
        if save_:
            save_fig('-'.join((stat, shortname)), save_dir)

    if plotting_:
        return fig, y, yerr
    return y, yerr

def pop_epistasis2(
        shortname, n_runs=50, n_mutations=16, generations=1e6,
        samples=None, start=1, n_datapoints=100, filter_=None,
        precision=True, pop_size=P, errors=pop_errors,
        colors=colors[None], axis_=None, x_scale='log',
        fig_width=fig_widths[None][2], suffix='', experiment='epistasis',
        save_dir=get_figs_save_dir(None), plotting_=True,
        verbose=True, save_=True, extension='.dat'):

    # init
    epistasis = []
    n = misc.comb(n_mutations, 2, 1)
    filename = atleast_1d(
        pop_filenames[shortname])[-1].replace('pop.evolution', experiment)
    x = get_pop_datapoints(shortname, generations, filter_, n_datapoints)[1:]
    if plotting_:
        fig = figure(figsize = (fig_width, fig_height))

    # load and evaluate
    for run in range(start, n_runs + 1):
        _dir = get_pop_dir(shortname, run, epi_dir)
        for g in x:
            data = load_file(
                _dir + filename + print_run(run) +
                print_generations(g, precision) + suffix + extension, verbose)
            if any(data):
                data = reshape(data, (pop_size, -1))
                data = ma.masked_invalid(data).mean(axis=1)
                epistasis.extend(data)
            else:
                epistasis.extend([NaN]*pop_size)

    # reshape and save output
    n_runs = n_runs + 1 - start
    filename += print_run(n_runs) + print_generations(g, precision)
    #return epistasis
    data = reshape(epistasis, (n_runs, x.size, pop_size))
    if save_: save(data_dir + filename, data)
    if plotting_:
        # mask an array where invalid values occur (NaNs or infs)
        data = ma.masked_invalid(data)
        y    = data.mean(axis=0).mean(axis=1)
        if experiment in errors.keys():
            yerr = between_and_within_variation_pop(data, errors[experiment])
        else:
            yerr = None
        plot_or_errorbar(
            x, y, experiment in errors.keys(), yerr, label=shortname)
        set_plot_kwargs(
            None, 'generations', experiment, filename, axis_, x_scale)
        if save_: save_fig('-'.join((experiment, shortname)), save_dir)
        return fig, data
    return data

def pop_robustness(
        shortname, stats_=['equal', 'intersect'],
        errors={'robustness': 'within'},
        n_datapoints=100, pop_size=P, normalize=n_perturb, samples=None,
        start=1, axis_=None, x_scale='log', fig_width=fig_widths[None][2],
        experiment='robustness', save_dir=get_figs_save_dir(None),
        plotting_=True, verbose=True, save_=True, extension='.dat',
        **args):

    # init
    data = {}
    for stat in stats_: data[stat] = []
    if plotting_: fig = figure(figsize = (fig_width, fig_height))

    # get parameters
    params = get_params(shortname, experiment)
    generations = params[2]
    n_runs, precision, filter_, suffix = params[9:13]
    runs = range(start, n_runs + 1)
    x = datapoints = get_pop_datapoints(
        shortname, generations, filter_, n_datapoints)[1:samples]

    # load and evaluate
    for run in runs:
        for g in datapoints:
            analysis = get_pop(
                shortname, run, g, suffix, extension, verbose, experiment)[0]
            for stat in data.keys():
                data[stat].extend(get_stat_function(stat, analysis, pop_size))

    # reshape and save output
    filename = atleast_1d(pop_filenames[shortname])[-1] + print_run(run)
    filename += print_generations(g, precision) + suffix
    for stat in data.keys():
        new_shape = get_stat_shape(stat, runs.size, datapoints.size, pop_size)
        data[stat] = reshape(data[stat], new_shape)

        if save_:
            save(data_dir + filename.replace('pop.evolution', stat), data[stat])

        if plotting_:
            # mask an array where invalid values occur (NaNs or infs)
            data[stat] = ma.masked_invalid(data[stat])
            # divide by number of perturbations (mutants)
            norm = float(normalize) if 'entropy' not in stat else -log(normalize)

            y = data[stat].mean(axis=0).mean(axis=1) / norm

            if experiment in errors.keys():
                yerr = between_and_within_variation_pop(
                    data[stat], errors[experiment]) / norm
            else:
                yerr = None

            plot_or_errorbar(
                x, y, experiment in errors.keys(), yerr, label = stat)

    if plotting_:
        if not axis_:
            axis_ = (1, generations, 0, 1)
        set_plot_kwargs(None, 'generations', '', filename, axis_, x_scale)

        if save_:
            save_fig('-'.join((experiment, shortname)), save_dir)

    if plotting_:
        return fig, data
    return data

def plot_stat_vs_parameter(
        shortnames, experiment='epistasis', xx='degree',
        errors={'epistasis': 'mean', 'equal': 'mean',
                'p': 'within', '1-p': 'within'},
        ge=1e5, le=1e6, threshold=.67, n_datapoints=100,
        runs_axis=0, axis_=None, x_scale='log', color='',
        save_=True, _dir=data_dir, extension='.npy', verbose=True,
        colors=colors[None], marker=markers[0], linestyle=linestyles[0],
        plotting_=False, plotting_2=True, journal=None, n_col=0,
        fig_width=fig_widths[None][2], fig_name=''):

    if plotting_:
        fig = figure(figsize=(fig_widths[None][2], fig_height))
        colors = iter(colors)

    norm = get_norm(experiment)

    stat = []
    error = []
    for shortname in shortnames:

        # get parameters
        params = get_params(shortname, experiment)
        generations = params[2]
        n_runs, precision, filter_, suffix = params[9:13]

        # load data
        data = get_pop(
            shortname, n_runs, generations, suffix, extension, verbose,
            experiment, _dir)[0]

        # mask an array where invalid values occur (NaNs or infs)
        data = ma.masked_invalid(data)

        if 'diff' in experiment:
            data = abs(data)

        y = data.mean(axis=0).mean(axis=1) / norm

        if experiment in errors.keys():
            yerr = between_and_within_variation_pop(
                data, errors[experiment]) / norm
        else:
            yerr = None

        x = get_pop_datapoints(
            shortname, generations, filter_, n_datapoints)[1:]

        # correct for small samples
        x, y, yerr = correct_sample_size(x, y, yerr, data, runs_axis, threshold)

        # plot
        if plotting_:
            plot(
                x, y, color=colors.next(), linestyle=linestyle, label=shortname)

        # find_nearest_i(pop_generations, 1e5) = 73
        # find_nearest_i(pop_generations, 1e6) = 89
        stat.append(y[(ge <= x) & (x <= le)].mean())

        if experiment in errors.keys():
            error.append(yerr[(ge <= x) & (x <= le)].mean())
        else:
            error.append([]*len(stat))

    # save experiment.mean() data
    if save_:
        data = array((shortnames, stat, error))
        fname = '-'.join((experiment, str(shortnames)))
        if experiment in errors.keys():
            suffix = '-error%s' %[errors[experiment]] + '-mean'
        else:
            suffix = ''
        suffix += print_generations(ge, False) + print_generations(le, False)
        savetxt(data_dir + fname + suffix + '.tsv', data.T, '%s\t%s\t%s')

    # set and save previous figure
    if plotting_:
        set_plot_kwargs(
            None, 'generations', experiment, str(shortnames), axis_, x_scale)

        if save_:
            save_fig(fname)

    # new fig
    if plotting_2:

        fig = figure(figsize=(fig_widths[journal][n_col], fig_height))
        title_ = 'Average between generations %.1e and %.1e' %(ge, le)

        i_x = pop_model.__doc__.split(', ').index(xx)
        x = [pop_model(shortname)[i_x] for shortname in shortnames]

        if xx == 'period' and 'm4r05p0' in shortnames:
            x[shortnames.index('m4r05p0')] = 2.3

        if xx == 'mut_bias':
            x = qp_biases#[::-1]

        plot_or_errorbar(
            x, stat, experiment in errors.keys(), error, color,
            label=get_stat_label(experiment))

        # set_plot_kwargs
        axis_ = list(axis())
        x_ticks = None
        y_label = get_stat_label(experiment)
        if   xx == 'degree':
            axis_[:2] = [ 1.5, 10.5]
            x_ticks   = x
        elif xx == 'rec_rate':
            axis_[:2] = [-.01,  .51]
        elif xx == 'period':
            axis_[:2] = [  .5,  6.5]
        set_plot_kwargs(x_ticks, xx, y_label, title_, axis_)
        if save_:
            if not fig_name:
                fig_name = '-'.join(
                    (experiment + '_vs_' + xx , str(shortnames)))
            save_fig(fig_name)

        return fig, x, stat, error

    return stat, error

def plot_stats_vs_parameter(
        stats_=['p', 'equal'], shortnames=all_shortnames[:-1],
        parameter='period', error='within', ge=1e5, le=1e6,
        threshold=.49, journal=None, n_col=2, save_=True,
        axis_=None, x_ticks=None, experiment='', fig_name=''):

    if parameter == 'period':
        axis_ = (.5,7.5,.1,1)
        x_ticks = [1, 2, 2.3] + range(3, 8)
        journal = 'plos'
        experiment = 'autoregulation'
        if stats_ == ['p', 'equal']:
            fig_name = 'p_r_vs_period'

    colors_ = iter(colors[journal])

    # get fig and x with first stat
    fig, x, y, yerr = plot_stat_vs_parameter(
        shortnames, stats_[0], parameter, {stats_[0]:error},
        ge, le, threshold, color=colors_.next(), save_=False,
        plotting_=False, plotting_2=True, journal=journal, n_col=n_col)

    # plot remaining stats
    for stat in stats_[1:]:
        # get stat
        y, yerr = plot_stat_vs_parameter(
            shortnames, stat, parameter, {stat:error}, ge, le, threshold,
            save_=False, plotting_=False, plotting_2=False)

        # plot
        plot_or_errorbar(
            x, y, error, yerr, colors_.next(), label=get_stat_label(stat))

    set_plot_kwargs(
        x_ticks, get_stat_label(parameter), get_y_label(stats_), '', axis_,
        ax=fig.gca(), x_ticklabels=x_ticks)

    if save_:
        save_fig(
            fig_name, get_figs_save_dir(journal, experiment),
            fig_prefixes[journal], dpis[journal], formats[journal])

    return fig

def plot_mean_median_mode(
        stat='p', data=None, shortname='m3r05b', vlines_=False,
        labels=['mean', 'median', 'mode'], screen='macbook', save_=True):

    if vlines_:
        fig = figure(figsize=screen_size[screen])
    else:
        fig = figure()

    generations = get_params(shortname, stat)[2]
    generations = get_pop_datapoints(shortname, generations)[1:]

    if stat == 'p':
        x_all = ps_all
        norm = 10.
        y0=.5;y1=1
    else:
        x_all = qs_all
        norm = 90.
        y0=.5;y1=.7

    if not any(data):
        data = get_stat_vs_g(shortname, stat)

    data = [data.mean(axis=1),
            ma.median(data, axis=1),
            stats.mstats.mode(data, axis=1)[0]]
            #ma.masked_equal([stats.mstats.mode(x)[0] for x in data], 0)]

    for j, y in enumerate(data):
        if y.mean() > 1: y /=norm
        plot(generations, y, label=labels[j])
        i = y.argmax()
        x = generations[i]
        if vlines_:
            z = find_nearest(x_all, y[i])
            label = 'max(x=%s)=%.2f(%s)' %(x, z, rint(z*norm).astype(int))
            vlines(x, y0,  y1, linestyle=linestyles[j+1], label=label)
            z = find_nearest(x_all, y[-10:].mean())
            label = 'mean=%.2f(%s)' %(z, rint(z*norm).astype(int))
            hlines(z,  1, 1e7, linestyle=linestyles[j+1], label=label)

    title_ = '' #shortname
    axis_ = None#(1, 1e7, y0, y1)
    set_plot_kwargs(None, 'generations', stat, title_, axis_, 'log',
                    _grid=True)

    if save_:
        fig_name = get_figname(stat, shortname) + 'mean_median_mode'
        save_fig(fig_name)

    return fig, data

def get_data_corrected_interval(
        shortname, stat, swapaxes_=True, reshape_=True,
        n_datapoints=15, threshold=.67):

    # load
    data = get_stat_vs_g(
        shortname, stat, swapaxes_=swapaxes_, reshape_=reshape_)

    if swapaxes_:
        interval = correct_sample_size2(data, 1, threshold)[-n_datapoints:]
        return data[interval]
    else:
        interval = correct_sample_size2(data, 0, threshold)[-n_datapoints:]
        return data[:,interval]

def get_g_threshold(
        shortnames=neutbias_shortnames, stat='p',
        n_datapoints=15, threshold=.67):

    g = []
    for shortname in atleast_1d(shortnames):
        # load
        data = get_stat_vs_g(shortname, stat)

        # get interval
        i = correct_sample_size2(data, 1, threshold)[-n_datapoints:]
        ge = print_generations(pop_generations[i[0]], True)[2:]
        le = print_generations(pop_generations[i[-1]], True)[2:]
        g.append([shortname, ge, le])

    return g

def get_qp_expection(
        stat='p', error=-1, shortnames=neutbias_shortnames, n_datapoints=15,
        threshold=.67):

    Y = []
    Yerr = []
    CIs = []
    for shortname in shortnames:

        if type(error) is int:
            # load, correct sample size and get interval
            data = get_data_corrected_interval(
                shortname, stat,
                n_datapoints=n_datapoints, threshold=threshold)
            Y.append(data.mean())
            Yerr.append(between_and_within_variation(data)[error])

        elif type(error) is str:
            # do not swapaxes or reshape
            data = get_data_corrected_interval(
                shortname, stat, False, False,
                n_datapoints=n_datapoints, threshold=threshold)
            y, yerr = get_y_yerr(data, stat, errors={stat:error})[:2]
            Y.append(y.mean())
            Yerr.append(yerr.mean())

    return array(Y), array(Yerr)

def get_mannwhitneyu(shortnames, stat):
    x = get_data_corrected_interval(shortnames[0], stat).flatten()
    y = get_data_corrected_interval(shortnames[1], stat).flatten()

    return stats.mstats.mannwhitneyu(x, y)

def get_q_mannwhitneyu(stat='1-p'):
    return [get_mannwhitneyu([m3, m2], stat)
            for m3, m2 in zip(mutbias_shortnames, neutbias_shortnames)]

def plot_q_histogram(shortnames, stat='1-p', statistic=ma.mean, fig=None):
    #for m3, m2, bias in zip(mutbias_shortnames, neutbias_shortnames, qp_biases)

    if not fig:
        fig = figure()
    colors_ = colors['plos'][:2]
    labels = ['no target', 'random']

    m3 = get_data_corrected_interval(shortnames[0], stat).flatten()
    m2 = get_data_corrected_interval(shortnames[1], stat).flatten()

    hist(
        [m3, m2], append(qs_all, [1.01]), align='left', color=colors_,
        label=labels)

    p_value = stats.mstats.mannwhitneyu(m3, m2)[1]
    q_m3, q_m2 = statistic(m3), statistic(m2)
    title_ = ('%s = %.2f, %s = %.2f\nP ~ %.1g Mann-Whitney U' %
              (labels[0], q_m3, labels[1], q_m2, p_value))
    axis_ = list(axis())
    axis_[0] = min(m3.min(), m2.min()) - .005
    axis_[1] = max(m3.max(), m2.max()) + .005
    set_plot_kwargs(None, 'q', '', title_, axis_)

    return fig

def pop_densities(
        degrees, threshold=0, runs=n_runs, generations=1e6,
        n_datapoints=100, filter_=None, rec_rate='05', models=(4,3,2),
        runs_axis=0, colors = colors[None], axis_=None, x_scale='log',
        fig_width=fig_widths[None][2],
        suffix='', experiment='density', save_dir=get_figs_save_dir(None),
        save_=True, extension='.npy', verbose=False):

    data = []
    for degree in degrees:
        fig = figure(figsize = (fig_width, fig_height))
        shortnames = ['m%sr%sk%s' %(model, rec_rate, degree)
                      for model in models]
        for shortname in shortnames:
            if verbose: print shortname
            # get x
            x = get_pop_datapoints(
                shortname, generations, filter_, n_datapoints)[1:]
            # get data
            filename  = atleast_1d(
                pop_filenames[shortname])[-1].replace(
                    'pop.evolution', experiment)
            densities = load(
                data_dir + filename + print_run(runs)
                + print_generations(generations) + suffix + extension)
            densities = ma.masked_invalid(densities)
            # get mean and variance
            y    = densities.mean(axis=0)/100/500
            yerr = densities.std (axis=0)/100/500
            # correct for small samples
            x, y, yerr = correct_sample_size(
                x, y, yerr, data, runs_axis, threshold)
            # plot
            errorbar(x, y, yerr, label = shortname)
            data.append((x,y,yerr))

        set_plot_kwargs(None, 'generations', experiment, 'K = %s' %degree,
                        axis_, x_scale)
        if save_:
            filename = '-'.join((experiment, str(shortnames)))
            save_fig(filename, save_dir)

    if save_:
        filename = (experiment + '-' + 'models%s' %str(models) + '_degrees%s' %
                    str(degrees))
        data = reshape(data, (len(degrees), len(models), -1))
        save(data_dir + filename, data)

    return fig, data

def get_average_correlation(data, window=2):
    stats_ = ['p', 'q', 'stability', 'robustness']
    shortname = 'm3r05b'
    x = get_pop_datapoints(shortname, pop_model(shortname)[2])[1:]
    p = array([])
    q = array([])
    for run in data.T:#numpy.rollaxis(data, 2):
        df = pd.DataFrame(run, columns=stats_, index=x)
        correls = pd.rolling_corr_pairwise(df, 2)
        p = append(p, correls.ix[:, 'robustness', 'p'].values)
        q = append(q, correls.ix[:, 'robustness', 'q'].values)
    return (ma.masked_invalid(p).reshape(104, -1),
            ma.masked_invalid(q).reshape(104, -1))

def plot_step_functions(dim=4, save_=True):
    # params
    bits = 4; density = c; binary = True

    # plot options
    fig_name = 'ultrasensitivity'
    x_label = r'$\displaystyle\sum_{j=1}^N w_{ij}\left(s_j(t)+\epsilon\right)$'
    y_label =  r'$s_i (t+1)$'
    steps = ('uniform', 'quantiles')
    titles = ('ultrasensitive', 'less sensitive')
    color = 'r'
    linewidth = 2
    n_samples = 1e4
    axis_ = (-2.66, 2.66, -1.33, 1.33)

    # figure
    # size
    journal = 'plos'; n_col = 2
    n_rows = 1; n_cols = 2
    fig_width = fig_widths[journal][n_col]
    fig = figure(figsize=(fig_width*n_cols, fig_width*n_rows))
    # change matplotlib's rc to use tex
    # (from http://matplotlib.org/users/usetex.html)
    rc('text', usetex=True)
    rc('font', family='serif')

    # plot
    for i, (step, title_) in enumerate(zip(steps, titles)):
        ax = fig.add_subplot(n_rows, n_cols, i+1)
        base, basename, jumps = get_base(dim, bits, steps=step)
        gpmap = get_gpmap(base, basename, density, bits, binary, jumps)
        x = linspace(axis_[0], axis_[1], n_samples)
        ax.plot(x, gpmap(x), color, lw=linewidth)
        set_plot_kwargs(None, x_label, y_label, title_, axis_)
        subplots_adjust(bottom=0.14)

    if save_:
        save_fig(fig_name)

    show()
    rcdefaults()
    return fig

def plot_step_functions2(dim=4, n_samples=1e5, save_=True):
    # params
    bits = 4; density = c; binary = True

    # plot options
    fig_name = 'ultrasensitivity-4panel'
    x_label = r'$\displaystyle\sum_{j=1}^N w_{ij}\left(s_j(t)+\epsilon\right)$'
    y_label =  r'$s_i (t+1)$'
    steps = ('uniform', 'quantiles')
    titles = (('extremes are favored', 'ultrasensitive'),
              ('states equally probable', 'less sensitive'))
    color = 'r'
    linewidth = 2
    a = -2.66; b = -a
    axis_ = (a, b, a/2, b/2)

    # figure
    # size
    journal = 'plos'; n_col = 2
    n_rows = 2; n_cols = 2
    fig_width = fig_widths[journal][n_col]
    fig = figure(figsize=(fig_width*n_cols, fig_width*n_rows))
    # change matplotlib's rc to use tex
    # (from http://matplotlib.org/users/usetex.html)
    rc('text', usetex=True)
    rc('font', family='serif')

    # plot
    for i, (step, title_) in enumerate(zip(steps, titles)):
        base, basename, jumps = get_base(dim, bits, steps=step)
        gpmap = get_gpmap(base, basename, density, bits, binary, jumps)

        # histogram
        ax = fig.add_subplot(n_rows, n_cols, i+1)
        x = (b - a) * random_sample(n_samples) + a
        ax.hist(gpmap(x), color='k')#color)
        set_plot_kwargs(base, y_label, '', title_[0], ax=ax,
                        panel_label=panel_labels[journal][i])

        # step function
        ax = fig.add_subplot(n_rows, n_cols, i+3)
        x = linspace(a, b, n_samples)
        ax.plot(x, gpmap(x), color, lw=linewidth)
        set_plot_kwargs(jumps, x_label, y_label, title_[1], axis_, ax=ax,
                        panel_label=panel_labels[journal][i+2])
        subplots_adjust(bottom=0.14)

    if save_:
        save_fig(fig_name)

    show()
    rcdefaults()
    return fig

'''
genotypes = [[generate_stable_genotype(bits, density, def_basename, gpmap, def_base, genotype_function, p, binary) for i in range(int(samples))] for p in arange(0, 1.1, .1)]
stability_mut = [[[test_stability(mutate(genotype, 1, change_sign = True), initial, basename, gpmap, base, bits)[0] for i in range(n_perturb)] for genotype, initial, phenotype in p] for p in genotypes]
stability_ini = [[[test_stability(genotype, i, basename, gpmap, base, bits)[0] for i in get_nearest_neighbours(initial, 3)] for genotype, initial, phenotype in p] for p in genotypes]

#genotypes = [[(genotype, initial, general_develop(genotype, initial, basename, gpmap, base, bits)[2]) for genotype, initial in p] for p in genotypes]
equal_mut = [[[test_stability(mutate(genotype, 1, change_sign = True), initial, basename, gpmap, base, bits)[2] == phenotype for i in range(n_perturb)] for genotype, initial, phenotype in p] for p in genotypes]
equal_ini = [[[test_stability(genotype, i, basename, gpmap, base, bits)[2] == phenotype for i in get_nearest_neighbours(initial, 3)]   for genotype, initial, phenotype in p] for p in genotypes]
'''


###########
## GRAPHS #
###########

def load_graph_degree(bit, alphas, gammas, deltas_in, deltas_out):
    data = []
    params = []
    lopp = itertools.product(alphas, gammas, deltas_in, deltas_out)
    for alpha, gamma, delta_in, delta_out in loop:
        params.append((bit, alpha, gamma, delta_in, delta_out))
        try:
            fname = ('N%sa%.2fg%.2fi%.2fo%.2f' %
                     (bit, alpha, gamma, delta_in, delta_out))
            k = loadtxt(gr_dir + fname + '.txt')
        except: #IOError ?
            k = NaN
        data.append(k)
    return ma.masked_invalid(data), array(params)

def get_minima(data, params, n_zeros=5, target=2, plotting_=True):
    zero = (target - data)
    if plotting_:
        fig = figure()
        plot(zero, 'o')
        [plot(param, '--') for param in params.T[1:4]]
    output = []
    for i in range(n_zeros):
        index = where(abs(zero) == abs(zero).min())[0][0]
        output.append((index, zero[index], params[index]))
        zero[index] = inf
    return output

def plot_zero_vs_bit(
        bits, alphas, gammas, deltas_in, deltas_out, n_zeros=3,
        color='', marker=markers[0], plotting_=True):

    output = []
    for bit in bits:
        data, params = load_graph_degree(
            bit, alphas, gammas, deltas_in, deltas_out)
        zeros = get_minima(data, params, plotting_=False)
        for i in range(n_zeros):
            output.append((bit, zeros[i][1], zeros[i][2][2]))

    if plotting_:
        fig = figure()
        output = array(output).T
        xx = []; yy = []
        for i in range(n_zeros):
            x = output[0][i::n_zeros]; xx.extend(x)
            y = output[2][i::n_zeros]; yy.extend(y)
            plot(x, y, color + marker)

    return output, xx, yy

def get_minimun_index(data, target = 2):
    zero = (target - data)
    return (array([where(abs(d) == abs(d).min())[0][0] for d in zero]),
            abs(zero).min(axis=1))

def parameters_dictionary(
        params, bits, index, alpha, gamma, delta_in, delta_out):
    d = [{'delta_out': delta_out, 'delta_in': delta, 'alpha': alpha,
          'gamma': gamma} for delta in deltas[index]]
    params = dict(zip(bits, d))


#################
## MULTI IMAGE2 #
#################

#data = [spin_genotype() for i in range(25)]
def multi_image(
        data, titles, n_rows=5, n_cols=5,
        cmap=None, norm=None, aspect='auto', interpolation='nearest',
        origin=None,
        fig_width=10, fig_height=12, w=0.15, h=0.14, figtitle=''):

    # makeshift init because of corn (pylab has problems there),
    # if sys.platform == 'darwin': from pylab import *
    if not cmap:
        cmap = cm.gray

    fig = figure(figsize = (fig_width, fig_height))
    t = fig.text(0.5, 0.96, figtitle,
                 horizontalalignment='center')
    cax = fig.add_axes([0.2, 0.05, 0.6, 0.04])
    ax = []
    images = []
    for i in range(n_rows):
        for j in range(n_cols):
            k = i*n_cols+j
            if k < len(data):
                pos = [0.075 + j*1.2*w, 0.12 + i*1.2*h, w, h]
                a = fig.add_axes(pos)
                if i > 0:
                    a.set_xticklabels([])
                if j > 0:
                    a.set_yticklabels([])
                a.set_title(titles[k])
                img = a.imshow(data[k], cmap, norm, aspect, interpolation,
                               origin = origin)
                images.append(img)
                ax.append(a)

    # The colorbar is based on the master (first/last) image.
    fig.colorbar(images[-1], cax, orientation='horizontal')

    return fig, images, ax, cax

# Set the first image as the master, with all the others
# observing it for changes in cmap or norm.

class ImageFollower:
    'update image in response to changes in clim or cmap on another image'
    def __init__(self, follower):
        self.follower = follower
    def __call__(self, leader):
        self.follower.set_cmap(leader.get_cmap())
        self.follower.set_clim(leader.get_clim())

def update_norm_and_cmap(images, norm=None, cmap=None):
    # makeshift init because of corn (pylab has problems there),
    # if sys.platform == 'darwin': from pylab import *
    if not cmap:
        cmap = cm.gray
    for i, im in enumerate(images):
        im.set_norm(norm)
        im.set_cmap(cmap)
        if i > 0:
            images[0].callbacksSM.connect('changed', ImageFollower(im))

#fig, images, ax, cax = multi_image(data, ['']*25, cmap = cm.jet)

#update_norm_and_cmap(images)

# We need the following only if we want to run this interactively and
# modify the colormap:

#axes(ax[0])     # Return the current axes to the first one,
#sci(images[0])  # because the current image must be in current axes.

#show()

def multi_image2(
        data, titles, n_rows = 2, n_cols = 2,
        cmap=None, norm=None, aspect='auto', interpolation='nearest',
        fig_width=12, fig=None, i_ax=1, x_label='j', y_label='i',
        x_step=1, x_ticks=None, y_ticks=None,
        x_ticklabels=None, y_ticklabels=None, cbar_label='', i_panel=None,
        axis_=(0.2, 0.02, 0.6, 0.03), journal=None, multi_image=False):

    # makeshift init because of corn (pylab has problems there),
    # if sys.platform == 'darwin': from pylab import *
    if not cmap:
        cmap = cm.RdBu

    if not i_panel:
        i_panel = i_ax-1

    # figsize = (width, height)  in inches
    if not fig:
        fig = figure(figsize=(fig_width*n_cols/2, fig_width*n_rows/2))

    cax = fig.add_axes(axis_)
    ax = []
    images = []
    for i, x in enumerate(data):
        a = fig.add_subplot(n_rows, n_cols, i+i_ax)
        ax.append(a)
        if not any(x_ticks):
            x_ticks = arange(0, shape(data[i])[1], x_step)
        if not any(y_ticks):
            y_ticks = arange(0, shape(data[i])[0])
        if not any(x_ticklabels):
            x_ticklabels = x_ticks + 1
        if not any(y_ticklabels):
            y_ticklabels = x_ticklabels
        set_plot_kwargs(
            x_ticks, x_label, y_label, titles[i], y_ticks=y_ticks, ax=a,
            x_ticklabels=x_ticklabels, y_ticklabels=y_ticklabels,
            multi_image=multi_image,
            panel_label=panel_labels[journal][i+i_panel])
        img = ax[i].imshow(data[i], cmap, norm, aspect, interpolation)
        images.append(img)
    cbar = fig.colorbar(images[0], cax, orientation='horizontal')
    cbar.set_label(cbar_label)
    #cbar.ax.set_xlabel(cbar_label)
    #cbar.ax.set_yticklabels(['Low', 'High'])
    #cbar.set_alpha(1)
    #cbar.draw_all()
    cbar.solids.set_edgecolor("face")
    draw()
    #show()
    return fig


def multi_analysis(
        data=[], sel=s, cut=1, min=zero_fit, start=0, end=G,
        step=50, interpolation='bilinear', cmap=None, origin='lower',
        norm=None):

    # makeshift init because of corn (pylab has problems there),
    # if sys.platform == 'darwin': from pylab import *
    if not cmap:
        cmap = cm.gray
    if not norm:
        norm = Normalize(vmin=1.0, vmax=1.1)

    param = '_s%s_cut%s_min%s' %(sel, cut, min)
    param_title = ('s = %s, cut = %s, min_fitness = %s\n500x1024x10' %
                   (sel, cut, min))

    for i in arange(start, end, step):
        filename = analysis_dir + 'evolved%s_analysis-G%s.dat' %(param, i)
        file     = open(filename, 'rb')
        analysis = cPickle.load(file)
        file.close()

        fig = figure()
        im  = imshow(
            analysis.period.T, interpolation=interpolation,
            cmap=cmap, origin=origin, norm=norm)

        title('G = %s\n%s' %(i, param_title))
        xlabel('population member #')
        ylabel('initial phenotype (decimal representation)')
        show()

        data.append(analysis.period)
        del analysis

    filename = data_dir + 'data%s.dat' %(sel, cut, min)
    file     = open(filename, 'wb')
    cPickle.dump(data, file)
    file.close()
    return data


def anovan(data):
    a, b, n = shape(data)
    Y    = data.mean() # 1 Grand mean
    YA   =                 data.reshape(a, -1).mean(axis=1).reshape(a, 1, 1)
    YB   = swapaxes(data, 0, 1).reshape(b, -1).mean(axis=1).reshape(1, b, 1)
    YAB  =                                data.mean(axis=2).reshape(a, b, 1)
    SSA  = n*b*sum((YA   - Y)**2) # 2
    SSB  = n*a*sum((YB   - Y)**2) # 3
    SSAB = n  *sum((YAB  - YA - YB + Y)**2) # 4 Interaction SS
    SSw  =     sum((data - YAB)**2) # 5 Within subgroups (error SS)
    dfA  = a-1
    dfB  = b-1
    dfAB = (a-1)*(b-1)
    dfw  = a*b*(n-1)
    MSA  = SSA  / dfA
    MSB  = SSB  / dfB
    MSAB = SSAB / dfAB
    MSw  = SSw  / dfw
    return [[SSA, SSB, SSAB, SSw], [MSA, MSB, MSAB, MSw], [dfA, dfB, dfAB, dfw]]
