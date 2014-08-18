import src.plotting_tools2
reload(src.plotting_tools2)
from src.plotting_tools2 import *

# 1) Stability and autoregulation are correlated

# fig 1
# P=1e4 * n=10 = samples=1e5 for each p
def fig_stability(
        data=None, boxplot_data=None, normed=True, param=5, loc=9, _save=True):

    #'p-boxplot_and_histogram-2panel'#'p-boxplot_and_stability_vs_p'
    fig_name = 'stability'
    journal = 'plos'; n_col = 2
    x_ticklabels = ['unstable', 'stable']

    n_rows = 1; n_cols = 2
    fig_width = fig_widths[journal][n_col]
    fig       = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    x = ps_all
    bins = array(range(12))/10.
    a_axis = [.5, 2.5, -.05, 1.05]
    b_axis = [-.05, 1.05, 0, 8e4]

    #fname = ('autoregulation-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary'
     #        + '-p1.0-diagonalp_P1e+04-n10.dat')
    #data.shape = (11, 100000)
    if data is None or boxplot_data is None:
        data, boxplot_data = get_boxplot_data()

    if normed:
        norm_ = float(data.shape[1])
        y_label = 'stability'
        reg_label = 'stability ~ exp(%dp)' %param
    else:
        norm_ = 1
        y_label = 'number of stable matrices, $n_f$'
        reg_label = '$n_f \sim exp(p)$'

    x_label = 'p'

    # panel A: p
    ax = fig.add_subplot(n_rows, n_cols, 1)

    # boxplot
    boxplot(boxplot_data)

    set_plot_kwargs(
        [1,2], '', x_label, _axis=a_axis, journal=journal, y_ticks=x, ax=ax,
        x_ticklabels=x_ticklabels,
        multi_image=True, panel_label=panel_labels[journal][0])

    # panel B: stability
    ax = fig.add_subplot(n_rows, n_cols, 2)

    # _vs_p
    #y, yerr = get_y_yerr(data)
    #plot_or_errorbar(x, y, True, yerr, colors[journal][0])
    #set_plot_kwargs(
     #   x, 'p', get_stat_label(stat), _axis = [-.05, 1.05, -.01, 1.200001],
      #  journal = journal, multi_image = True,
       # panel_label = panel_labels[journal][1]))

    # histogram
    hist(boxplot_data[1], bins, align='left')

    y_ticklabels = ax.get_yticks()/_norm

    nf = array([where(stable==True)[0].size for stable in data])

    xx, yy, ya, yb = polyfit_regression(
        x, nf, x_step=.1, x_label=x_label, y_label='n_f', label=reg_label,
        color=colors[journal][1])

    set_plot_kwargs(
        x, x_label, y_label, _axis=b_axis, journal=journal, ax=ax,
        y_ticklabels=y_ticklabels, loc=loc,
        panel_label=panel_labels[journal][1])

    #del data, boxplot_data

    if _save:
        save_fig(
            fig_name, get_figs_save_dir(journal, 'autoregulation'),
            fig_prefixes[journal], dpis[journal], formats[journal])

    return fig, data, boxplot_data

# 2) Evolution of autoregulation

# fig 2
# n_runs = [pop_model(shortname)[9] for shortname in six_shortnames]
# = [300, 200, 200, 100, 100, 100]
def fig_evolution(shortnames=six_shortnames, _save=True):

    titles = [get_model_label_text('m3r05b'),
              get_model_label_text('m3r05p2'),
              get_model_label_text('m4r05p0').split(', ')[0]]

    return multiple_plot(
        shortnames, fig_name='p_evolution', titles=titles, _save=_save)


# 3) Autoregulatory inputs are more conserved over time

# fig 3
# generation = 48 is maximum p
def fig_viability(
        data=None, generations=[1e6, 48], positions=[1,2],
        journal='plos', n_col=2, _save=True):

    fig_name = 'viability'
    experiment = 'autoregulation'

    cmap = cm.RdBu
    _norm = Normalize(vmin=0, vmax=1)

    n_rows = 2; n_cols = 2
    fig_width = fig_widths[journal][n_col]
    fig       = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    fname = ('survivability-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary'
             + '-p0.5-model3_G1e+07_period1_u0.1-binary_r0.5-n300_G1.0e+07')

    x_ticklabels = ['diagonal', 'off-diagonal']
    _axis = [positions[0] - 0.5, positions[-1] + 0.5, 0, 1]

    if not any(data):
        data  = ma.masked_invalid(load(data_dir + fname + '.npy'))
        #data.shape = (300, 104, 10, 10)
    # [89, 19], data.mean(axis=1).argmax() = 19 (max p)
    i_x = [where(pop_generations == g)[0][0] for g in generations]

    # panels A and B: diagonal vs off-diagonal boxplots
    for k, i in enumerate(i_x):
        a = fig.add_subplot(n_rows, n_cols, k+1)
        p = ma.masked_invalid([    ma.diag(run[i]) for run in data])
        q = ma.masked_invalid([flatnondiag(run[i]) for run in data])
        boxplot([p, q], positions=positions)
        set_plot_kwargs(
            positions, '', fig_name,
            'generation %s' %print_generations(generations[k], False)[2:],
            _axis=_axis, journal=journal, y_ticks=ps_all, ax=a,
            x_ticklabels=x_ticklabels, multi_image=True,
            panel_label=panel_labels[journal][k])

    # panels C and D: average matrix heatmaps
    w = data.mean(axis=0)
    fig = multi_image2(
        w[i_x], ['']*2, n_rows, n_cols, cmap, _norm, fig=fig, i_ax=3,
        cbar_label=fig_name)

    if _save:
        save_fig(fig_name, get_figs_save_dir(journal, experiment),
                 fig_prefixes[journal], dpis[journal], formats[journal])
    return fig, data


# n_runs = [pop_model(shortname)[9] for shortname in six_shortnames]
# = [300, 100]
def conservation_boxplots_2panel(
        fig=None, generations=[1e6, 103309], positions=[1,2], verbose=False):

    shortnames = ['m3r05b', 'm2r05b']
    titles = ['no target', 'random']
    x_ticklabels = ['diagonal', 'off-diagonal']
    _axis = [positions[0] - 0.5, positions[-1] + 0.5, 0, 1]
    i,j = [where(pop_generations == g)[0][0] for g in generations]

    matrices = []; boxplots = []
    for k, shortname in enumerate(shortnames):

        # load data
        data = get_pop(
            shortname, extension='.npy', experiment='a.matrix')[0]
        data = data.swapaxes(0,1)
        data = ma.masked_invalid(data)

        # analysis
        # ma.count_masked(y) == 0 for generations = [1e6, 103309]
        conservation = 1-abs(data[i]-data[j])
        p = ma.masked_invalid([    ma.diag(run) for run in conservation])
        q = ma.masked_invalid([flatnondiag(run) for run in conservation])

        # save
        matrices.append(conservation.mean(axis=0))
        boxplots.append([p, q])

        # plotting
        if fig:
            ax = fig.add_subplot(n_rows, n_cols, k+1)
            boxplot([p, q], positions=positions)
            set_plot_kwargs(
                positions, '', fig_name, titles[k], _axis=_axis,
                journal=journal, y_ticks=ps_all, ax=ax,
                x_ticklabels=x_ticklabels, multi_image=True,
                panel_label=panel_labels[journal][k])

        # print statistics
        if verbose:
            print p.mean(), q.mean(), ma.median(p), ma.median(q)
            print stats.mannwhitneyu(ravel(p), ravel(q))

        del data, conservation, p, q

    return fig, matrices, array(boxplots)


# fig 4
def fig_conservation(
        generations=[1e6, 103309], positions=[1, 2, 4, 5],
        journal='plos', n_col=2, _save=True):

    fig_name = 'conservation'

    titles = ['no target', 'random']
    x_ticklabels = ['diagonal', 'off-diagonal']
    x_ticklabels = [', '.join(labels)
                    for labels in itertools.product(x_ticklabels, titles)]
    _axis = [positions[0] - 0.5, positions[-1] + 0.5, 0, 1]

    cmap = cm.RdBu
    _norm = Normalize(vmin=0, vmax=1)

    n_rows = 2; n_cols = 2
    fig_width = fig_widths[journal][n_col]
    fig       = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    # panels A and B: diagonal vs off-diagonal boxplots - do not plot
    matrices, boxplots = conservation_boxplots_2panel()[1:]
    # plot just one panel across the 2 columns
    ax = subplot2grid((n_rows, n_cols), (0,0), colspan=n_cols)

    boxplot(
        [boxplots[0,0], boxplots[1,0], boxplots[0,1], boxplots[1,1]],
        positions=positions)

    set_plot_kwargs(
        positions, '', fig_name, '', _axis=_axis, journal=journal,
        y_ticks=ps_all, ax=ax, x_ticklabels=x_ticklabels,
        multi_image=True, panel_label=panel_labels[journal][0])

    # panels C and D: average matrix heatmaps
    fig = multi_image2(
        matrices, titles, n_rows, n_cols, cmap, _norm, fig=fig, i_ax=3,
        cbar_label=fig_name, i_panel=1)

    if _save:
        save_fig(fig_name)#, get_figs_save_dir(journal, 'autoregulation'),
                 #fig_prefixes[journal], dpis[journal], formats[journal])

    return fig, matrices, boxplots


# 4) Coevolution of stability, robustness and autoregulation

# fig 5
def fig_coevolution(vspan=[[52, 2779], [1e5, 1e7]], threshold=.1, save_=True):
    return all_models_all_stats(
        ['p', 's', 'robustness', 'q', '2p'], 'm3r05b', threshold=threshold,
        y_ticks=ps_all, vspan=vspan, journal='plos',
        experiment='autoregulation', fig_name='coevolution-new', save_=save_)[0]


# panels 6A and 6B: boxplots
# P=1e3 * n=10 = samples=1e4 for each p
def fig_robustness(
        _step=.2, stability='stable', n_runs=10, p=.1, q=.5,
        fig=None, n_rows=1, n_cols=2, i_ax=0, title_='', _save=True):

    if stability == 'stable':
        fig_name = 'robustness'
    else:
        fig_name = 'robustness-random' #binary

    fname = ('robustness-qp.pop-%s-p%.1f-q%.2f_P1000-n%d' %
             (stability, p, q, n_runs))

    journal    = 'plos'; n_col = 2
    experiment = 'autoregulation'

    if not fig:
        fig_width = fig_widths[journal][n_col]
        fig       = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    # vs_p
    if p == .1:
        datapoints = ps_all
        x_ticks = arange(11)
        x_label = 'p'
        prefix  = ''
    # vs_q
    else:
        datapoints = qs_all
        x_ticks = arange(11)*9
        x_label = 'q'
        prefix  = 'q-'

    # panel A
    a = fig.add_subplot(n_rows, n_cols, 1+i_ax)

    data = ma.masked_invalid(load(data_dir + fname + '.npy'))
    data = data.reshape(datapoints.size, -1)
    #data.shape = (11, 10000)

    boxplot([x.compressed() for x in data])

    set_plot_kwargs(
        x_ticks+1, x_label, 'robustness', title_, journal=journal,
        y_ticks=ps_all, ax=a, x_ticklabels=datapoints[x_ticks],
        multi_image=True, panel_label=panel_labels[journal][i_ax])

    # panel B
    a = fig.add_subplot(n_rows, n_cols, 2+i_ax)
    n = shape(data)[1]
    p = [[x]*n for x in datapoints]
    all_data = ma.masked_invalid(append(p, data))
    all_data = all_data.reshape((2,-1))
    boxplot_data = get_percentiles(all_data, step=_step)
    # [ma.count(p) for p in all_data]
    # or [sum([where(x==p)[0].size for x in boxplot_data]) for p in datapoints]
    boxplot(boxplot_data)
    set_plot_kwargs(
        None, 'bins of increasing robustness', x_label, title_,
        journal=journal, y_ticks=ps_all, ax=a, x_ticklabels=['']*5,#None,
        multi_image=True, panel_label=panel_labels[journal][i_ax+1])

    if _save:
        save_fig(
            prefix + fig_name, get_figs_save_dir(journal, experiment),
            fig_prefixes[journal], dpis[journal], formats[journal])

    return fig, data, boxplot_data


# fig 6
def fig_robustness_vs_p(
        generations=[48, 1e6], title_='stable, non-evolved', _save=True):

    fig_name = 'robustness'

    journal = 'plos'; n_col = 2
    n_rows = 2; n_cols = 2

    fig_width = fig_widths[journal][n_col]
    fig       = figure(figsize = (fig_width*n_cols, fig_width*n_rows))

    # panels A and B:
    fig = fig_robustness(
        fig=fig, n_rows=n_rows, n_cols=n_cols, title_=title_, _save=False)[0]

    # panels C and D:
    fig = fig_robustness_vs_p_evolved(
        generations, n_rows, n_cols, fig=fig, i_ax=3, _save=False)


    if _save:
        save_fig(
            fig_name, get_figs_save_dir(journal, 'autoregulation'),
            fig_prefixes[journal], dpis[journal], formats[journal])

    return fig


# 5) Evolved networks are a distinct subset of stable networks

# fig 7
def fig_robustness_stable_vs_evolved(
        min_n=30, q=[.63, .55], generations=[956, 2.1e5], data=None,
        samples=100, loc=3, _save=True):

    fig_name = 'robustness_oversampling'

    journal = 'plos'; n_col = 2
    n_rows = 1; n_cols = 2

    fig_width = fig_widths[journal][n_col]
    fig       = figure(figsize = (fig_width*n_cols, fig_width*n_rows))

    # panel A:
    ax = fig.add_subplot(n_rows, n_cols, 1)
    fig, data = plot_evolution_oversampling(
        q[0], generations[0], data, min_n, samples, fig=fig, loc=loc,
        multi_image=True, panel=0, _save=False)[:2]

    # panel B:
    ax = fig.add_subplot(n_rows, n_cols, 2)
    fig, data = plot_evolution_oversampling(
        q[1], generations[1], data, min_n, samples, fig=fig, loc=loc,
        multi_image=True, panel=1, _save=False)[:2]

    if _save:
        save_fig(
            fig_name, get_figs_save_dir(journal, 'autoregulation'),
            fig_prefixes[journal], dpis[journal], formats[journal])

    return fig, data


# 6) Off-diagonal elements are nearly neutral at equilibrium

# fig 8
def fig_neutrality(
        error='within', threshold=.67, n_datapoints=15, ge=1e5, le=1e6,
        _save=True):

    _axis = (0, 1, 0, 1)
    ticks = arange(0, 1.1, .1)
    fig_name = 'neutrality'
    #fig_name = '-error_'.join(('neutrality', error))
    shortnames = mutbias_shortnames
    stats = ['p', '1-p']

    # no target
    fig = plot_stats_vs_parameter(
        stats, shortnames, 'mut_bias', error, ge, le, threshold, 'plos',
        _save=False)

    # random
    #plot(ticks, ticks, 'k--', label = 'neutrality')
    x = qp_biases

    y, yerr = get_qp_expection(
        'p', error, neutbias_shortnames, n_datapoints, threshold)

    fill_between(x, y-yerr, y+yerr, facecolor='k', alpha=.5)

    y, yerr = get_qp_expection(
        '1-p', error, neutbias_shortnames, n_datapoints, threshold)

    fill_between(x, y-yerr, y+yerr, facecolor='r', alpha=.5)

    set_plot_kwargs(
        ticks, 'no selection', 'selection', '', _axis, y_ticks=ticks,
        _grid=True)

    if _save:
        save_fig(fig_name)

    return fig

def fig_neutrality2(error=-1, n_datapoints=15, threshold=.67, _save=False):

    stats = ['p', '1-p']
    x = qp_biases

    _axis = (0, 1, 0, 1)
    ticks = arange(0, 1.1, .1)
    marker = markers[0]
    fig_name = '-error_std'.join(('neutrality', str(error)))
    journal = 'plos'; n_col = 2

    n_rows = 1; n_cols = 1
    fig_width = fig_widths[journal][n_col]
    fig       = figure(figsize = (fig_width*n_cols, fig_width*n_rows))

    for stat, color in zip(stats, colors[journal]):

        # m3: selection (no target)
        y, yerr = get_qp_expection(
            stat, error, mutbias_shortnames, n_datapoints, threshold)

        errorbar(
            x, y, yerr, fmt=marker, mfc=color, ecolor=color,
            label = get_stat_label(stat))

        # m2: no selections (random)
        y, yerr = get_qp_expection(
            stat, error, neutbias_shortnames, n_datapoints, threshold)

        fill_between(x, y-yerr, y+yerr, facecolor=color, alpha=.5)

    set_plot_kwargs(
        ticks, 'no selection', 'selection', '', _axis, y_ticks=ticks,
        _grid = True)

    if _save:
        save_fig(fig_name)

    return fig



#############
## SUP FIGS #
#############


# fig S1
# p_vs_stability - dont have the code..

# fig S2
def fig_p_evolution_all_models():
    return all_models_all_stats(
        'p', all_shortnames, fig_name='p_evolution-all_models', _save=True)[0]

# fig S3
# average_matrix - multi_heatmap_conservation() in plotting_tools2.py

# fig S4
# robustness-random
def fig_robustness_random(_step=.2, n_runs=30):
    return fig_robustness(_step, 'binary', n_runs)

# fig S5
# P=1e3 * n=10 = samples=1e4 for each p and q
def fig_robustness_vs_p_vs_q(
        qs=[.5, .55, .65, .7], n_cols=2, n_rows=2, _plot=True, _save=True):

    fig_name = 'robustness_vs_p_vs_q'
    journal  = 'plos'; n_col = 2
    experiment = 'autoregulation'

    if not _plot:
        fig_width = fig_widths[journal][n_col]
        fig       = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    for i, q in enumerate(qs):
        fname = 'robustness-qp.pop-stable-p0.1-q%.2f_P1000-n10' %q
        try:
            data = ma.masked_invalid(load(data_dir + fname + '.npy'))
            data = data.reshape(11, -1)
        except:
            pass
        else:
            if _plot:
                plot(ps_all, data.mean(axis=1), label='q = %.2f'%q)
            else:
                a = fig.add_subplot(n_rows, n_cols, i+1)
                boxplot(list(data))
                set_plot_kwargs(
                    range(1,12), 'p', 'robustness', 'q = %.2f'%q,
                    y_ticks=None, ax=a, x_ticklabels=ps_all,
                    multi_image=True, panel_label=panel_labels[journal][i])
                suffix = '-boxplots'
    if _plot:
        set_plot_kwargs(
            ps_all, 'p', 'robustness', y_ticks=None, ax=gca(),
            x_ticklabels = ps_all)

        suffix = '-plot'

    if _save:
        save_fig(
            fig_name + suffix, get_figs_save_dir(journal, experiment),
            fig_prefixes[journal], dpis[journal], formats[journal])

# fig S6 A and B
# P=1e3 * n=30 = samples=3e4 for each q
# p = 0.5
def fig_q_stability(
        fig=None, n_rows=1, n_cols=2, i_ax=0, _save=True, normed=True):

    fig_name = 'q-stability'
    journal = 'plos'; n_col = 2
    x_ticklabels = ['unstable', 'stable']

    a_axis = [.5, 2.5, -.05, 1.05]
    b_axis = [-.01, 1.01, 0, 3e4]

    if not fig:
        fig_width = fig_widths[journal][n_col]
        fig       = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    # load data
    fname = 's-qp.pop-binary-p0.5-q0.01_P1000-n30'
    data = ma.masked_invalid(load(data_dir + fname + '.npy'))
    data = data.reshape(qs_all.size, -1)
    #data.shape = (91, 30000)
    data, boxplot_data = get_boxplot_data(data=data, norm=90.)

    if normed:
        _norm = float(data.shape[1])
        y_label = 'stability'
    else:
        _norm = 1
        y_label = 'number of stable matrices, $n_f$'

    # panel A: q
    ax = fig.add_subplot(n_rows, n_cols, 1+i_ax)

    boxplot(boxplot_data)

    set_plot_kwargs(
        [1,2], '', 'q', _axis=a_axis, journal=journal, y_ticks=ps_all, ax=ax,
        x_ticklabels=x_ticklabels,
        multi_image=True, panel_label=panel_labels[journal][0])

    # panel B: stability
    ax = fig.add_subplot(n_rows, n_cols, 2+i_ax)

    hist(boxplot_data[1], append(qs_all, 1.01), align='left')

    y_ticklabels = ['%.1f'%y for y in ax.get_yticks()/_norm]

    set_plot_kwargs(
        ps_all, 'q', y_label, _axis=b_axis, journal=journal, ax=ax,
        y_ticklabels=y_ticklabels,
        multi_image=True, panel_label=panel_labels[journal][1])

    del data, boxplot_data

    if _save:
        save_fig(
            fig_name, get_figs_save_dir(journal, 'autoregulation'),
            fig_prefixes[journal], dpis[journal], formats[journal])

    return fig

# fig S6
def fig_q_stability_robustness(_save=True):

    fig_name = 'q-stability-robustness'
    journal = 'plos'; n_col = 2
    n_rows = 2; n_cols = 2

    fig_width = fig_widths[journal][n_col]
    fig       = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    # panel A and B
    fig = fig_q_stability(fig, n_rows, n_cols, 0, False)
    # panel C and D
    # fig_robustness(step, stability, n_runs, p, q, fig, n_rows, n_cols, i_ax)
    fig = fig_robustness(
        .1, 'stable', 10, .5, .01, fig, n_rows, n_cols, 2, '', False)[0]

    if _save:
        save_fig(
            fig_name, get_figs_save_dir(journal, 'autoregulation'),
            fig_prefixes[journal], dpis[journal], formats[journal])

    return fig

# fig S7
def q_bimodal(
        data=None, i_generations=[63, 94], _axis=(0,1,0,1e4), _vlines=False,
        _save=True):

    fig_name = 'q-bimodal'

    shortname = 'm3r05b'
    y_labels = ['number of evolved matrices', '']

    journal = 'plos'; n_col = 2

    n_rows = 1; n_cols = 2
    fig_width = fig_widths[journal][n_col]
    fig       = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    if not any(data):
        data = get_stat_vs_g(shortname, '1-p')

    for i, g in enumerate(i_generations):
        ax = fig.add_subplot(n_rows, n_cols, i+1)

        n, bins, pa = hist(data[g], qs_all)

        if _vlines:
            statistics = [data[g].mean(), stats.mstats.mode(data[g])[0]]
            _colors = colors[journal][:2]
            labels = ['mean', 'mode']
            for x, color, label in zip(statistics, _colors, labels):
                ax.axvline(x, c=color, ls='solid', label=label)

        title_ = 'generation %s' %print_generations(
            pop_generations[g], False)[2:]

        set_plot_kwargs(
            ps_all, 'q', y_labels[i], title_, _axis, ax=ax,
            panel_label=panel_labels[journal][i])

    if _save:
        save_fig(
            fig_name, get_figs_save_dir(journal, 'autoregulation'),
            fig_prefixes[journal], dpis[journal], formats[journal])

    return fig, data

# fig S8
def fig_q_mean_median_mode():
    return plot_mean_median_mode('q')

# fig S9
def fig_p_q_vs_period():
    return plot_stats_vs_parameter(['p','1-p'], fig_name='p_q_vs_period')

# fig S10
def fig_engineering_robustness(_save=True):

    fig_name = 'engineering_robustness'
    pp = [ (.1,.9), (.1, 1),   (0,.9),   (0,1)]
    pq = [(.9,.53), (1,.53), (.9,.54), (1,.54)]

    journal = 'plos'; n_col = 2
    n_rows = 1; n_cols = 2

    fig_width = fig_widths[journal][n_col]
    fig       = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    # panel A: superfrakenmatrix
    a = fig.add_subplot(n_rows, n_cols, 1)

    fnames = ('robustness-stable.pop-binary-'
             + 'p%.1f-p%.1f-diagonalp_P1000-n200.npy' %(x[0], x[1])
             for x in pp)

    data = [load(data_dir + fname).mean(axis=1) for fname in fnames]

    boxplot(data)

    set_plot_kwargs(
        range(1,5), '(p, q)', 'robustness', 'engineered networks',
        journal=journal, y_ticks=ps_all, ax=a, x_ticklabels=map(str, pq),
        multi_image=True, panel_label=panel_labels[journal][0])

    # panel B: random stable
    a = fig.add_subplot(n_rows, n_cols, 2)
    data = [get_non_evolved_pop(
        'robustness', False, 200, None, p, q, extension='.npy')[0]
            for p, q in pq]

    boxplot(data)

    set_plot_kwargs(
        range(1,5), '(p, q)', 'robustness', 'random networks',
        journal=journal, y_ticks=ps_all, ax=a, x_ticklabels=map(str,pq),
        multi_image=True, panel_label=panel_labels[journal][1])

    if _save:
        save_fig(
            fig_name, get_figs_save_dir(journal, 'autoregulation'),
            fig_prefixes[journal], dpis[journal], formats[journal])

    return fig

# fig S11
def fig_qp_starting_conditions(threshold=.61, _save=True):
    _axis = (1, 3.2e6, 0, 1.02)
    fig_name = 'starting_conditions'
    stats = ['p', '1-p', 's', 'robustness']
    x_label = 'generations'
    fig = one_stat_per_panel(
        start_shortnames, stats, 2, 2, journal='plos', n_col=2,
        shortname_label=get_model_label, threshold=threshold, _axis=_axis,
        x_label=x_label, y_ticks=ps_all, _grid=True, _save=_save,
        fig_name=fig_name)[0]
    return fig

# fig S12
def fig_p_r_vs_period(_save=False):
    return plot_stats_vs_parameter(
        ['p','equal'], _save=_save, fig_name='p_r_vs_period')

# fig S13
def fig_recomb(threshold=.49, _save=True):
    return fig_mutation_rate(
        shortnames=['m3r05b', 'm3r0b'],
        _stats=['p', '1-p', 's', 'robustness'],
        threshold=threshold, shortname_label=get_model_label_rec_text,
        _grid=True, n_col=2, fig_name='no-recomb', _save=_save)

# fig S14
def fig_density(threshold=.49, _save=True):
    return fig_mutation_rate(
        shortnames=['m3r05b', 'm3r05k2', 'm3r0k2'],
        _stats=['p', '1-p', 's', 'robustness'],
        threshold=threshold, shortname_label=get_model_label_density_text,
        _grid=True, n_col=2, fig_name='k2', _save=_save)

# fig S15
def cycling_genes(samples=1e5, _save=True):

    fig_name = 'cycling_genes'
    labels = ['cycles of size 2', 'all cycles']

    journal = 'plos'; n_col = 2

    n_rows = 1; n_cols = 1
    fig_width = fig_widths[journal][n_col]
    fig       = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    cycling_genes = get_number_cycling_genes(int(samples))
    x = [cycling_genes[2], list(itertools.chain.from_iterable(cycling_genes))]
    bins = arange(1,12)

    hist(x, bins, align='left', color=colors[journal][:2], label=labels)

    set_plot_kwargs(
        range(1,11), 'number of cycling genes', 'number of matrices',
        _axis = (.5, 10.5, 0, axis()[3]),
        #labels = labels, n_hist = 2, n_bins = bins.size,
        journal = journal)

    if _save:
        save_fig(
            fig_name, get_figs_save_dir(journal, 'autoregulation'),
            fig_prefixes[journal], dpis[journal], formats[journal])

    return fig

# sup fig
def fig_mutation_rate(
        shortnames=m3u_shortnames, _stats=['p', 'q', 's', 'robustness'],
        threshold=.49, shortname_label=get_model_label_mut_text,
        _grid=True, n_col=2, fig_name='mutation_rate', _save=True):

    return one_stat_per_panel(
        shortnames, _stats, 2, 2, journal='plos', n_col=n_col,
        stat_label=get_stat_label, shortname_label=shortname_label,
        threshold=threshold, _vlines=None, _axis=None,
        x_label='', y_ticks=ps_all, _grid=_grid, _save=_save,
        fig_name=fig_name)[0]

# sup fig
def fig_selection_rate(threshold=.49, _save=True):
    return fig_mutation_rate(
        shortnames=m4s_shortnames, _stats=['p', 'q', 's', 'robustness'],
        threshold=threshold, shortname_label=get_model_label_sel_text,
        _grid=True, n_col=2, fig_name='selection_rate', _save=_save)

# sup fig
def fig_qp_initial_var(threshold=.49, _save=True):
    return fig_mutation_rate(
        shortnames=['m3r05b', 'b05p05q05', 'm3r05d'],
        _stats=['p', '1-p', 's', 'robustness'],
        threshold=threshold, shortname_label=get_model_label_init_var,
        _grid=True, n_col=2, fig_name='qp_initial_var', _save=_save)

# sup fig
def fig_q_vs_random(i=89, _norm=90., _save=True):

    fig_name = 'q_vs_random'
    labels = ['no target', 'random']

    journal = 'plos'; n_col = 2
    n_rows = 1; n_cols = 2

    fig_width = fig_widths[journal][n_col]
    fig       = figure(figsize = (fig_width*n_cols, fig_width*n_rows))

    # data
    m3 = get_stat_vs_g('m3r05b', 'q')[i]/_norm
    m2 = get_stat_vs_g('m2r05b', 'q')[i]/_norm

    # panel A: hist
    a = fig.add_subplot(n_rows, n_cols, 1)
    n, bins, pa = hist([m3, m2], normed=True)
    set_plot_kwargs(
        None, 'q', '', labels=labels, n_hist=2, n_bins=10, journal='plos',
        ax=a, multi_image=True, panel_label=panel_labels[journal][0])

    # panel B: boxplot
    a = fig.add_subplot(n_rows, n_cols, 2)
    b = boxplot([m3, m2])
    set_plot_kwargs(
        [1,2], '', 'q', _axis=[.5, 2.5, 0, 1], journal=journal,
        y_ticks=ps_all, ax=a, x_ticklabels=labels,
        multi_image=True, panel_label=panel_labels[journal][1])

    print ma.median(m3), ma.median(m2), stats.mannwhitneyu(m3, m2)

    if _save:
        save_fig(fig_name)#, get_figs_save_dir(journal, 'autoregulation'),
                 #fig_prefixes[journal], dpis[journal], formats[journal])

    return fig, [m3, m2]

# sup fig
def fig_robustness_vs_p_evolved(
        generations=[48, 1e6], n_rows=1, n_cols=2, _stats=['p', 'robustness'],
        min_n=100, fig=None, i_ax=1, _plotting=True, _save=True):

    shortname = 'm3r05b'

    data = get_stats_vs_g(shortname, _stats).swapaxes(0,1)

    if not _plotting:
        return data

    fig_name = '%s_vs_%s-evolved' %(_stats[1], _stats[0])
    journal = 'plos'; n_col = 2

    if _stats[0] == 'p':
        bins = ps_all
        x_step = 1
        _norm = 1
        f = lambda x: '%.1f' %x

    elif _stats[0] == 'q':
        bins = arange(91)
        x_step = 4
        _norm = 90.
        f = lambda x: '%.2f' %x

    if not fig:
        fig_width = fig_widths[journal][n_col]
        fig       = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    i_x = [where(pop_generations == g)[0][0] for g in generations]

    # panels A and B:
    for k, i in enumerate(i_x):
        a = fig.add_subplot(n_rows, n_cols, k+i_ax)

        percentiles = get_percentiles(data[i], bins, discrete=True)

        j = where(array(map(len,percentiles)) > min_n)[0]

        b = boxplot(percentiles[j])

        set_plot_kwargs(
            range(1,j.size+1,x_step), _stats[0], _stats[1],
            'generation %s' %print_generations(generations[k], False)[2:],
            journal=journal, ax=a, x_ticklabels=map(f,bins[j][::x_step]/_norm),
            multi_image=True, panel_label=panel_labels[journal][k+i_ax-1])

        del a, percentiles, j, b

    del data

    if _save:
        save_fig(
            fig_name, get_figs_save_dir(journal, 'autoregulation'),
            fig_prefixes[journal], dpis[journal], formats[journal])
    return fig

def qpgr(
        fname='qpr.matrix-m3r05b', mask='qpr.matrix.masks-m3r05b',
        cmap=None, _norm=None, fig=None, generations=[48, 1e6],
        n_rows=1, n_cols=2, i_ax=1, _save=True):

    if not cmap:
        cmap = cm.jet
    if not _norm:
        norm = Normalize(vmin=0, vmax=1)

    fig_name = fname.split('-')[0]
    journal = 'plos'; n_col = 2

    if not fig:
        fig_width = fig_widths[journal][n_col]
        fig       = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    # get qp data
    data = ma.masked_invalid(load(data_dir + fname + '.npy'))
    if mask: data.mask = load(data_dir + mask + '.npy')
    i_x = [where(pop_generations == g)[0][0] for g in generations]

    # panels A and B
    fig = multi_image2(
        data[i_x], ['']*2, n_rows, n_cols, cmap, _norm, fig=fig, i_ax=i_ax,
        x_label='q', y_label='p', x_step=9,
        x_ticklabels=qs_all[arange(11)*9], y_ticklabels=ps_all,
        journal=journal)

    if _save:
        save_fig(
            fig_name, get_figs_save_dir(journal, 'autoregulation'),
            fig_prefixes[journal], dpis[journal], formats[journal])

    return fig

# sup fig 11
def fig_robustness_vs_qp(generations=[48, 1e6], _save=True):

    fig_name = 'robustness_vs_qp'

    cmap = cm.jet
    _norm = Normalize(vmin=0, vmax=1)

    journal = 'plos'; n_col = 2
    n_rows = 2; n_cols = 2

    fig_width = fig_widths[journal][n_col]
    fig       = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    # panels A and B:
    fig = fig_robustness_vs_p_evolved(
        generations, n_rows, n_cols, fig=fig, _save=False)

    # panels C and D:
    fig = fig_robustness_vs_p_evolved(
        generations, n_rows, n_cols, _stats=['q', 'robustness'],
        fig=fig, i_ax=3, _save=False)

    # panels E and F:
    #fig = qpgr(fig = fig, generations = generations, n_rows = n_rows,
     #          n_cols = n_cols, i_ax = 5, _save = False)

    if _save:
        save_fig(fig_name, get_figs_save_dir(journal, 'autoregulation'),
                 fig_prefixes[journal], dpis[journal], formats[journal])
    return fig

# sup fig
def fig_stability_vs_p_evolved(
        _stats=['s', 'p'], generations=[48, 1e6],
        n_rows=1, n_cols=2, fig=None, i_ax=1, _plotting=True, _save=True):

    shortname = 'm3r05b'

    data = get_stats_vs_g(shortname, _stats).swapaxes(0,1)

    if not _plotting:
        return data

    fig_name = '%s_vs_%s-evolved' %(
        get_stat_label(_stats[0]), get_stat_label(_stats[1]))
    journal = 'plos'; n_col = 2

    if not fig:
        fig_width = fig_widths[journal][n_col]
        fig       = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    i_x = [where(pop_generations == g)[0][0] for g in generations]

    # panels A and B:
    for k, i in enumerate(i_x):
        a = fig.add_subplot(n_rows, n_cols, k+i_ax)
        boxplot_data = get_boxplot_data2(data=data[i])
        n, bins, pa = hist(boxplot_data[1], array(range(12))/10., align='left')

        y_label = 'number of stable matrices' if k==0 else ''
        set_plot_kwargs(
            ps_all, get_stat_label(_stats[1]), y_label,
            'generation %s' %print_generations(generations[k], False)[2:],
            (0, 1.05, 0, n.max()+n.min()),
            journal=journal, ax=a, x_ticklabels=ps_all,
            multi_image=True, panel_label=panel_labels[journal][k])

        del a, boxplot_data, n, bins, pa

    del data, i_x

    if _save:
        save_fig(fig_name)#, get_figs_save_dir(journal, 'autoregulation'),
                 #fig_prefixes[journal], dpis[journal], formats[journal])
    return fig

# sup fig
def fig_gradient(shortname='m3r05b', f=ma.median, _save=True):

    fig_name = 'gradient'
    journal  = 'plos'
    qps, qpr = load_qpgsr_matrices(shortname)[1:]
    s = [f(gradient_norm(x)) for x in qps]
    r = [f(gradient_norm(x)) for x in qpr]
    fig = figure()
    plot(pop_generations, transpose([s,r]))
    set_plot_kwargs(
        None, 'generations', f.__name__ + ' gradient', x_scale='log',
        labels=['stability', 'robustness'], journal='plos')

    if _save:
        save_fig(
            fig_name, get_figs_save_dir(journal, 'autoregulation'),
            fig_prefixes[journal], dpis[journal], formats[journal])
    return fig

# sup fig
def fig_qp_mut_bias_2panel():
    _axis      = (1, 1e7, 0, 1)
    y_ticks    = arange(0,1.1,.1)
    fig_name   = 'qp.mut_bias-2panel'
    shortnames = mutbias_shortnames
    stats      = ['p', '1-p']
    fig, axess = one_stat_per_panel(
        shortnames, stats, 2, 1, journal=None, n_col=1,
        shortname_label=get_model_label_mut_bias, _axis=_axis, y_ticks=y_ticks,
        _grid=True, _save=True, fig_name=fig_name)
    return fig

def fig_q_histograms():

    fig_name = 'q_histograms'

    journal = 'plos'; n_col = 2

    n_rows = 5; n_cols = 3
    fig_width = fig_widths[journal][n_col]
    fig       = figure(figsize = (fig_width*n_cols, fig_width*n_rows))

    for i, shortnames in enumerate(zip(
            mutbias_shortnames, neutbias_shortnames)[::-1]):
        a = fig.add_subplot(n_rows, n_cols, i+1)
        plot_q_histogram(shortnames, '1-p', ma.mean, fig)

    save_fig(fig_name)
    return fig
