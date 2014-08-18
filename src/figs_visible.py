import src.plotting_tools2
reload(src.plotting_tools2)
from src.plotting_tools2 import *

# from https://github.com/yhat/ggplot
#from ggplot import *

reload(plos)
journal = plos
experiment = 'visible'

################
# OLD FIGURES ##
################

def old_figure1():
    return stable_states_vs_N()

def old_figure2AB(_save=True):

    fig_name = 'stable_states_vs_k'
    x_label = '# gene expression states, k'
    y_label = '# phenotypes'
    label = 'visible'
    titles = ['binary matrices', 'real matrices']
    fnames = ['N4bun8inf', 'N4rsn8inf']

    leg_loc = 2 # upper left
    panel_loc = 1 # upper right
    journal = 'plos'; n_col = 2
    fig_width = fig_widths[journal][n_col]
    n_rows = 1; n_cols = 2
    fig = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    i = 0
    ax = fig.add_subplot(n_rows, n_cols, i+1)

    fig, y = new_stable_states_vs_k(
        fnames[i], fig=fig, potential=True, ls='o', label=label)

    set_plot_kwargs(
        None, x_label, y_label, titles[i], journal=journal, ax=ax, loc=leg_loc,
        panel_label=panel_labels[journal][i], panel_label_loc=panel_loc)

    i = 1
    ax = fig.add_subplot(n_rows, n_cols, i+1)

    fig, y = new_stable_states_vs_k(
        fnames[i], fig=fig, potential=True, ls='o', label=label)

    set_plot_kwargs(
        None, x_label, y_label, titles[i], journal=journal, ax=ax, loc=leg_loc,
        panel_label=panel_labels[journal][i], panel_label_loc=panel_loc)

    if _save:
        save_fig(
            fig_name, get_figs_save_dir(journal, experiment),
            fig_prefixes[journal], dpis[journal], formats[journal])

    return fig

#def old_figure2():
#    fig_name = 'stable_states_vs_k'
#    fnames = ['N4bun8inf', 'N4rsn8inf']
#    labels = ['binary', 'real']
#    return stable_states_vs_k_vs_all(
#        fnames, colors[None], labels=labels, _title=' ',
#        fig_name=fig_name)

def old_figure3(
        fnames=myrun_noiseN4s, n=None, devo_times=[], samples=[],
        fig=None, noise_lbl='noise', pot_lbl='potential',
        verbose=False, save_=False):

    if n: n = 10**n

    # set_plot_kwargs
    experiment = 'transition'
    #fig_name = ''
    myrun = cPickle.load(open(logs_dir + fnames[0] + '.run'))
    title_ = get_myrun_title(myrun, experiment, n)

    if not fig:
        fig = figure()

    fig, y = stable_states_vs_k_vs_all(
        fnames, n, devo_times, samples, noise_label=noise_lbl, pot_lbl=pot_lbl,
        title_=title_, fig=fig, verbose=verbose, save_=save_)

    return fig

def old_figure4(_save=True):

    fig_name = 'meta'
    x_label = '# gene expression states, k'
    y_label = '# phenotypes'
    #titles = []
    fnames = [
        myrun_noiseN4u, myrun_noiseN4f, myrun_noiseN4k2,
        myrun_noiseN2s]#, myrun_noiseN3s, myrun_noiseN5s]

    leg_loc = 2 # upper left
    panel_loc = 1 # upper right
    journal = 'plos'; n_col = 2
    fig_width = fig_widths[journal][n_col]
    n_rows = 3; n_cols = 2
    fig = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    for i, fname in enumerate(fnames):
        ax = fig.add_subplot(n_rows, n_cols, i+1)

        stable_states_vs_k_vs_all(fname, fig=fig, _save=False)

        set_plot_kwargs(
            None, x_label, y_label, journal=journal, ax=ax, loc=leg_loc,
            panel_label=panel_labels[journal][i], panel_label_loc=panel_loc)

    if _save:
        save_fig(
            fig_name, get_figs_save_dir(journal, experiment),
            fig_prefixes[journal], dpis[journal], formats[journal])

    return fig

def figureS1(_save=True):
    #fig_name =
    fnames = ['N4bun8inf']#myrun_fnamesN4binary
    x_label = '# gene expression states, k'
    labels = ['']
    _title = ' '
    color = colors[None]
    return percent_visible_states(
        fnames, color, x_label, labels, _title, _save=_save)

def figureS2(_save=True):
    fig_name = 'stability_vs_k_vs_noise'
    fnames = myrun_stabilityN4
    return stable_states_vs_k_vs_all(fnames, fig_name=fig_name, _save=_save)

def figureS3(_save=True):

    fig_name = 'overrepresentation'
    x_label = r'gene expression state, $s_i$'
    y_label = '% visible'
    #titles = []
    ks = [3, 4, 6, 7]

    #leg_loc = 2 # upper left
    panel_loc = 1 # upper right
    journal = 'plos'; n_col = 2
    fig_width = fig_widths[journal][n_col]
    n_rows = 2; n_cols = 2
    fig = figure(figsize=(fig_width*n_cols, fig_width*n_rows))

    for i, k in enumerate(ks):
        ax = fig.add_subplot(n_rows, n_cols, i+1)

        stable_states_overrepresentation(k, fig=fig, _save=False)

        set_plot_kwargs(
            None, x_label, y_label, journal=journal, ax=ax,# loc=leg_loc,
            panel_label=panel_labels[journal][i], panel_label_loc=panel_loc)

    if _save:
        save_fig(
            fig_name, get_figs_save_dir(journal, experiment),
            fig_prefixes[journal], dpis[journal], formats[journal])

    return fig


#########################
# NEW FIGURES - VICTOR ##
#########################

def savefig(fig, figname):
    # NOTE: savefig is ignoring rcParams, don't know why..
    # switching to manual as a workaround
    _dir = matplotlib.rcParams['savefig.directory']
    dpi = matplotlib.rcParams['savefig.dpi']
    _format = matplotlib.rcParams['savefig.format']
    _bbox = matplotlib.rcParams['savefig.bbox']
    _pad_inches = matplotlib.rcParams['savefig.pad_inches']
    fname = '.'.join((figname, _format))
    fig.savefig(_dir + fname, dpi=dpi, format=_format, bbox_inches=_bbox,
                pad_inches=_pad_inches)

def figure1(fnames=[myrun_noiseN4m0u, myrun_noiseN4m0s], n=log10(102353),
            devo_times=[], verbose=False, save_=True):

    experiment = 'transition'
    figname = 'figure1'
    x_label = 'N'
    y_label = r'$\mathcal{U}$'
    #noise_lbl = '\sigma' # no need to add r'$$' (already in get_noise_label())
    samp_lbl = r'$\lambda$'
    pot_lbl = r'$\mathcal{U}_{max} = 2^N$'
    titles = ['random', 'stable']

    leg_loc = 'best'
    panel_loc = -1 # upper left outside the box
    n_rows = 1; n_cols = 2

    set_custom_rcParams(journal, experiment, n_cols)
    fig = figure()

    for i, fname in enumerate(fnames):
        ax = fig.add_subplot(n_rows, n_cols, i+1)

        fig = stable_states_vs_N_vs_noise(
            fnames[i], 10**n, devo_times=devo_times, fig=fig,
            pot_lbl='', noise_lbl='', samp_lbl='', verbose=verbose)

        axis((3,20,10,1e6))
        xlabel(x_label)
        ylabel(y_label)
        title(titles[i])
        legend(loc='best', fancybox=True, title=r'$\sigma$')
        set_panel_labels(ax, uppercase[i], panel_loc)
        remove_spines(ax)
        remove_ticks(ax)
        ax.set_axis_bgcolor("#E5E5E5")

        # annotated label of 'potential' dashed line
        x = ax.get_children()[2].get_xydata()[-1]
        ax.annotate(pot_lbl, (x[0],x[1]), xytext=(0.7, 0.92),
                    textcoords='axes fraction', va='center')

        # annotated label of 'n_samples' dashed line
        x = ax.get_children()[3].get_xydata()[-1]
        y = find_nearest(get_int_datapoints(0, log10(1e8), 100), 10**n)
        ax.annotate(samp_lbl, (x[0],x[1]), xytext=(10, y),
                    textcoords='data', va='bottom')

    if save_:
        savefig(fig, figname)

    reset_default_rcParams()
    return fig

def figure2(
        fnames=[myrun_noiseN4m0u, myrun_noiseN4m0s], dim=3,
        devo_times=[], hline=True, verbose=True, save_=True):

    experiment = 'transition'
    figname = 'figure2'
    _format = '.eps'

    panel_loc = -1 # upper left outside the box
    n_rows = 1; n_cols = 2

    set_custom_rcParams(journal, experiment, n_cols)
    matplotlib.rcParams['font.size'] = 14

    fig = figure()

    for i, fname in enumerate(fnames):
        ax = fig.add_subplot(n_rows, n_cols, i+1)

        fig = plot_access_vs_samples_vs_noise(
            fname, dim, devo_times, hline, fig, verbose)

        legend(loc='upper center', fancybox=True, title=r'$\sigma$')
        set_panel_labels(ax, uppercase[i], panel_loc)
        ggplot.theme_gray().apply_theme(ax)
        ax.set_axis_bgcolor("#E5E5E5")
        remove_spines(ax)
        remove_ticks(ax)

    if save_:
        #fname = '.'.join((figname, _format))
        savefig(fig, figname)
        #fig.savefig(_dir + fname, dpi=dpi, format=_format, bbox_inches=_bbox,
         #           pad_inches=_pad_inches)

    reset_default_rcParams()

    return fig

def figure3(
        fnames=[myrun_noiseN4m0u, myrun_noiseN4m0s], n=log10(2.42013e+06),
        devo_times=[], verbose=False, save_=True):

    experiment = 'transition'
    figname = 'figure3'
    titles = ['random', 'stable']
    _format = '.eps'

    panel_loc = -1 # upper left outside the box
    n_rows = 1; n_cols = 2

    set_custom_rcParams(journal, experiment, n_cols)
    matplotlib.rcParams['font.size'] = 14

    fig = figure()

    for i, fname in enumerate(fnames):
        ax = fig.add_subplot(n_rows, n_cols, i+1)
        fig = old_figure3(
            fname, n, devo_times,
            fig=fig, noise_lbl='', pot_lbl='', verbose=verbose)

        # set_plot_kwargs
        # set panel B's axis as A's
        if i == 0:
            _axis = axis()
        else:
            axis(_axis)
        xlabel(r'$k$')
        ylabel(r'$\mathcal{U}$')
        title(titles[i])
        legend(loc='best', fancybox=True, title=r'$\sigma$')
        set_panel_labels(ax, uppercase[i], panel_loc)
        ggplot.theme_gray().apply_theme(ax)
        ax.set_axis_bgcolor("#E5E5E5")
        remove_spines(ax)
        remove_ticks(ax)

    if save_:
        #fname = '.'.join((figname, _format))
        savefig(fig, figname)
        #fig.savefig(_dir + fname, dpi=dpi, format=_format, bbox_inches=_bbox,
         #           pad_inches=_pad_inches)

    reset_default_rcParams()

    return fig

def figure4(
        fnames=[myrun_noiseN4m0realu, myrun_noiseN4m0reals],
        n=log10(1.6681e+06), verbose=False, save_=True):

    experiment = 'transition'
    figname = 'figure4'
    titles = ['random', 'stable']
    _format = '.eps'

    panel_loc = -1 # upper left outside the box
    n_rows = 1; n_cols = 2

    set_custom_rcParams(journal, experiment, n_cols)

    # this seems to only affect font size in legend
    matplotlib.rcParams['font.size'] = 14

    fig = figure()

    for i, fname in enumerate(fnames):
        ax = fig.add_subplot(n_rows, n_cols, i+1)

        fig = stable_states_vs_k_vs_noise(
            fname, n, fig=fig, noise_lbl='', pot_lbl='', verbose=verbose)

        # annotate potential curve
        k,N = 9,4
        umax = k**N
        label = r'$\mathcal{U}_{max} = k^N$'
        annotate(label, (k,umax), xytext=(6.9, 5000),
                 textcoords='data', va='bottom')

        # set_plot_kwargs
        ax.axis((1,10,10,1e4))
        ax.set_xticks(range(2,10))
        xlabel(r'$k$')
        ylabel(r'$\mathcal{U}$')
        title(titles[i])
        legend(loc='best', fancybox=True, title=r'$\sigma$')
        set_panel_labels(ax, uppercase[i], panel_loc)
        ggplot.theme_gray().apply_theme(ax)
        ax.set_axis_bgcolor("#E5E5E5")
        remove_spines(ax)
        remove_ticks(ax)

    if save_:
        #fname = '.'.join((figname, _format))
        savefig(fig, figname)
        #fig.savefig(_dir + fname, dpi=dpi, format=_format, bbox_inches=_bbox,
         #           pad_inches=_pad_inches)

    reset_default_rcParams()

    return fig

def figure5(
        fnames=[myrun_noiseN4m0nfu, myrun_noiseN4m0rfu],
        devo_times=[], verbose=False, save_=True):

    experiment = 'transition'
    figname = 'figure5'
    _format = '.eps'

    samp_lbl = r'$\lambda$'
    pot_lbl = r'$\mathcal{U}_{max} = 2^N$'
    titles = ['neighbor flip', 'random or no flip']

    leg_loc = 'best'
    panel_loc = -1 # upper left outside the box
    n_rows = 2; n_cols = 2

    set_custom_rcParams(journal, experiment, n_cols)
    matplotlib.rcParams['font.size'] = 14
    figsize = [12.0, 6.0*n_rows]
    fig = figure(figsize=figsize)

    # from figure1()
    n = log10(102353)
    for i, fname in enumerate(fnames):
        ax = fig.add_subplot(n_rows, n_cols, i+1)

        fig = stable_states_vs_N_vs_noise(
            fnames[i], 10**n, devo_times=devo_times, fig=fig,
            pot_lbl='', noise_lbl='', samp_lbl='', verbose=verbose)

        axis((3,20,10,1e6))
        xlabel('N')
        ylabel(r'$\mathcal{U}$')
        title(titles[i])
        legend(loc='best', fancybox=True)#, title=r'$\sigma$')
        set_panel_labels(ax, uppercase[i], panel_loc)
        remove_spines(ax)
        remove_ticks(ax)
        ax.set_axis_bgcolor("#E5E5E5")

        # annotated label of 'potential' dashed line
        x = ax.get_children()[2].get_xydata()[-1]
        ax.annotate(pot_lbl, (x[0],x[1]), xytext=(0.7, 0.92),
                    textcoords='axes fraction', va='center')

        # annotated label of 'n_samples' dashed line
        x = ax.get_children()[3].get_xydata()[-1]
        y = find_nearest(get_int_datapoints(0, log10(1e8), 100), 10**n)
        ax.annotate(samp_lbl, (x[0],x[1]), xytext=(10, y),
                    textcoords='data', va='bottom')


    # from figure3()
    n = log10(2.42013e+06)
    for i, fname in enumerate(fnames):
        ax = fig.add_subplot(n_rows, n_cols, i+3)
        fig = old_figure3(
            fname, n, devo_times,
            fig=fig, noise_lbl='', pot_lbl='', verbose=verbose)

        # set_plot_kwargs
        # set panel B's axis as A's
        if i == 0:
            _axis = axis()
        else:
            axis(_axis)
        xlabel(r'$k$')
        ylabel(r'$\mathcal{U}$')
        title(titles[i])
        legend(loc='best', fancybox=True)#, title=r'$\sigma$')
        set_panel_labels(ax, uppercase[i], panel_loc)
        ggplot.theme_gray().apply_theme(ax)
        ax.set_axis_bgcolor("#E5E5E5")
        remove_spines(ax)
        remove_ticks(ax)

    if save_:
        #fname = '.'.join((figname, _format))
        savefig(fig, figname)
        #fig.savefig(_dir + fname, dpi=dpi, format=_format, bbox_inches=_bbox,
         #           pad_inches=_pad_inches)

    reset_default_rcParams()

    return fig
