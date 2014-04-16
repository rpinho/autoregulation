import src.plotting_tools2
reload(src.plotting_tools2)
from src.plotting_tools2 import *

experiment = 'visible'

def figure1():
    return stable_states_vs_N()

def figure2AB(_save=True):

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

#def figure2():
#    fig_name = 'stable_states_vs_k'
#    fnames = ['N4bun8inf', 'N4rsn8inf']
#    labels = ['binary', 'real']
#    return stable_states_vs_k_vs_all(
#        fnames, colors[None], labels=labels, _title=' ',
#        fig_name=fig_name)

def figure3(_save=True):
    #fig_name = ''
    fnames = myrun_noiseN4s
    fig, y = stable_states_vs_k_vs_all(fnames, _save=_save)

def figure4(_save=True):

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
