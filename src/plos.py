### MATPLOTLIBRC FORMAT for PLoS journals
# from http://matplotlib.org/users/customizing.html

#TODO: consider creating local matplotlibrc in the current working directory

### LINES
lines = {'linewidth'   : 1.5, # line width in points
         'markersize'  : 8    # markersize, in points
}
### PATCHES
patch = {#'facecolor'   : "348ABD"
}
### FONT
font = {'family'       : 'Times New Roman', #sans-serif
        'weight'       : 'medium',#normal',#bold'
        'size'         : 16 #12.0
        #'serif'        : 'Times'
}
### TEXT
### AXES
axes = {#'facecolor'    : '#E5E5E5' #white   # axes background color
}
### TICKS
### GRIDS
grid = {#'color'        : 'white', #black   # grid color
        #'linestyle'    : 'solid', # dotted
        #'linewidth'    :  1.4# 0.5     # in points
#grid.alpha       :   1.0     # transparency, between 0.0 and 1.0
}
### LEGEND
legend = {'fancybox'   : True, #False  # if True, use a rounded box for the
                               # legend, else a rectangle
          'fontsize'   : 'medium'
# Special text sizes can be defined
# relative to font.size, using the following values: xx-small, x-small,
# small, medium, large, x-large, xx-large, larger, or smaller
}
### FIGURE
figure = {'figwidth'   : (8, 12), # (1-column, 2-column)
          'figheight'  : 6, #9.19
          'aspct_ratio': 3./4, #golden_mean
          'facecolor'  : "1.0",
          'edgecolor'  : "0.5"
}
### SAVING FIGURES
savefig = {'dpi'       : 300,      # figure dots per inch
           'format'    : 'png',    # png, ps, pdf, svg
           'bbox'      : 'tight',  # 'tight' or 'standard'.
           'pad_inches': 0.1      # Padding to be used when bbox is set to 'tight'
#savefig.directory   : ~        # default directory in savefig dialog box,
                                # leave empty to always use current working directory
}
