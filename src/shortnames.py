## VISIBLE STATES
# no noise
myrun_fnamesN4real = ['N4rsn8inf', 'N4run8inf']
myrun_fnamesN4binary = ['N4bsn8inf', 'N4bun8inf']
myrun_fnamesN5 = ['N5run8inf']

# vs. N, vs. matrices,
myrun_fnamesN4 = myrun_fnamesN4real + myrun_fnamesN4binary
# vs. N, vs. matrices, min = 0
myrun_fnamesN4m0 = ['N4m0rsn8inf', 'N4m0run8inf',
                    'N4m0bsn8inf', 'N4m0bun8inf']

## vs. noise

##########
## REAL ##
##########
# unstable
myrun_noiseN4realu = ['N4ruo1e-100n8100', 'N4ruo0.01n8100', 'N4ruo0.05n8100',
                      'N4ruo0.1n8100', 'N4ruo0.15n8100',
                      'N4ruo0.2n8100', 'N4ruo0.25n8100']

# min = 0
myrun_noiseN4m0realu = ['N4m0ruo1e-100n8100', 'N4m0ruo0.01n8100',
                        'N4m0ruo0.05n8100', 'N4m0ruo0.1n8100',
                        'N4m0ruo0.15n8100', 'N4m0ruo0.2n8100',
                        'N4m0ruo0.25n8100']

# devo_time = 1000
myrun_noiseN4m0realu1000 = ['N4m0ruo1e-100n81000', 'N4m0ruo0.01n81000',
                            'N4m0ruo0.05n81000', 'N4m0ruo0.1n81000',
                            'N4m0ruo0.15n81000', 'N4m0ruo0.2n81000',
                            'N4m0ruo0.25n81000']

# stable
myrun_noiseN4reals = ['N4rso1e-100n8100', 'N4rso0.01n8100', 'N4rso0.05n8100',
                      'N4rso0.1n8100', 'N4rso0.15n8100',
                      'N4rso0.2n8100', 'N4rso0.25n8100']

# min = 0
myrun_noiseN4m0reals = ['N4m0rso1e-100n8100', 'N4m0rso0.01n8100',
                        'N4m0rso0.05n8100', 'N4m0rso0.1n8100',
                        'N4m0rso0.15n8100', 'N4m0rso0.2n8100',
                        'N4m0rso0.25n8100']


##############
## UNSTABLE ##
##############
# unstable, devo_time = 100
myrun_noiseN4u = ['N4buo1e-100n8100', 'N4buo0.01n8100', 'N4buo0.05n8100',
                  'N4buo0.1n8100', 'N4buo0.15n8100',
                  'N4buo0.2n8100', 'N4buo0.25n8100']
# + noise = 0
myrun_noiseN4ux = ['N4bun8inf'] + myrun_noiseN4u

# unstable, min = 0
myrun_noiseN4m0u = ['N4m0buo1e-100n8100', 'N4m0buo0.01n8100',
                    'N4m0buo0.05n8100', 'N4m0buo0.1n8100', 'N4m0buo0.15n8100',
                    'N4m0buo0.2n8100', 'N4m0buo0.25n8100']
# + noise = 0
myrun_noiseN4m0ux = ['N4m0bun8inf'] + myrun_noiseN4m0u

# unstable, devo_time = 1000, min = 0
myrun_noiseN4m0u1000 = ['N4m0buo1e-100n81000', 'N4m0buo0.01n81000',
                        'N4m0buo0.05n81000', 'N4m0buo0.1n81000',
                        'N4m0buo0.15n81000', 'N4m0buo0.2n81000',
                        'N4m0buo0.25n81000']
# + noise = 0
myrun_noiseN4m0u1000x = ['N4m0bun8inf'] + myrun_noiseN4m0u1000

# unstable, devo_time = 10,000, min = 0
myrun_noiseN4m0u10000 = ['N4m0buo1e-100n810000', 'N4m0buo0.01n810000',
                         'N4m0buo0.05n810000', 'N4m0buo0.1n810000',
                         'N4m0buo0.15n810000', 'N4m0buo0.2n810000',
                         'N4m0buo0.25n810000']
# + noise = 0
myrun_noiseN4m0u10000x = ['N4m0bun8inf'] + myrun_noiseN4m0u10000


############
## STABLE ##
############
# stable, devo_time = 100
myrun_noiseN4s = ['N4bso1e-100n8100', 'N4bso0.01n8100', 'N4bso0.05n8100',
                  'N4bso0.1n8100', 'N4bso0.15n8100',
                  'N4bso0.2n8100', 'N4bso0.25n8100']

# + noise = 0
myrun_noiseN4sx = ['N4bsn8inf'] + myrun_noiseN4s

# stable, devo_time = 101 (NOTE: only k = 4)
myrun_noiseN4k4s101 = ['N4k4bso1e-100n8101', 'N4k4bso0.01n8101',
                       'N4k4bso0.05n8101', 'N4k4bso0.1n8101',
                       'N4k4bso0.15n8101', 'N4k4bso0.2n8101',
                       'N4k4bso0.25n8101']

# stable, devo_time = 200
myrun_noiseN4s200 = ['N4bso1e-100n8200', 'N4bso0.01n8200', 'N4bso0.05n8200',
                     'N4bso0.1n8200', 'N4bso0.15n8200',
                     'N4bso0.2n8200', 'N4bso0.25n8200']

# stable, min = 0, devo_time = 100
myrun_noiseN4m0s = ['N4m0bso1e-100n8100', 'N4m0bso0.01n8100',
                    'N4m0bso0.05n8100', 'N4m0bso0.1n8100', 'N4m0bso0.15n8100',
                    'N4m0bso0.2n8100', 'N4m0bso0.25n8100']

# + noise = 0
myrun_noiseN4m0sx = ['N4m0bsn8inf'] +  myrun_noiseN4m0s

# stable, min = 0, devo_time = 200
myrun_noiseN4m0s200 = ['N4m0bso1e-100n8200', 'N4m0bso0.01n8200',
                       'N4m0bso0.05n8200', 'N4m0bso0.1n8200',
                       'N4m0bso0.15n8200', 'N4m0bso0.2n8200',
                       'N4m0bso0.25n8200']

# stable, min = 0, devo_time = 1000
myrun_noiseN4m0s1000 = ['N4m0bso1e-100n81000', 'N4m0bso0.01n81000',
                        'N4m0bso0.05n81000', 'N4m0bso0.1n81000',
                        'N4m0bso0.15n81000', 'N4m0bso0.2n81000',
                        'N4m0bso0.25n81000']
# + noise = 0
myrun_noiseN4m0s1000x = ['N4m0bsn8inf'] + myrun_noiseN4m0s1000

# stable, min = 0, devo_time = 10000
myrun_noiseN4m0s10000 = ['N4m0bso1e-100n810000', 'N4m0bso0.01n810000',
                         'N4m0bso0.05n810000', 'N4m0bso0.1n810000',
                         'N4m0bso0.15n810000', 'N4m0bso0.2n810000',
                         'N4m0bso0.25n810000']
# + noise = 0
myrun_noiseN4m0s10000x = ['N4m0bsn8inf'] + myrun_noiseN4m0s10000


#############################
# ALTERNATIVE NOISE MODELS ##
#############################
# N = 4
# random_flip
myrun_noiseN4rf = ['N4bsrf1e-100n8100', 'N4bsrf1n8100', 'N4bsrf2n8100',
                   'N4bsrf3n8100', 'N4bsrf4n8100']

# random_flip, unstable, min = 0
myrun_noiseN4m0rfu = ['N4m0burf1e-100n8100', 'N4m0burf1n8100',
                      'N4m0burf2n8100', 'N4m0burf3n8100',
                      'N4m0burf4n8100']

# noise_function = force_different_flip
# noise_time = 'before'
myrun_noiseN4df = ['N4bsdf1e-100n8100', 'N4bsdf1n8100', 'N4bsdf2n8100',
                   'N4bsdf3n8100', 'N4bsdf4n8100']

# noise_function = neighbour_flip
# noise_time = 'before'
myrun_noiseN4nf = ['N4bsnf1e-100n8100', 'N4bsnf1n8100', 'N4bsnf2n8100',
                   'N4bsnf3n8100', 'N4bsnf4n8100']

# neighbour_flip, unstable, min = 0
myrun_noiseN4m0nfu = ['N4m0bunf1e-100n8100', 'N4m0bunf1n8100',
                      'N4m0bunf2n8100', 'N4m0bunf3n8100',
                      'N4m0bunf4n8100']

# neighbour_flip, stable, min = 0
myrun_noiseN4m0nfs = ['N4m0bsnf1e-100n8100', 'N4m0bsnf1n8100',
                      'N4m0bsnf2n8100', 'N4m0bsnf3n8100',
                      'N4m0bsnf4n8100']

# force_different_flip, unstable, min = 0
myrun_noiseN4m0dfu = ['N4m0budf1e-100n8100', 'N4m0budf1n8100',
                      'N4m0budf2n8100', 'N4m0budf3n8100',
                      'N4m0budf4n8100']

# force_different_flip, stable, min = 0
myrun_noiseN4m0dfs = ['N4m0bsdf1e-100n8100', 'N4m0bsdf1n8100',
                      'N4m0bsdf2n8100', 'N4m0bsdf3n8100',
                      'N4m0bsdf4n8100']

###############
## SHMULEVICH #
###############
# noise_time = 'shmulevich'
# neighbour_flip, unstable, min = 0
myrun_noiseN4m0nfus = ['N4m0bunf1e-100sn8100', 'N4m0bunf0.01sn8100',
                       'N4m0bunf0.05sn8100', 'N4m0bunf0.1sn8100',
                       'N4m0bunf0.15sn8100', 'N4m0bunf0.2sn8100',
                       'N4m0bunf0.25sn8100']

# neighbour_flip, stable
myrun_noiseN4nfss = ['N4bsnf1e-100sn8100', 'N4bsnf0.01sn8100',
                     'N4bsnf0.05sn8100', 'N4bsnf0.1sn8100',
                     'N4bsnf0.15sn8100', 'N4bsnf0.2sn8100',
                     'N4bsnf0.25sn8100', 'N4bsnf0.3sn8100',
                     'N4bsnf0.35sn8100', 'N4bsnf0.4sn8100',
                     'N4bsnf0.45sn8100']

# force_different_flip, stable
myrun_noiseN4dfss = ['N4bsdf1e-100sn8100', 'N4bsdf0.01sn8100',
                     'N4bsdf0.05sn8100', 'N4bsdf0.1sn8100',
                     'N4bsdf0.15sn8100', 'N4bsdf0.2sn8100',
                     'N4bsdf0.25sn8100']

# neighbour_flip, stable, min = 0
myrun_noiseN4m0nfss = ['N4m0bsnf1e-100sn8100', 'N4m0bsnf0.01sn8100',
                       'N4m0bsnf0.05sn8100', 'N4m0bsnf0.1sn8100',
                       'N4m0bsnf0.15sn8100', 'N4m0bsnf0.2sn8100',
                       'N4m0bsnf0.25sn8100']

# force_different_flip, stable, min = 0
myrun_noiseN4m0dfss = ['N4m0bsdf1e-100sn8100', 'N4m0bsdf0.01sn8100',
                       'N4m0bsdf0.05sn8100', 'N4m0bsdf0.1sn8100',
                       'N4m0bsdf0.15sn8100', 'N4m0bsdf0.2sn8100',
                       'N4m0bsdf0.25sn8100']

# noise_time = 'after'
myrun_noiseN4a = ['N4bso1e-100an8100', 'N4bso0.01an8100', 'N4bso0.05an8100',
                  'N4bso0.1an8100', 'N4bso0.15an8100',
                  'N4bso0.2an8100', 'N4bso0.25an8100']


# full_enum
myrun_noiseN4f = ['N4bso1e-100100', 'N4bso0.01100', 'N4bso0.05100',
                  'N4bso0.1100', 'N4bso0.15100',
                  'N4bso0.2100', 'N4bso0.25100']


#######################
# NON-ULTRASENSITIVE ##
#######################
## steps = 'quantiles'
# stable
myrun_noiseN4q = ['N4bso1e-100n8q100', 'N4bso0.01n8q100', 'N4bso0.05n8q100',
                  'N4bso0.1n8q100', 'N4bso0.15n8q100',
                  'N4bso0.2n8q100', 'N4bso0.25n8q100']
# unstable
myrun_noiseN4qu = ['N4buo1e-100n8q100', 'N4buo0.01n8q100', 'N4buo0.05n8q100',
                   'N4buo0.1n8q100', 'N4buo0.15n8q100',
                   'N4buo0.2n8q100', 'N4buo0.25n8q100']
# full_enum
myrun_noiseN4qf = ['N4o1e-100q100', 'N4o0.01q100', 'N4o0.05q100',
                   'N4o0.1q100', 'N4o0.15q100',
                   'N4o0.2q100', 'N4o0.25q100']

############
# DENSITY ##
############
# density = 0.5
myrun_noiseN4k2 = ['N4k2bso1e-100n8100', 'N4k2bso0.01n8100', 'N4k2bso0.05n8100',
                   'N4k2bso0.1n8100', 'N4k2bso0.15n8100',
                   'N4k2bso0.2n8100', 'N4k2bso0.25n8100']
#############
# stability #
myrun_stabilityN4 = ['N4buo1e-100n8s', 'N4buo0.01n8s', 'N4buo0.05n8s',
                     'N4buo0.1n8s', 'N4buo0.15n8s',
                     'N4buo0.2n8100s', 'N4buo0.25n8100s']

############
# OTHER N ##
############
# N < 4
myrun_noiseN2s = ['N2bso1e-100n8100', 'N2bso0.01n8100', 'N2bso0.05n8100',
                  'N2bso0.1n8100', 'N2bso0.15n8100',
                  'N2bso0.2n8100', 'N2bso0.25n8100']
myrun_noiseN3s = ['N3bso1e-100n8100', 'N3bso0.01n8100', 'N3bso0.05n8100',
                  'N3bso0.1n8100', 'N3bso0.15n8100',
                  'N3bso0.2n8100', 'N3bso0.25n8100']
# N > 4
#myrun_noiseN5u = ['N5buo1e-100n832', 'N5buo0.01n832',
#                  'N5buo0.05n832', 'N5buo0.1n832', 'N5buo0.15n832']
myrun_noiseN5s = ['N5bso1e-100n8100', 'N5bso0.01n8100',
                  'N5bso0.05n8100', 'N5bso0.1n8100', 'N5bso0.15n8100',
                  'N5bso0.2n8200', 'N5bso0.25n8200']
myrun_noiseN6s = ['N6bso1e-100n8100', 'N6bso0.01n8100',
                  'N6bso0.05n8100', 'N6bso0.1n8100', 'N6bso0.15n8100',
                  'N6bso0.2n8200', 'N6bso0.25n8200']
myrun_noiseN7s = ['N7bso1e-100n8100', 'N7bso0.01n8100',
                  'N7bso0.05n8100', 'N7bso0.1n8100', 'N7bso0.15n8100',
                  'N7bso0.2n8100', 'N7bso0.25n8100']
myrun_noiseN8s = ['N8bso1e-100n8100', 'N8bso0.01n8100',
                  'N8bso0.05n8100', 'N8bso0.1n8100', 'N8bso0.15n8100',
                  'N8bso0.2n8100', 'N8bso0.25n8100']
myrun_noiseN9s = ['N9bso1e-100n8100', 'N9bso0.01n8100',
                  'N9bso0.05n8100', 'N9bso0.1n8100', 'N9bso0.15n8100',
                  'N9bso0.2n8100', 'N9bso0.25n8100']
myrun_noiseN10s = ['N10bso1e-100n8100', 'N10bso0.01n8100',
                   'N10bso0.05n8100', 'N10bso0.1n8100', 'N10bso0.15n8100',
                   'N10bso0.2n8200', 'N10bso0.25n8200']

# meta
meta_quantiles = [myrun_noiseN4q, myrun_noiseN4qu, myrun_noiseN4qf]
#meta_shmulevich = [myrun_noiseN4dfs, myrun_noiseN4nfs]
meta_flips = [myrun_noiseN4a, myrun_noiseN4rf, myrun_noiseN4df, myrun_noiseN4nf]
meta_myrunsN4 = [myrun_noiseN4f, myrun_noiseN4u, myrun_noiseN4s]
meta_myruns_all = meta_myrunsN4 + [
    myrun_noiseN5s, myrun_noiseN6s, myrun_noiseN10s]

## EVOLUTION
# mutation rate shortnames
m4u_shortnames = ['m4r05b', 'm4r05u1', 'm4r05u001']
m3u_shortnames = ['m3r05b', 'm3r05u1', 'm3r05u001']
m2u_shortnames = ['m2r05b', 'm2r05u1', 'm2r05u001']

# selection strength shortnames
m4s_shortnames = ['m4r05b', 'm4r05s1', 'm4r05s001']

# recombination shortnames
m4r_shortnames = ['m4r05b', 'm4r0b']

# Fig 2 shortnames
six_shortnames = ['m3r05b', 'm4r05b', 'm3r05p2', 'm3r05p3', 'm2r05b', 'm4r05p0']

# Fig S2 shortnames
all_shortnames = ['m4r05b', 'm3r05b', 'm4r05p0', 'm3r05p2', 'm3r05p3',
                  'm3r05p4', 'm3r05p5', 'm3r05p6', 'm3r05p7', 'm2r05b']

# starting conditions shortnames
start_shortnames = ['b05p05q05', 'm3p05q09', 'm3p05q01', 'm3p09q09', 'm3p09q01']

# density shortnames
k2_shortnames = ['m3r05b', 'm3r05k2', 'm3r0k2', 'm3r05g', 'm3r0g']

# mut bias shortnames
mutbias_shortnames = ['b099p05q05', 'b095p05q05', 'b09p05q05', 'b085p05q05',
                      'b08p05q05', 'b075p05q05', 'b072p05q05',
                      'b07p05q05', 'b065p05q05', 'b06p05q05',
                      'b059p05q05', 'b058p05q05', 'b057p05q05', 'b056p05q05',
                      'b055p05q05', 'b054p05q05', 'b053p05q05', 'b052p05q05',
                      'b051p05q05', 'b05p05q05', 'b049p05q05',
                      'b045p05q05', 'b04p05q05', 'b035p05q05', 'b03p05q05',
                      'b025p05q05', 'b02p05q05',
                      'b015p05q05', 'b01p05q05', 'b005p05q05', 'b001p05q05']
neutbias_shortnames = ['m2b099', 'm2b095', 'm2b09', 'm2b085', 'm2b08',
                       'm2b075', 'm2b072', 'm2b07', 'm2b065', 'm2b06',
                       'm2b059', 'm2b058', 'm2b057', 'm2b056',
                       'm2b055', 'm2b054', 'm2b053', 'm2b052', 'm2b051',
                       'm2b05', 'm2b049', 'm2b045', 'm2b04', 'm2b035', 'm2b03',
                       'm2b025', 'm2b02',
                       'm2b015', 'm2b01', 'm2b005', 'm2b001']
