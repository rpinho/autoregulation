import combinations
reload(combinations)
from combinations import *
import auxiliary_functions
reload(auxiliary_functions)
from auxiliary_functions import *
import labels
reload(labels)
from labels import *

pop_filenames = {'m4r05b':   ['pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+07_period1_s0.1_u0.1-binary_r0.5'],
                 'm4r04b':    'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.4',
                 'm4r03b':    'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.3',
                 'm4r02b':    'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.2',
                 'm4r01b':    'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.1',
                 'm4r009b':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.09',
                 'm4r008b':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.08',
                 'm4r007b':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.07',
                 'm4r006b':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.06',
                 'm4r005b':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.05',
                 'm4r004b':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.04',
                 'm4r003b':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.03',
                 'm4r002b':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.02',
                 'm4r001b':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.01',
                 'm4r001k2':  'pop.evolution-N_[10]_c_[ 0.2]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.01',
                 'm4r0b':     'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary',
                 'm3r05b':   ['pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+07_period1_u0.1-binary_r0.5'],
                 'm3r04b':    'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.4',
                 'm3r03b':    'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.3',
                 'm3r02b':    'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.2',
                 'm3r01b':    'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.1',
                 'm3r005b':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.05',
                 'm3r004b':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.04',
                 'm3r003b':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.03',
                 'm3r002b':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.02',
                 'm3r001b':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.01',
                 'm3r0b':     'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary',
                 'm4r05r':   ['pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-random-n1-model4_G1.0e+05_period1_u0.1-random_r0.5',
                              'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-random-model4_G1e+05_period1_s0.1_u0.1-random_r0.5'],
                 'm3r05r':   ['pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-random-n1-model3_G1.0e+05_period1_u0.1-random_r0.5',
                              'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-random-model3_G1e+05_period1_u0.1-random_r0.5'],
                 'm3r0r':    ['pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-random-n1-model3_G1.0e+05_period1_u0.1-random',
                              'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-random-model3_G1e+05_period1_u0.1-random'],
                 'm4r05s001': 'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.01_u0.1-binary_r0.5',
                 'm4r05s1':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s1_u0.1-binary_r0.5',
                 'm4r05u1':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u1-binary_r0.5',
                 'm4r05u001': 'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.01-binary_r0.5',
                 'm3r05u001':['pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.01-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+07_period1_u0.01-binary_r0.5'],
                 'm3r05u1':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u1-binary_r0.5',
                 'm4r05g':   ['pop.evolution-N_[10]_c_[ 0.2]_min_[-1.0]_dim_[2]-noise0-graph-exp_pow-n1e+06-p0.5-csign-model4_G1e+06_period1_s0.1_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 0.2]_min_[-1.0]_dim_[2]-noise0-graph-exp_pow-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.5'],
                 'm3r05g':   ['pop.evolution-N_[10]_c_[ 0.2]_min_[-1.0]_dim_[2]-noise0-graph-exp_pow-n1e+06-p0.5-csign-model3_G1e+06_period1_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 0.2]_min_[-1.0]_dim_[2]-noise0-graph-exp_pow-p0.5-model3_G1e+06_period1_u0.1-binary_r0.5'],
                 'm4r0g':     'pop.evolution-N_[10]_c_[ 0.2]_min_[-1.0]_dim_[2]-noise0-graph-exp_pow-p0.5-model4_G1e+04_period1_s0.1_u0.1-binary',
                 'm3r0g':     'pop.evolution-N_[10]_c_[ 0.2]_min_[-1.0]_dim_[2]-noise0-graph-exp_pow-p0.5-model3_G1e+04_period1_u0.1-binary',
                 'm4r05k2':  [#'pop.evolution-N_[10]_c_[ 0.2]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+04_period1_s0.1_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 0.2]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.5'],
                 'm4r05k3':   'pop.evolution-N_[10]_c_[ 0.3]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.5',
                 'm4r05rk3':  'pop.evolution-N_[10]_c_[ 0.3]_min_[-1.0]_dim_[2]-noise0-random-model4_G1e+06_period1_s0.1_u0.1-random_r0.5',
                 'm4r05k4':   'pop.evolution-N_[10]_c_[ 0.4]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.5',
                 'm4r05k5':   'pop.evolution-N_[10]_c_[ 0.5]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.5',
                 'm4r05k6':   'pop.evolution-N_[10]_c_[ 0.6]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.5',
                 'm4r05k7':   'pop.evolution-N_[10]_c_[ 0.7]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.5',
                 'm4r05k8':   'pop.evolution-N_[10]_c_[ 0.8]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.5',
                 'm4r05k9':   'pop.evolution-N_[10]_c_[ 0.9]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary_r0.5',
                 'm4r0k2':    'pop.evolution-N_[10]_c_[ 0.2]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary',
                 'm4r0k3':    'pop.evolution-N_[10]_c_[ 0.3]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary',
                 'm4r0k4':    'pop.evolution-N_[10]_c_[ 0.4]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary',
                 'm4r0k5':    'pop.evolution-N_[10]_c_[ 0.5]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary',
                 'm4r0k6':    'pop.evolution-N_[10]_c_[ 0.6]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary',
                 'm4r0k7':    'pop.evolution-N_[10]_c_[ 0.7]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary',
                 'm4r0k8':    'pop.evolution-N_[10]_c_[ 0.8]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary',
                 'm4r0k9':    'pop.evolution-N_[10]_c_[ 0.9]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_period1_s0.1_u0.1-binary',
                 'm3r0k2':    'pop.evolution-N_[10]_c_[ 0.2]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary',
                 'm3r05k2':   'pop.evolution-N_[10]_c_[ 0.2]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.5',
                 'm3r05k3':   'pop.evolution-N_[10]_c_[ 0.3]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.5',
                 'm3r05k4':   'pop.evolution-N_[10]_c_[ 0.4]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.5',
                 'm3r05k5':   'pop.evolution-N_[10]_c_[ 0.5]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.5',
                 'm3r05k6':   'pop.evolution-N_[10]_c_[ 0.6]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.5',
                 'm3r05k7':   'pop.evolution-N_[10]_c_[ 0.7]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.5',
                 'm3r05k8':   'pop.evolution-N_[10]_c_[ 0.8]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.5',
                 'm3r05k9':   'pop.evolution-N_[10]_c_[ 0.9]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.5',
                 'm2r05k2':  ['pop.evolution-N_[10]_c_[ 0.2]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_period1_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 0.2]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_u0.1-binary_r0.5'],
                 'm2r05k3':  ['pop.evolution-N_[10]_c_[ 0.3]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_period1_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 0.3]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_u0.1-binary_r0.5'],
                 'm2r05k4':  ['pop.evolution-N_[10]_c_[ 0.4]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_period1_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 0.4]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_u0.1-binary_r0.5'],
                 'm2r05k5':  ['pop.evolution-N_[10]_c_[ 0.5]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_period1_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 0.5]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_u0.1-binary_r0.5'],
                 'm2r05k6':  ['pop.evolution-N_[10]_c_[ 0.6]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_period1_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 0.6]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_u0.1-binary_r0.5'],
                 'm2r05k7':  ['pop.evolution-N_[10]_c_[ 0.7]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_period1_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 0.7]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_u0.1-binary_r0.5'],
                 'm2r05k8':  ['pop.evolution-N_[10]_c_[ 0.8]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_period1_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 0.8]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_u0.1-binary_r0.5'],
                 'm2r05k9':  ['pop.evolution-N_[10]_c_[ 0.9]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_period1_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 0.9]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_u0.1-binary_r0.5'],
                 'm2r05b':   ['pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_period1_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_u0.1-binary_r0.5'],
                 'm2r0b':     'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_u0.1-binary',
                 'm2r05u1':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_u1-binary_r0.5',
                 'm2r05u001': 'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model2_G1e+06_u0.01-binary_r0.5',
                 'm4r05n':    'pop.evolution-N_[20]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+05_period1_s0.1_u0.1-binary_r0.5',
                 'm3r05n':    'pop.evolution-N_[20]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+05_period1_u0.1-binary_r0.5',
                 'm4r0n':     'pop.evolution-N_[20]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+05_period1_s0.1_u0.1-binary',
                 'm3r0n':     'pop.evolution-N_[20]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+05_period1_u0.1-binary',
                 'm4r05p0':  ['pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+06_s0.1_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model4_G1e+07_s0.1_u0.1-binary_r0.5'],
                 'm3r05p2':  ['pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period2_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+07_period2_u0.1-binary_r0.5'],
                 'm3r05p3':  ['pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period3_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+07_period3_u0.1-binary_r0.5'],
                 'm3r05p4':  ['pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period4_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+07_period4_u0.1-binary_r0.5'],
                 'm3r05p5':  ['pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period5_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+07_period5_u0.1-binary_r0.5'],
                 'm3r05p6':  ['pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period6_u0.1-binary_r0.5',
                              'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+07_period6_u0.1-binary_r0.5'],
                 'm3r05p7':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+07_period7_u0.1-binary_r0.5',
                 'm3p05q01':  'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.10-model3_G1e+07_period1_u0.1-binary_r0.5',
                 'm3p05q09':  'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.90-model3_G1e+07_period1_u0.1-binary_r0.5',
                 'm3p09q01':  'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.9-q0.10-model3_G1e+07_period1_u0.1-binary_r0.5',
                 'm3p09q09':  'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.9-q0.90-model3_G1e+07_period1_u0.1-binary_r0.5',
                 'b099p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias2.33_r0.5',
                 'b095p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias1.64_r0.5',
                 'b09p05q05': 'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias1.28_r0.5',
                 'b085p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias1.037_r0.5',
                 'b08p05q05': 'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias0.84_r0.5',
                 'b075p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias0.674_r0.5',
                 'b072p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias0.583_r0.5',
                 'b07p05q05': 'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias0.52_r0.5',
                 'b065p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias0.385_r0.5',
                 'b06p05q05': 'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias0.25_r0.5',
                 'b059p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias0.228_r0.5',
                 'b058p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias0.2_r0.5',
                 'b057p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias0.176_r0.5',
                 'b056p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias0.151_r0.5',
                 'b055p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias0.126_r0.5',
                 'b054p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias0.1_r0.5',
                 'b053p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias0.076_r0.5',
                 'b052p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias0.05_r0.5',
                 'b051p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias0.025_r0.5',
                 'b05p05q05': 'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-binary_r0.5',
                 'b049p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias-0.025_r0.5',
                 'b045p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias-0.126_r0.5',
                 'b04p05q05': 'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias-0.25_r0.5',
                 'b035p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias-0.385_r0.5',
                 'b03p05q05': 'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias-0.52_r0.5',
                 'b025p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias-0.674_r0.5',
                 'b02p05q05': 'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias-0.84_r0.5',
                 'b015p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias-1.037_r0.5',
                 'b01p05q05': 'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias-1.28_r0.5',
                 'b005p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias-1.64_r0.5',
                 'b001p05q05':'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model3_G1e+07_period1_u0.1-bias-2.33_r0.5',
                 'b08p05q09': 'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.90-model3_G1e+07_period1_u0.1-bias0.84_r0.5',
                 'b08p05q01': 'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.10-model3_G1e+07_period1_u0.1-bias0.84_r0.5',
                 'b02p05q09': 'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.90-model3_G1e+07_period1_u0.1-bias-0.84_r0.5',
                 'b02p05q01': 'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.10-model3_G1e+07_period1_u0.1-bias-0.84_r0.5',
                 'm2b099':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias2.33_r0.5',
                 'm2b095':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias1.64_r0.5',
                 'm2b09':     'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias1.28_r0.5',
                 'm2b085':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias1.037_r0.5',
                 'm2b08':     'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias0.84_r0.5',
                 'm2b075':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias0.674_r0.5',
                 'm2b072':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias0.583_r0.5',
                 'm2b07':     'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias0.52_r0.5',
                 'm2b065':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias0.385_r0.5',
                 'm2b06':     'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias0.25_r0.5',
                 'm2b059':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias0.228_r0.5',
                 'm2b058':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias0.2_r0.5',
                 'm2b057':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias0.176_r0.5',
                 'm2b056':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias0.151_r0.5',
                 'm2b055':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias0.126_r0.5',
                 'm2b054':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias0.1_r0.5',
                 'm2b053':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias0.076_r0.5',
                 'm2b052':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias0.05_r0.5',
                 'm2b051':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias0.025_r0.5',
                 'm2b05':     'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-binary_r0.5',
                 'm2b049':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias-0.025_r0.5',
                 'm2b045':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias-0.126_r0.5',
                 'm2b04':     'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias-0.25_r0.5',
                 'm2b035':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias-0.385_r0.5',
                 'm2b03':     'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias-0.52_r0.5',
                 'm2b025':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias-0.674_r0.5',
                 'm2b02':     'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias-0.84_r0.5',
                 'm2b015':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias-1.037_r0.5',
                 'm2b01':     'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias-1.28_r0.5',
                 'm2b005':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias-1.64_r0.5',
                 'm2b001':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-model2_G1e+07_u0.1-bias-2.33_r0.5',
                 'm3r05d':    'qp.pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise0-binary-p0.5-q0.50-normal0.145-model3_G1e+06_period1_u0.1-binary_r0.5',
                 'm3r05m':    'pop.evolution-N_[10]_c_[ 1.]_min_[0.0]_dim_[2]-noise0-binary-p0.5-model3_G1e+06_period1_u0.1-binary_r0.5',
                 'm3r05a1':   'pop.evolution-N_[10]_c_[ 1.]_min_[-1.0]_dim_[2]-noise1e-100-binary-p0.5-cSigmoid2_a1.00T10E1e-04M100-model3_G1e+06_period1_u0.1-binary_r0.5',
                 'm3r05a1m':  'pop.evolution-N_[10]_c_[ 1.]_min_[0.0]_dim_[2]-noise1e-100-binary-p0.5-cSigmoid_a1.00T10E1e-04M100-model3_G1e+06_period1_u0.1-binary_r0.5',
                 'm4r05a1m':  'pop.evolution-N_[10]_c_[ 1.]_min_[0.0]_dim_[2]-noise1e-100-binary-p0.5-cSigmoid_a1.00T10E1e-04M100-model4_G1e+06_period1_s0.1_u0.1-binary_r0.5',
                 'm3r05a10m': 'pop.evolution-N_[10]_c_[ 1.]_min_[0.0]_dim_[2]-noise1e-100-binary-p0.5-cSigmoid_a10.00T10E1e-04M100-model3_G1e+06_period1_u0.1-binary_r0.5'
                 }

delete_x_index = [42, 43, 45, 46, 48, 49, 52, 58, 59, 61, 62, 64, 65, 67, 74, 75, 77, 78, 79, 81, 83]#None

def pop_model(shortname):
    """model, period, generations, sel_strength, mut_rate, rec_rate, bit, density, degree, n_runs, precision, _filter, suffix, binary, graph, rec_function, p, q, mut_bias, min, devo_time, gpmap, a"""
    if shortname == 'm4r05n':
        return 4, 1, 1e5,   s, u,  r, 20,  c, 20, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r0n':
        return 4, 1, 1e5,   s, u,  0, 20,  c, 20, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r05r':
        return 4, 1, 1e5,   s, u,  r,  N,  c,  N, n_runs,  True,           None,                '', False, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r05rk3':
        return 4, 1, 1e6,   s, u,  r,  N, .3,  3, n_runs,  True,           None,                '', False, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r05g':
        return 4, 1, 1e6,   s, u,  r,  N, .2,  2, n_runs,  True,           None,                '',  True,  True, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r0g':
        return 4, 1, 1e4,   s, u,  0,  N, .2,  2, n_runs,  True,           None,                '',  True,  True, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r05p0':
        return 4, 0, 1e7,   s, u,  r,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r05s1':
        return 4, 1, 1e6,   1, u,  r,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r05s001':
        return 4, 1, 1e6, .01, u,  r,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r05u1':
        return 4, 1, 1e6,   s, 1,  r,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r05u001':
        return 4, 1, 1e6,   s,.01, r,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r05a1m':
        return 4, 1, 1e6,   s, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None, None, None, None,  0., 100, 'cSigmoid',  1
    if shortname == 'm4r05b':
        return 4, 1, 1e7,   s, u,  r,  N,  c,  N,    200,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r04b':
        return 4, 1, 1e6,   s, u, .4,  N,  c,  N, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r03b':
        return 4, 1, 1e6,   s, u, .3,  N,  c,  N, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r02b':
        return 4, 1, 1e6,   s, u, .2,  N,  c,  N, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r01b':
        return 4, 1, 1e6,   s, u, .1,  N,  c,  N, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r009b':
        return 4, 1, 1e6,   s, u, .09, N,  c,  N,     30,  True,           None, '-flat_recombine',  True, False, 'flat_recombine', None, None, None, -1., inf, None, None
    if shortname == 'm4r008b':
        return 4, 1, 1e6,   s, u, .08, N,  c,  N,     30,  True,           None, '-flat_recombine',  True, False, 'flat_recombine', None, None, None, -1., inf, None, None
    if shortname == 'm4r007b':
        return 4, 1, 1e6,   s, u, .07, N,  c,  N,     30,  True,           None, '-flat_recombine',  True, False, 'flat_recombine', None, None, None, -1., inf, None, None
    if shortname == 'm4r006b':
        return 4, 1, 1e6,   s, u, .06, N,  c,  N,     30,  True,           None, '-flat_recombine',  True, False, 'flat_recombine', None, None, None, -1., inf, None, None
    if shortname == 'm4r005b':
        return 4, 1, 1e6,   s, u, .05, N,  c,  N,     30,  True,           None, '-flat_recombine',  True, False, 'flat_recombine', None, None, None, -1., inf, None, None
    if shortname == 'm4r004b':
        return 4, 1, 1e6,   s, u, .04, N,  c,  N,     30,  True,           None, '-flat_recombine',  True, False, 'flat_recombine', None, None, None, -1., inf, None, None
    if shortname == 'm4r003b':
        return 4, 1, 1e6,   s, u, .03, N,  c,  N,     30,  True,           None, '-flat_recombine',  True, False, 'flat_recombine', None, None, None, -1., inf, None, None
    if shortname == 'm4r002b':
        return 4, 1, 1e6,   s, u, .02, N,  c,  N,     30,  True,           None, '-flat_recombine',  True, False, 'flat_recombine', None, None, None, -1., inf, None, None
    if shortname == 'm4r001b':
        return 4, 1, 1e6,   s, u, .01, N,  c,  N,     30,  True,           None, '-flat_recombine',  True, False, 'flat_recombine', None, None, None, -1., inf, None, None
    if shortname == 'm4r0b':
        return 4, 1, 1e6,   s, u,  0,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r001k2':
        return 4, 1, 1e6,   s, u, .01, N, .2,  2,     30,  True,           None, '-flat_recombine',  True, False, 'flat_recombine', None, None, None, -1., inf, None, None
    if shortname == 'm4r05k2':
        return 4, 1, 1e6,   s, u,  r,  N, .2,  2, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r05k3':
        return 4, 1, 1e6,   s, u,  r,  N, .3,  3,     40,  True,           None,           '-cfix',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r05k4':
        return 4, 1, 1e6,   s, u,  r,  N, .4,  4,     40,  True,           None,           '-cfix',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r05k5':
        return 4, 1, 1e6,   s, u,  r,  N, .5,  5,     40,  True,           None,           '-cfix',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r05k6':
        return 4, 1, 1e6,   s, u,  r,  N, .6,  6,     40,  True,           None,           '-cfix',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r05k7':
        return 4, 1, 1e6,   s, u,  r,  N, .7,  7,     40,  True,           None,           '-cfix',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r05k8':
        return 4, 1, 1e6,   s, u,  r,  N, .8,  8,     40,  True,           None,           '-cfix',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r05k9':
        return 4, 1, 1e6,   s, u,  r,  N, .9,  9,     40,  True,           None,           '-cfix',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r0k2':
        return 4, 1, 1e6,   s, u,  0,  N, .2,  2, n_runs,  True,           None,           '-cfix',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r0k3':
        return 4, 1, 1e6,   s, u,  0,  N, .3,  3, n_runs,  True,           None,           '-cfix',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r0k4':
        return 4, 1, 1e6,   s, u,  0,  N, .4,  4, n_runs,  True,           None,           '-cfix',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r0k5':
        return 4, 1, 1e6,   s, u,  0,  N, .5,  5, n_runs,  True,           None,           '-cfix',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r0k6':
        return 4, 1, 1e6,   s, u,  0,  N, .6,  6, n_runs,  True,           None,           '-cfix',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r0k7':
        return 4, 1, 1e6,   s, u,  0,  N, .7,  7, n_runs,  True,           None,           '-cfix',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r0k8':
        return 4, 1, 1e6,   s, u,  0,  N, .8,  8, n_runs,  True,           None,           '-cfix',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm4r0k9':
        return 4, 1, 1e6,   s, u,  0,  N, .9,  9, n_runs,  True,           None,           '-cfix',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05n':
        return 3, 1, 1e5, inf, u,  r, 20,  c, 20, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r0n':
        return 3, 1, 1e5, inf, u,  0, 20,  c, 20, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05m':
        return 3, 1, 1e6, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None, None, None, None,  0., inf, None, None
    if shortname == 'm3r05a1':
        return 3, 1, 1e6, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None, None, None, None, -1., 100, 'cSigmoid2', 1
    if shortname == 'm3r05a1m':
        return 3, 1, 1e6, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None, None, None, None,  0., 100, 'cSigmoid',  1
    if shortname == 'm3r05a10m':
        return 3, 1, 1e6, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None, None, None, None,  0., 100, 'cSigmoid',  10
    if shortname == 'm3r05r':
        return 3, 1, 1e5, inf, u,  r,  N,  c,  N, n_runs,  True,           None,                '', False, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r0r':
        return 3, 1, 1e5, inf, u,  0,  N,  c,  N, n_runs,  True,           None,                '', False, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05g':
        return 3, 1, 1e6, inf, u,  r,  N, .2,  2, n_runs,  True,           None,                '',  True,  True, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r0g':
        return 3, 1, 1e4, inf, u,  0,  N, .2,  2, n_runs,  True,           None,                '',  True,  True, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05p2':
        return 3, 2, 1e7, inf, u,  r,  N,  c,  N,    200,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05p3':
        return 3, 3, 1e7, inf, u,  r,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05p4':
        return 3, 4, 1e7, inf, u,  r,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05p5':
        return 3, 5, 1e7, inf, u,  r,  N,  c,  N,    200,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05p6':
        return 3, 6, 1e7, inf, u,  r,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05p7':
        return 3, 7, 1e7, inf, u,  r,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05u1':
        return 3, 1, 1e6, inf, 1,  r,  N,  c,  N, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05u001':
        return 3, 1, 1e7, inf,.01, r,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3p09q01':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3p05q01':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if 'm3p' in shortname:
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'b099p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, 2.33, -1., inf, None, None
    if shortname == 'b095p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, 1.64, -1., inf, None, None
    if shortname == 'b09p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, 1.28, -1., inf, None, None
    if shortname == 'b085p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,1.037, -1., inf, None, None
    if shortname == 'b08p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,  .84, -1., inf, None, None
    if shortname == 'b075p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .674, -1., inf, None, None
    if shortname == 'b072p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .583, -1., inf, None, None
    if shortname == 'b07p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,  .52, -1., inf, None, None
    if shortname == 'b065p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .385, -1., inf, None, None
    if shortname == 'b06p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,  .25, -1., inf, None, None
    if shortname == 'b059p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .228, -1., inf, None, None
    if shortname == 'b058p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,   .2, -1., inf, None, None
    if shortname == 'b057p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .176, -1., inf, None, None
    if shortname == 'b056p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .151, -1., inf, None, None
    if shortname == 'b055p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .126, -1., inf, None, None
    if shortname == 'b054p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,   .1, -1., inf, None, None
    if shortname == 'b053p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .076, -1., inf, None, None
    if shortname == 'b052p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,  .05, -1., inf, None, None
    if shortname == 'b051p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .025, -1., inf, None, None
    if shortname == 'b05p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,    0, -1., inf, None, None
    if shortname == 'b049p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,-.025, -1., inf, None, None
    if shortname == 'b045p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,-.126, -1., inf, None, None
    if shortname == 'b04p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, -.25, -1., inf, None, None
    if shortname == 'b035p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,-.385, -1., inf, None, None
    if shortname == 'b03p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, -.52, -1., inf, None, None
    if shortname == 'b025p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,-.674, -1., inf, None, None
    if shortname == 'b02p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, -.84, -1., inf, None, None
    if shortname == 'b015p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,-1.037, -1., inf, None, None
    if shortname == 'b01p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,-1.28, -1., inf, None, None
    if shortname == 'b005p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,-1.64, -1., inf, None, None
    if shortname == 'b001p05q05':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,-2.33, -1., inf, None, None
    if shortname == 'b08p05q09':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .9,  .84, -1., inf, None, None
    if shortname == 'b08p05q01':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,    100,  True,           None,                '',  True, False, None,   .5,   .1,  .84, -1., inf, None, None
    if shortname == 'b02p05q09':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .9, -.84, -1., inf, None, None
    if shortname == 'b02p05q01':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,    100,  True,           None,                '',  True, False, None,   .5,   .1, -.84, -1., inf, None, None
    if shortname == 'm2b099':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, 2.33, -1., inf, None, None
    if shortname == 'm2b095':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, 1.64, -1., inf, None, None
    if shortname == 'm2b09':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, 1.28, -1., inf, None, None
    if shortname == 'm2b085':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,1.037, -1., inf, None, None
    if shortname == 'm2b08':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,  .84, -1., inf, None, None
    if shortname == 'm2b075':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .674, -1., inf, None, None
    if shortname == 'm2b072':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .583, -1., inf, None, None
    if shortname == 'm2b07':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,  .52, -1., inf, None, None
    if shortname == 'm2b065':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .385, -1., inf, None, None
    if shortname == 'm2b06':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,  .25, -1., inf, None, None
    if shortname == 'm2b059':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .228, -1., inf, None, None
    if shortname == 'm2b058':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,   .2, -1., inf, None, None
    if shortname == 'm2b057':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .176, -1., inf, None, None
    if shortname == 'm2b056':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .151, -1., inf, None, None
    if shortname == 'm2b055':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .126, -1., inf, None, None
    if shortname == 'm2b054':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,   .1, -1., inf, None, None
    if shortname == 'm2b053':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .076, -1., inf, None, None
    if shortname == 'm2b052':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,  .05, -1., inf, None, None
    if shortname == 'm2b051':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, .025, -1., inf, None, None
    if shortname == 'm2b05':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,    0, -1., inf, None, None
    if shortname == 'm2b049':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,-.025, -1., inf, None, None
    if shortname == 'm2b045':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,-.126, -1., inf, None, None
    if shortname == 'm2b04':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, -.25, -1., inf, None, None
    if shortname == 'm2b035':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,-.385, -1., inf, None, None
    if shortname == 'm2b03':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, -.52, -1., inf, None, None
    if shortname == 'm2b025':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,-.674, -1., inf, None, None
    if shortname == 'm2b02':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5, -.84, -1., inf, None, None
    if shortname == 'm2b015':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,-1.037, -1., inf, None, None
    if shortname == 'm2b01':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,-1.28, -1., inf, None, None
    if shortname == 'm2b005':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,-1.64, -1., inf, None, None
    if shortname == 'm2b001':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,     50,  True,           None,                '',  True, False, None,   .5,   .5,-2.33, -1., inf, None, None
    if shortname == 'm3r05d':
        return 3, 1, 1e6, inf, u,  r,  N,  c,  N,    100,  True,           None,                '',  True, False, None,   .5,   .5, None, -1., inf, None, None
    if shortname == 'm3r05b':
        return 3, 1, 1e7, inf, u,  r,  N,  c,  N,    300,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r04b':
        return 3, 1, 1e6, inf, u, .4,  N,  c,  N,     30,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r03b':
        return 3, 1, 1e6, inf, u, .3,  N,  c,  N,     30,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r02b':
        return 3, 1, 1e6, inf, u, .2,  N,  c,  N,     30,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r01b':
        return 3, 1, 1e6, inf, u, .1,  N,  c,  N,     30,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r005b':
        return 3, 1, 1e6, inf, u, .05, N,  c,  N,     30,  True,           None, '-flat_recombine',  True, False, 'flat_recombine', None, None, None, -1., inf, None, None
    if shortname == 'm3r004b':
        return 3, 1, 1e6, inf, u, .04, N,  c,  N,     30,  True,           None, '-flat_recombine',  True, False, 'flat_recombine', None, None, None, -1., inf, None, None
    if shortname == 'm3r003b':
        return 3, 1, 1e6, inf, u, .03, N,  c,  N,     30,  True,           None, '-flat_recombine',  True, False, 'flat_recombine', None, None, None, -1., inf, None, None
    if shortname == 'm3r002b':
        return 3, 1, 1e6, inf, u, .02, N,  c,  N,     30,  True,           None, '-flat_recombine',  True, False, 'flat_recombine', None, None, None, -1., inf, None, None
    if shortname == 'm3r001b':
        return 3, 1, 1e6, inf, u, .01, N,  c,  N,     30,  True,           None, '-flat_recombine',  True, False, 'flat_recombine', None, None, None, -1., inf, None, None
    if shortname == 'm3r0b':
        return 3, 1, 1e6, inf, u,  0,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r0k2':
        return 3, 1, 1e6, inf, u,  0,  N, .2,  2, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05k2':
        return 3, 1, 1e6, inf, u,  r,  N, .2,  2, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05k3':
        return 3, 1, 1e6, inf, u,  r,  N, .3,  3,     44,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05k4':
        return 3, 1, 1e6, inf, u,  r,  N, .4,  4,     44,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05k5':
        return 3, 1, 1e6, inf, u,  r,  N, .5,  5,     44,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05k6':
        return 3, 1, 1e6, inf, u,  r,  N, .6,  6,     44,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05k7':
        return 3, 1, 1e6, inf, u,  r,  N, .7,  7,     44,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05k8':
        return 3, 1, 1e6, inf, u,  r,  N, .8,  8,     44,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm3r05k9':
        return 3, 1, 1e6, inf, u,  r,  N, .9,  9,     44,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm2r05b':
        return 2, 0, 1e7, inf, u,  r,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm2r0b':
        return 2, 0, 1e6, inf, u,  0,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm2r05u1':
        return 2, 0, 1e6, inf, 1,  r,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm2r05u001':
        return 2, 0, 1e6, inf,.01, r,  N,  c,  N,    100,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm2r05k2':
        return 2, 0, 1e6, inf, u,  r,  N, .2,  2, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm2r05k3':
        return 2, 0, 1e6, inf, u,  r,  N, .3,  3, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm2r05k4':
        return 2, 0, 1e6, inf, u,  r,  N, .4,  4, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm2r05k5':
        return 2, 0, 1e6, inf, u,  r,  N, .5,  5, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm2r05k6':
        return 2, 0, 1e6, inf, u,  r,  N, .6,  6, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm2r05k7':
        return 2, 0, 1e6, inf, u,  r,  N, .7,  7, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm2r05k8':
        return 2, 0, 1e6, inf, u,  r,  N, .8,  8, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None
    if shortname == 'm2r05k9':
        return 2, 0, 1e6, inf, u,  r,  N, .9,  9, n_runs,  True,           None,                '',  True, False, None, None, None, None, -1., inf, None, None

def get_pop_datapoints(
        shortname=None, generations=1e6, _filter=None, n=100, start=-0.1):
#    if shortname == 'm4r05k2':
#        return append(get_pop_datapoints(generations=1e4),
#                      get_pop_datapoints(generations=pop_model(shortname)[2],
#                                         _filter = delete_x_index)[51:])
    # shortname == 'm3r05p2', 'm3r05p3', 'm3r05p5':
    if shortname == 'm3r05d' and n == 'all1000':
        return append(range(1000),
                      get_pop_datapoints(
                          generations=pop_model(shortname)[2])[42:])
    elif shortname and generations > 1e6:
        return append(get_pop_datapoints(generations=1e6),
                      get_pop_datapoints(
                          generations=pop_model(shortname)[2])[79:])
    else:
        return get_int_datapoints(start, log10(generations), n, _filter)

def get_runs(shortname, n_runs=None):
    if not n_runs:
        n_runs = pop_model(shortname)[9]
    return xrange(1, n_runs+1)

# shortname == 'm3r05d' does not use parameter generations - see above
def get_generations(shortname, n=None, generations=1e7):
    return iter(get_pop_datapoints(shortname, generations)[:n])

def get_params(shortname, stat):
    params = array(pop_model(shortname))
    if stat == 'pop.evolution.u0':
        pass
        #params[12] = '-u0g10'
    if shortname == 'm3r0b':# and (stat == 's' or stat == 'f'):
        pass
        #params[9]  = 50
        #params[11] = delete_x_index
    if shortname == 'b08p05q01' and (stat == '1-p' or stat == 'p'):
        pass
        #params[2] = 1e6
        #params[9] = 50
    if shortname == 'm2r05b' and (stat == 'a.matrix'):# or stat == '1-p'):
        params[2] = 1e6
    return params

def get_experiment_extension(experiment, shortname=''):
    if 'pop.evolution' in experiment:
        return '.dat'
    else:
        return '.npy'

def get_experiment_dir(experiment, shortname=''):
    if 'pop.evolution' in experiment:
        return pop_dir
    if experiment == 'pop.perturb_and_analyse':
        return pa_dir
    if experiment == 'pop.epistasis':
        return epi_dir
    if experiment == 'pop.robustness':
        return pa_dir
    if experiment == 'pop.survivability':
        return surv_dir
    if 'stable.pop' in experiment:
        return pop_dir + 'stable/'
    if 'qp.pop' in experiment:
        return pop_dir + 'qp/'
    if 'autoregulation' in experiment:
        return data_dir + 'autoregulation/'
    if experiment == 'qp' and 'm3r05' in shortname:
        return data_dir + 'qp/' + shortname + '/'
    if experiment == 'phenotypes':
        return data_dir + 'phenotypes/'
    if experiment == 'stability':
        return data_dir + 'stability/'
    else:
        return data_dir

def get_pop_dir(shortname, run, _dir=pop_dir, n_runs=300, step=10):
    if _dir == data_dir:
        return _dir
    a = arange(1,  n_runs+1, step)
    b = arange(10, n_runs+1, step)
    i = (run-1)/step
    return _dir + shortname + '/n%d-%d/' %(a[i], b[i])

def _get_figs_save_dir(journal, experiment='', shortname=''):
    if journal:
        return report_dir + experiment + '/figs/'
    if 'qp' in experiment:
        return figs_dir   + experiment.split('.')[0] + '/' + shortname + '/'
    if 'transition' in experiment:
        return figs_dir   + experiment + '/'
    else:
        return figs_dir

def get_figs_save_dir(journal, experiment='', shortname=''):
    dir = _get_figs_save_dir(journal, experiment, shortname)
    ensure_dir_exists(dir)
    return dir

def get_noise_figname(fname):
    return fname.replace('o1e-100','').replace('1e-100','')

def get_figname(
        _stats, shortnames,
        stat_label=get_stat_label, shortname_label=get_model_label,
        diversity={}, errors={}, suffix2=''):

    fname = str(_stats) if len(_stats) > 1 else stat_label(_stats[0])
    if 'diversity' in _stats:
        fname = fname.replace('diversity', str(diversity.keys()))
    suffix = str(shortnames) if len(shortnames) > 1 else shortnames[0]
    if errors.values():
        suffix += '-errorbars%s' %errors.values()
    return '-'.join((fname, suffix, suffix2))

def load_qpgsr_matrices(shortname):
    fnames = ['qpg', 'qps', 'qpr']
    masks  = [False, True, True]
    return [load_masked_npy('.matrix-'.join((fname, shortname)), mask)
            for fname, mask in zip(fnames, masks)]

def get_non_evolved_pop(
        experiment, stable, run, act_fraction=p, p=None, q=None, period=1,
        pop_size=1000, binary=True, genotype_func='diagonal_p_genotype',
        extension='.dat', verbose=True):

    filename = get_experiment_dir(experiment) + experiment + get_filename()
    filename += get_suffix(
        stable=stable, binary=binary, act_fraction=act_fraction, p=p, q=q,
        period=period, genotype_func=genotype_func, pop_size=pop_size,
        run=run, experiment=experiment)

    _load = load if extension == '.npy' else cPickle.load

    return load_file(filename + extension, verbose, _load), filename

def get_pop(
        shortname, run=None, generations=None, suffix='',
        extension=None, verbose=True, experiment='pop.evolution',
        _dir=None, dry_run=False, G_string=1000):

    ## default behaviors
    if not run:
        run = get_params(shortname, experiment)[9]
    if generations == None:
        generations = get_params(shortname, experiment)[2]
    if not _dir:
        _dir = get_experiment_dir(experiment)
    if not extension:
        extension = get_experiment_extension(experiment)

    filenames = atleast_1d(pop_filenames[shortname])
    # for stats, start with the most recent filename/params
    if experiment != 'pop.evolution': filenames = filenames[::-1]

    for filename in filenames:
        if shortname not in _dir:
            _dir = get_pop_dir(shortname, run, _dir)
        filename = _dir + filename.replace('pop.evolution', experiment)
        filename += print_run(run)
        fname1 = filename + print_generations(generations, precision=True)
        fname1 += suffix
        try:
            file = open(fname1 + extension)
        except IOError:
            if verbose:
                print 'cannot open', fname1 + extension
            try:
                file = open(fname1) # try no extension
            except IOError:
                if verbose:
                    print 'cannot open', fname1
#                if generations > G_string and shortname == 'm3r0b':
#                    fname2 = filename + print_generations(
#                        generations, precision=False) + suffix
#                    try:
#                        file = open(fname2 + extension)
#                    except IOError:
#                        if verbose: print 'cannot open', fname2 + extension
#                    else:
#                        return handle_insecure_string_pickle(
#                            fname2, file, extension, verbose, dry_run)
            else:
                return handle_insecure_string_pickle(
                    fname1, file, '', verbose, dry_run)
        else:
            return handle_insecure_string_pickle(
                fname1, file, extension, verbose, dry_run)
    return 0, 0


def ensure_dir_exists(dir):
    if (not os.path.exists(dir)):
        os.makedirs(dir)

def handle_insecure_string_pickle(
        filename, file, extension='.dat', verbose=True, dry_run=False):
    if dry_run:
        return filename + extension
    try:
        data = load(file) if extension == '.npy' else cPickle.load(file)
    except ValueError or EOFError:
        if verbose:
            print 'cannot open', filename + extension
    else:
        if verbose:
            print 'opening', filename + extension
        return data, filename
    return 0,0

#def rename_and_rm_files(shortname, runs)

def mv_files(
        shortname, start, stop, datapoints, step, dry_run=True, verbose=True):
    fnames = []
    for run in arange(start, stop + 1):
        for generation in datapoints:
            filename = atleast_1d(
                get_pop(
                    shortname, run, generation, verbose=False, dry_run=True))[-1]
            if filename:
                new = filename.replace(
                    get_pop_dir(shortname, run),
                    get_pop_dir(shortname, run + step)).replace(
                        print_run(run), print_run(run + step))
                fnames.append((filename, new))
                if not dry_run: mv_file(filename, new, verbose)
    return fnames

def rm_files_consecutive(
        shortname, start, stop, datapoints, dry_run=True, verbose=True,
        get_pop_verbose=False):
    fnames = []
    for run in arange(start, stop + 1):
        # g1 > g2 because list is in reverse
        for g1, g2 in pairwise(datapoints[::-1]):
            file1 = atleast_1d(
                get_pop(
                    shortname, run, g1, verbose=get_pop_verbose, dry_run=True))[-1]
            file2 = atleast_1d(
                get_pop(
                    shortname, run, g2, verbose=get_pop_verbose, dry_run=True))[-1]
            if file1 and not file2:
                fnames.append(file1)
                if not dry_run: rm_file(file1, verbose)
    return fnames

# used ctime on corn and mtime on evo
def rm_files_time(
        shortname, start, stop, datapoints, _time,
        time_function=os.path.getctime, dry_run=True, verbose=True):
    fnames = []
    for run in arange(start, stop + 1):
        for generation in datapoints:
            filename = atleast_1d(
                get_pop(
                    shortname, run, generation, verbose=False, dry_run=True))[-1]
            if filename and time_function(filename) < _time:
                fnames.append(filename)
                if not dry_run: rm_file(filename, verbose)
    return array(fnames)

# [rm_file(file, True) for file in rm_files_time()]
def rm_file(filename, verbose=False):
    if verbose:
        print 'removing', filename
    os.remove(filename)

def mv_file(old, new, verbose=False):
    if verbose:
        print 'moving', old, 'to', new
    os.renames(old, new)

def load_masked_npy(fname, mask=False):
    data = ma.masked_invalid(load(data_dir + fname + '.npy'))
    if mask:
        data.mask = load(data_dir + fname + '.mask.npy')
    return data

def load_file(filename, verbose=False, _load=cPickle.load, **args):
    try:
        file = open(filename)
    except IOError:
        print 'cannot open', filename
        return {}
    else:
        if verbose:
            print 'opening', filename
        data = _load(file)
        file.close()
        return data

# for missing files
def load_nan_obj(filename, size, attr, verbose=False):
    class Empty:
        pass
    try:
        file = open(filename)
    except IOError:
        empty = Empty()
        setattr(empty, attr, array([NaN]*size))
        return empty
    else:
        if verbose:
            print 'opening', filename
        data = cPickle.load(file)
        file.close()
        return data

# stability
def load_last_line_of_txt(filename, n=None, verbose=False, dtype=float):
    try:
        data = loadtxt(filename, dtype)
    except IOError:
        if verbose:
            print 'cannot open', filename
        return 0, 0
#    except UserWarning:
#        if verbose:
#            print 'Empty input file:', filename
#        return 0, 0
    else:
        # load at closest point to n_samples == n
        if n:
            i = find_nearest_i(data.T[0], n)
        # load last line by default
        else:
            i = -1
        try:
            line = data[i]
        except IndexError:
            if verbose:
                print 'Empty input file:', filename
            return 0, 0
        else:
            if verbose:
                print 'opening', filename
            return line

def load_stability_from_txt_file(filename, n=None, verbose=True, *args):
    samples, n = load_last_line_of_txt(filename, n, verbose)
    if samples:
        return n / samples
    else:
        return 0

# always call with samples = 'n1e+08'
def load_stable_states_from_txt_file(
        filename, n=None, samples=[1e8, 1e6], devo_times=[16, 32, 100, 200],
        verbose=False):

    # get all filenames with parameter substitution
    filenames = get_replaced_suffix(filename, samples, devo_times)
    # load all runs
    data = [load_last_line_of_txt(filename, n, verbose)
            for filename in filenames]

    # data transformation
    data = array(data)
    # get x
    n_samples = data.T[0]
    # get y
    n_states = data.T[1]

    # get position of max visible states
    imax = n_states.argmax()

    # get number of runs
    n_runs = where(n_samples == 1)[0].size

    if verbose:
        if n_runs > 1:
            print 'ATTENTION: using ', n_runs, ' runs!'
        if devo_times:
            print 'T = ', devo_times
        if samples:
            print 'samples = ', samples
        print 'samples = %g\t\tn = %d' %(n_samples[imax], n_states[imax])

    return n_states[imax]

def save_file(filename, data, flag=True):
    if flag:
        try:
            file = open(filename)
        except IOError:
            pass
        else:
            print 'rewrote', filename

    file = open(filename, 'wb')
    cPickle.dump(data, file)
    file.close()

def format_output(data, padding, converter=repr):
    if len(padding) < len(data):
        return 'padding has to be at least as big as data'
    else:
        output = [converter(x).rjust(x_max)
                  for x, x_max in zip(data, padding[:len(data)])] + ['\n']
        return ' '.join(output)

def save_txt(filename, data, padding, mode='ab'):
    file = open(filename, mode)
    file.write(format_output(data, padding))
    file.close()

def save_txt2(filename, data, mode='ab'):
    file = open(filename, mode)
    x, y = data
    savetxt(file, array([[x],[y]]).T, '%d\t%d')
    file.close()

def save_txt3(filename, data, mode='ab', header=''):
    file = open(filename, mode)
    try:
        savetxt(file, array([data]), '%d', header=header) #newline = ', ')
    except TypeError:
        # numpy.__version__ < 1.7.0
        file.write('# %s\n'%header)
        savetxt(file, array([data]), '%d')
    file.close()

def save_phenotype_sample(
        filename, unique_phenotypes, phenotypes, samples, only_unique):

    save_txt2(filename, (samples, unique_phenotypes.size))

    if only_unique:
        data = unique_phenotypes
    else:
        data = phenotypes
        # Note: phenotypes are decimals ONLY!
    save_txt3(filename.replace('.tsv', '.txt'), data, 'w', str(samples))

def split_filename_from_workdir(filename):
    return filename.split('model/')[-1]

# get all filename variations
def get_replaced_filename(
        filename='', myrun=None,
        bits=None, min_=None, dim=None, noise=None,
        samples=None, devo_time=None,
        default_samples=1e8, default_devo_time=100):

    suffix = ''

    # myrun is fname
    if not filename and myrun:
        if isinstance(myrun, str):
            myrun = cPickle.load(open(logs_dir + myrun + '.run'))
        filename = split_filename_from_workdir(myrun.filename)

    # myrun is object
    if bits and myrun:
        old = print_N(myrun.bits)
        new = print_N(bits)

    elif min_ is not None and myrun:
        old = print_min(myrun.min)
        new = print_min(min_)

    elif dim and myrun:
        old = print_dim(myrun.dim)
        new = print_dim(dim)

    elif noise is not None and myrun:
        old = print_noise(myrun.noise_function, myrun.noise, myrun.noise_time)
        new = print_noise(myrun.noise_function, noise, myrun.noise_time)

        if myrun.noise == 0 or myrun.devo_time == inf:
            suffix = print_devo_time(default_devo_time)

    elif samples:
        old = print_samples(default_samples)
        new = print_samples(samples)

    elif devo_time:
        old = print_devo_time(default_devo_time)
        new = print_devo_time(devo_time)

    else:
        return filename

    return filename.replace(old, new) + suffix

def get_replaced_suffix(filename, samples=[], devo_times=[]):

    filenames = [get_replaced_filename(filename, samples=sample)
                 for sample in samples]

    filenames += [get_replaced_filename(filename, devo_time=devo_time)
                  for devo_time in devo_times]

    return set(filenames) if filenames else [filename]

def get_suffix(
        samples=None, full_enum=False, symm=False, stable=False, binary=False,
        _filter=None, fix_gpmap=False, gpmap=None, steepness=a, steps='uniform',
        converge_time=T, converge_threshold=converge, devo_time=inf,
        act_fraction=p, p=None, q=None, period=1,
        n_mutants=None, all_mutants=None, deletion=False, change_sign=False,
        mut_bias=False, network=None, knockouts=None, fixed_states=None,
        random=False, get_initials=True, graph=False, scale_free=False,
        genotype_func='', pop_size=None, run=None, experiment='',
        deterministic=True):

    if full_enum:
        suffix = 'binary-enum'
        if symm:
            suffix += '-symm'
        else:
            suffix += '-full'
    else:
        if network:
            suffix = network
            if knockouts:
                suffix += str(knockouts)
            elif fixed_states:
                suffix += str(fixed_states)
            elif random:
                suffix += '-random'
        elif graph:
            suffix = 'graph'
            if scale_free:
                suffix += '-scale_free'
            else:
                suffix += '-exp_pow'
        elif stable:
            suffix = 'stable'
        elif binary:
            suffix = 'binary'
        else:
            suffix = 'random'
        if samples:
            suffix += print_samples(samples)
        if binary:
            if q != None:
                if p != None:
                    suffix += '-p%.1f' %p if type(p) != str else p
                suffix += '-q%.2f' %q if type(q) != str else q
            else:
                if act_fraction != None:
                    suffix += '-p%.1f' %act_fraction
                if p != None:
                    suffix += '-p%.1f' %p

    if n_mutants or all_mutants:
        suffix += get_mut_suffix(
            deletion, change_sign, mut_bias, n_mutants, all_mutants)

    if fix_gpmap:
        gpmap = gpmap.__name__ if gpmap.__name__ != 'gpmap' else gpmap.__doc__
        suffix += '-%s' %gpmap
        if 'cSigmoid' in gpmap:
            suffix += ('_a%.2fT%dE%.0eM%s' %
                       (steepness, converge_time, converge_threshold,
                        devo_time))

    if steps == 'quantiles':
        suffix += '-quantiles'

    if _filter is not None:
        suffix += '-filter%g' %_filter

    if not (fix_gpmap and 'cSigmoid' in gpmap):
        if type(devo_time) is str:
            suffix += '-' + devo_time +'_devo_times'
        elif devo_time != inf:
            suffix += '-devo_time%g' %devo_time

    if not get_initials:
        suffix += '-initials'

    if 'diagonal' in genotype_func:
        suffix += '-%s' %''.join(genotype_func.split('_')[:2])

    if 'qp' in genotype_func and period and period > 1:
        suffix += '-period%s' %period

    if (('qp' in genotype_func or 'diagonal' in genotype_func) and
        'evolution' not in experiment):
        if pop_size:
            suffix += print_generations(pop_size, False).replace('G', 'P')
        if run:
            suffix += print_run(run)

    if (('qp' in genotype_func or 'diagonal' in genotype_func) and not
        deterministic):
        suffix += '-normal%s'%qp_scale

    return suffix

def get_basename(base, precision=2):
    return array2string(base, precision=precision)

def print_density(density, precision=2):
    if precision or density == 1.:
        return get_basename(array([density]), precision)
    else:
        return '[ %.1g]' %density

def print_int_or_float(
        x, precision=True, double_precision=False, threshold=1000):
    if x <= threshold:
        return '%d' %int(x)
    elif double_precision:
        return '%.2e' %x
    elif precision:
        return '%.1e' %x
    else:
        return '%.0e' %x

def print_generations(
        generations, precision=True, double_precision=False, G_string=1000):
    return '_G' + print_int_or_float(
        generations, precision, double_precision, G_string)

def print_N(bits):
    return '-N_%s' %[bits]

def print_c(density, precision):
    return '_c_%s' %print_density(density, precision)

def print_min(min_):
    return '_min_%s' %[min_]

def print_dim(dim):
    return '_dim_%s' %[dim]

def print_noise(noise_function, noise, noise_time):
    filename = '-'

    if noise_function != 'random_noise':
        filename += noise_function
    else:
        filename += 'noise'

    filename += str(noise)

    if noise_time != 'before':
        filename += noise_time

    return filename + '-'
    #'-%s%s%s-' %(myrun.noise_function, myrun.noise, myrun.noise_time)

def print_devo_time(devo_time):
    return '-devo_time%g' %devo_time

def print_run(run):
    return '-n%s' %run

def print_samples(samples):
    return print_generations(samples, False, False).replace('_G', '-n')

def get_filename(bit=N, density=c, min_=-1., dim=2, noise=0, precision=2):
    return ('-N_%s_c_%s_min_%s_dim_%s-noise%s-' %
            ([bit], print_density(density, precision), [min_], [dim], noise))

def get_noise_run_prefix(noise_function, noise, noise_time):
    prefix = ''

    if noise_function == 'random_noise':
        prefix += 'o'
    elif noise_function == 'random_flip':
        prefix += 'rf'
    elif noise_function == 'force_different_flip':
        prefix += 'df'
    elif noise_function == 'neighbour_flip':
        prefix += 'nf'

    prefix += str(noise)

    # default noise_time = 'before' has no prefix
    if noise_time == 'after':
        prefix += 'a'
    elif noise_time == 'shmulevich':
        prefix += 's'

    return prefix

def get_mut_suffix(
        deletion=False, change_sign=False, mut_bias=False,
        n_mutants=None, all_mutants=None):

    if all_mutants:
        suffix = '-all_mut%d' %all_mutants
    elif n_mutants:
        suffix = 'x%dmut' %n_mutants
    else:
        suffix = ''
    if deletion:
        suffix += '-deletion'
    elif change_sign:
        suffix += '-binary'
    elif mut_bias:
        suffix += '-bias%g' %mut_bias
    else:
        suffix += '-random'
    return suffix

def get_evolution_filename(
        model, generations, period, n_neighbours=0, initials=False,
        mutations=False, sel_strength=None, mut_rate=None, rec_rate=None,
        run=None, mut_suffix=None, rec_function=None, precision=False,
        G_string=1000):

    filename = 'model%s'  %model
    filename += print_generations(generations, precision, False, G_string)
    if period:
        filename += '_period%s' %period
    if model == 4:
        filename += '_s%.1g' %sel_strength
    if initials:
        filename += '_initials'
    elif mutations:
        filename += '_mutations'
    if n_neighbours:
        filename += '%s' %n_neighbours
    if mut_rate:
        filename += '_u%.1g' %mut_rate
    if mut_suffix:
        filename += mut_suffix
    if rec_rate:
        filename += '_r%.1g' %rec_rate
    if rec_function:
        filename += '-' + rec_function
    if run:
        filename += print_run(run)
    return filename

def get_evolution_shortname(
        model, generations=None, period=None,
        sel_strength=None, mut_rate=None, rec_rate=None,
        shortname=None, run=None, g=None, precision=False,
        G_string=1000):

    filename = 'm%s' %model
    if generations:
        filename += 'G%s' %print_generations(
            generations, False, False, G_string)[2:]
    if period:
        filename += 'p%s' %period
    if model == 4 and sel_strength:
        filename += 's%.1g' %sel_strength
    if mut_rate:
        filename += 'u%.1g' %mut_rate
    if rec_rate:
        filename += 'r%.1g' %rec_rate
    if shortname:
        filename = shortname
    if run:
        filename += print_run(run)[1:]
    if g:
        filename += 'g%s' %print_generations(g, precision, False, G_string)[2:]
    return filename


############################
## RegulonDB's E. coli TRN #
############################

# ('source', 'target', 'link', 'comments')
network_files = {
    'network_tf_tf8.2':   {'columns': (0,1,2), 'comments': 37,
                           'release': '7.4 Date: April-22-2013'},
    'network_tf_tf7.4':   {'columns': (0,1,2), 'comments': 31+1,
                           'release': '7.4 Date: MARCH-29-12'},
    'network_tf_tf7.3':   {'columns': (0,1,2), 'comments': 30+1,
                           'release': '7.3 Date: 01-NOV-11'},
    'network_tf_tf7.2':   {'columns': (0,1,2), 'comments': 30+1,
                           'release': '7.2 Date: 06-MAY-11'},
    'network_tf_tf7.0':   {'columns': (0,1,2), 'comments': 34+1,
                           'release': '7.0 Date: Jan-26-11'},
    'network_tf_gene':    {'columns': (0,1,2), 'comments': 31,
                           'release': '7.3 Date: 01-NOV-11'},
    'network_tf_gene7.2': {'columns': (1,3,5), 'comments': 34,
                           'release': '7.2 Date: 06-MAY-11'},
    'network_tf_gene7.0': {'columns': (1,3,5), 'comments': 39,
                           'release': '7.0 Date: Jan-26-11'},
    'sRNADataSet':        {'columns': (1,2,3), 'comments': 38,
                           'release': '7.2 Date: 06-MAY-11'},
    'TFSet':              {'columns': (1,2),   'comments': 31,
                           'release': '7.3 Date: 01-NOV-11'},
    'TFSet7.2':           {'columns': (1,2)},
    'S1':                 {'columns': (0,5,10,11)}}

net_short_names = {'tf_tf':      ['network_tf_tf8.2'],
                   'tf_tf7.4':   ['network_tf_tf7.4'],
                   'tf_tf7.3':   ['network_tf_tf7.3'],
                   'tf_tf7.2':   ['network_tf_tf7.2'],
                   'tf_tf7.0':   ['network_tf_tf7.0'],
                   'tf_gene':    ['network_tf_gene'],
                   'tf_gene7.2': ['network_tf_gene7.2'],
                   'tf_gene7.0': ['network_tf_gene7.0'],
                   'sRNA':       ['sRNADataSet'],
                   'meta2':      ['network_tf_tf', 'network_tf_gene'],
                   'meta3':      ['network_tf_tf', 'network_tf_gene',
                                  'sRNADataSet']}

regulon_log_file = 'RegulonDB.log'
regulon_err_file = 'RegulonDB.err'

def read_network(fname='network_tf_tf', subdir='ecoli/', skiprows=0):
    return loadtxt(
        nets_dir + subdir + fname + '.txt', str, skiprows=skiprows,
        usecols=network_files[fname]['columns'])

def filter_network(network, _filter):
    return delete(network, where(network == _filter)[0], axis=0)

def print_some_network_statistics(
        short_name, genotype, network, n_tfs, not_regulated, candidates, n, r,
        knockouts = None):
    # number of nodes/genes, i.e., size of the network
    bits = len(genotype)
    # number of edges/interactions
    M = nonzero(genotype)[0].size
    # density of the network
    density = float(M)/bits**2
    # average degree of the network
    degree = float(M)/bits
    # number of regulators
    regulators = n['regulators']; del n['regulators']
    # self-regulation
    self_reg = [where(genotype.diagonal() == link)[0].size for link in (-1, 1, 0)]
    # error control
    lines   = n['lines']; del n['lines']
    repeats = n['repeats']; del n['repeats']
    errors  = n['errors']; del n['errors']
    # file output
    fnames = net_short_names[short_name]
    file = open(nets_dir + regulon_log_file,  'a')
    print >> file, '# Release: ' + network_files[fnames[0]]['release']
    print >> file, "printing to log (from %s)" %fnames
    print >> file, '%d lines in file - %d repeats = %d interactions ... %s! %d errors' %(lines, repeats, M, lines - repeats == M, errors)
    print >> file, 'N = %d genes, M = %d interactions, density = %.1g, degree = %.1f' %(bits, M, density, degree)
    print >> file, '%d genes coding for TFs, %d genes actually regulating' %(n_tfs, regulators)
    print >> file, '%d genes are not regulated; self-regulation [-,+,0] = %s' %(not_regulated, str(self_reg))
    if n['+-']:   print >> file, 'attention: self-regulation may include +- dual interactions'
    print >> file, 'genes that are only self-regulated and regulate at least another gene:'
    print >> file, candidates
    if knockouts: print >> file, 'knockouts: %s' %knockouts
    print >> file, n
    print >> file, r
    print >> file, 'There is a total number of %.2g possible combinations of this network' %(2**(n['+-'] + n['?']))
    print >> file, ''
    file.close()


##################
## OLD FUNCTIONS #
##################

# save ensemble to a binary file
def save_ensemble(filename, ensemble):
    file = open(filename, 'wb')
    for triplet in ensemble:
        np.save(file, triplet[0]) # initial phenotype
        np.save(file, triplet[1]) # optimum phenotype
        pop_size = array(len(triplet[2]))
        np.save(file, pop_size)
        for individual in triplet[2]: # population of genotypes
            np.save(file, individual)
    file.close()

# load ensemble to a binary file
def load_ensemble(filename, size):
    ensemble = []
    file = open(filename)
    for i in range(size):
        initial = np.load(file)
        optimum = np.load(file)
        pop_size = np.load(file)
        population = []
        for i in range(pop_size):
            population.append(np.load(file))
        ensemble.append((initial, optimum, population))
    file.close()
    return ensemble

# read phenotypes from compressed array file (.npz)
def load_phenotypes(filename):
    file = open(filename)
    phenotypes = np.load(file)['arr_0']
    file.close()
    return phenotypes

#
def save_list(filename, list):
    file = open(filename, 'w')
    file.writelines(', '.join(map(str,list)))
    file.close()

#
def load_list(filename):
    reader = csv.reader(open(filename))
    list = []
    list.extend(reader)
    new_list = []
    for string in list[0]:
        new_list.append(float(string))
    return new_list

# correct for small samples
def correct_sample_size(x, y, yerr, data, runs_axis=0, threshold=.5):
    masked  = ma.count_masked(data, axis=runs_axis)
    if shape(masked) != shape(x):
        masked = masked.T[0]
    runs = shape(data)[runs_axis]
    samples = runs - masked
    i = where(samples < runs*threshold)[0]
    # don't delete everything:
    if i.size == x.size:
        i = where(samples < samples.max())[0]
    return delete(x, i), delete(y, i), delete(yerr, i)

def correct_sample_size2(data, runs_axis=0, threshold=.49):
    non_masked = ma.count(data, axis=runs_axis)
    if non_masked.ndim != 1:
        non_masked = non_masked.T[0]
    n_runs = shape(data)[runs_axis]
    return where(non_masked > n_runs*threshold)[0]
