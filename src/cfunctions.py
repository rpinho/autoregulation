import header
reload(header)
from header import *

# TODO: consider using python ctypes instead of scypy.weave

# C's implementation of sigmoid: 1 / ( 1 + exp(-a*x) )
def cSigmoid(input, steepness):
    """cSigmoid"""
    size    = input.size
    output  = input.copy()
    support = """#include <math.h>"""
    code = """
           for (int i=0; i<size; ++i)
               output(i) = 1 / (1 + exp(-steepness*output(i)));
           return_val = py_output;
           """
    return weave.inline(
        code, ['output','size','steepness'],
        type_converters=weave.converters.blitz,
        support_code=support, libraries=['m'])

# C's implementation of sigmoid2: 2 / ( 1 + exp(-a*x) ) - 1
def cSigmoid2(input, steepness):
    """cSigmoid2"""
    size    = input.size
    output  = input.copy()
    support = """#include <math.h>"""
    code = """
           for (int i=0; i<size; ++i)
               output(i) = 2 / (1 + exp(-steepness*output(i))) - 1;
           return_val = py_output;
           """
    return weave.inline(
        code, ['output','size','steepness'],
        type_converters=weave.converters.blitz,
        support_code=support, libraries=['m'])

# 0 -> -1 or 1 randomly
def csign_r(input):
    """csign_r"""
    size    = input.size
    output  = input.copy()
    support = """
              #include <stdio.h>
              #include <time.h>
              """
    code    = """
              unsigned int iseed = (unsigned int)time(NULL);
              srandom(iseed);

              for (int i=0; i<size; ++i)
                 if (output(i) < 0)
	            output(i) = -1;
	            else if (output(i) > 0)
	               output(i) = 1;
	               else
                          output(i) = random()&01 == 1 ? -1 : 1;
              return_val = py_output;
              """
    return weave.inline(
        code, ['output','size'], type_converters=weave.converters.blitz)

# 0 -> 1
def csign(input):
    """csign"""
    size   = input.size
    output = input.copy()
    code = """
           for (int i=0; i<size; ++i)
               output(i) = output(i) < 0 ? -1 : 1;
           return_val = py_output;
           """
    return weave.inline(
        code, ['output','size'], type_converters=weave.converters.blitz)

# 0 -> -1
def csign0(input):
    """csign0"""
    size   = input.size
    output = input.copy()
    code = """
           for (int i=0; i<size; ++i)
               output(i) = output(i) <= 0 ? -1 : 1;
           return_val = py_output;
           """
    return weave.inline(
        code, ['output','size'], type_converters=weave.converters.blitz)

# 0 -> 0
def csign101(input):
    """csign101"""
    return sign(input)
'''
    size    = input.size
    output  = input.copy()
    code    = """
              for (int i=0; i<size; ++i)
                 if (output(i) < 0)
	            output(i) = -1;
	            else if (output(i) > 0)
	               output(i) = 1;
	               else
                          output(i) = 0;
              return_val = py_output;
              """
    return weave.inline(code, ['output','size'], type_converters = weave.converters.blitz)
'''
# 0 -> s(t)
def csign_t(input):
    """csign_t"""
    size    = input.size
    output  = input.copy()
    code    = """
              for (int i=0; i<size; ++i)
                 if (output(i) < 0)
	            output(i) = -1;
	            else if (output(i) > 0)
	               output(i) = 1;
	               else
                          output(i) = output(i);
              return_val = py_output;
              """
    return weave.inline(
        code, ['output','size'], type_converters=weave.converters.blitz)

# 0 -> 1
def cHeaviside(input):
    """cHeaviside"""
    size   = input.size
    output = input.copy()
    code = """
           for (int i=0; i<size; ++i)
               output(i) = output(i) < 0 ? 0 : 1;
           return_val = py_output;
           """
    return weave.inline(
        code, ['output','size'], type_converters=weave.converters.blitz)

# 0 -> 0
def cHeaviside0(input):
    """cHeaviside0"""
    size   = input.size
    output = input.copy()
    code = """
           for (int i=0; i<size; ++i)
               output(i) = output(i) <= 0 ? 0 : 1;
           return_val = py_output;
           """
    return weave.inline(
        code, ['output','size'], type_converters=weave.converters.blitz)

def cstep(input, jump, state):
    n_states = state.size
    n_jumps  =  jump.size
    n_bits   = input.size
    output   = input.copy()
    support = """
           #include <iostream>

           int Heaviside(double x) {
             return x < 0 ? 0 : 1;
           }
           double step(double x, int n, double *state, double *jump) {
             return n == 0 ? state[0] : (state[n] - state[n-1]) * Heaviside(x-jump[n-1]) + step(x, n-1, state, jump);
           }
           /*void print_matrix(double *matrix, int n) {
             for (int i=0; i<n; i++)
               std::cout << matrix[i] << std::endl;
           }*/
           """
    code = """
           double states[n_states];
           for (int i=0; i<n_states; ++i)
             states[i] = state(i);
           //print_matrix(states, n_states);

           double jumps[n_jumps];
           for (int i=0; i<n_jumps; ++i)
             jumps[i] = jump(i);
           //print_matrix(jumps, n_jumps);

           for (int i=0; i<n_bits; ++i)
             output(i) = step(output(i), n_states-1, states, jumps);

           return_val = py_output;
           """
    return weave.inline(
        code, ['jump','state','output','n_states', 'n_jumps', 'n_bits'],
        type_converters=weave.converters.blitz, support_code=support)

# k = 2, same as csign or cHeaviside: 0 -> 1
def cstep1(input, jump, state):
    size   = input.size
    output = input.copy()
    code = """
           for (int i=0; i<size; ++i)
             if (output(i) < jump(0))
	       output(i) = state(0);
	       else
	         output(i) = state(1);
           return_val = py_output;
           """
    return weave.inline(
        code, ['jump','state','output','size'],
        type_converters=weave.converters.blitz)

# k = 3
def cstep2(input, jump, state):
    size   = input.size
    output = input.copy()
    code = """
           for (int i=0; i<size; ++i)
             if (output(i) < jump(0))
	       output(i) = state(0);
	       else if (output(i) < jump(1))
	         output(i) = state(1);
	         else
		   output(i) = state(2);
           return_val = py_output;
           """
    return weave.inline(
        code, ['jump','state','output','size'],
        type_converters=weave.converters.blitz)

def cstep3(input, jump, state):
    size   = input.size
    output = input.copy()
    code = """
           for (int i=0; i<size; ++i)
             if (output(i) < jump(0))
	       output(i) = state(0);
	       else if (output(i) < jump(1))
	         output(i) = state(1);
	         else if (output(i) < jump(2))
		   output(i) = state(2);
                   else
		     output(i) = state(3);
           return_val = py_output;
           """
    return weave.inline(
        code, ['jump','state','output','size'],
        type_converters=weave.converters.blitz)

def cstep4(input, jump, state):
    size   = input.size
    output = input.copy()
    code = """
           for (int i=0; i<size; ++i)
             if (output(i) < jump(0))
	       output(i) = state(0);
	       else if (output(i) < jump(1))
	         output(i) = state(1);
	         else if (output(i) < jump(2))
		   output(i) = state(2);
                   else if (output(i) < jump(3))
		     output(i) = state(3);
                     else
		       output(i) = state(4);
           return_val = py_output;
           """
    return weave.inline(
        code, ['jump','state','output','size'],
        type_converters=weave.converters.blitz)

def cstep5(input, jump, state):
    size   = input.size
    output = input.copy()
    code = """
           for (int i=0; i<size; ++i)
             if (output(i) < jump(0))
	       output(i) = state(0);
	       else if (output(i) < jump(1))
	         output(i) = state(1);
	         else if (output(i) < jump(2))
		   output(i) = state(2);
                   else if (output(i) < jump(3))
		     output(i) = state(3);
                     else if (output(i) < jump(4))
		       output(i) = state(4);
                       else
		         output(i) = state(5);
           return_val = py_output;
           """
    return weave.inline(
        code, ['jump','state','output','size'],
        type_converters=weave.converters.blitz)

def cstep6(input, jump, state):
    size   = input.size
    output = input.copy()
    code = """
           for (int i=0; i<size; ++i)
             if (output(i) < jump(0))
	       output(i) = state(0);
	       else if (output(i) < jump(1))
	         output(i) = state(1);
	         else if (output(i) < jump(2))
		   output(i) = state(2);
                   else if (output(i) < jump(3))
		     output(i) = state(3);
                     else if (output(i) < jump(4))
		       output(i) = state(4);
                       else if (output(i) < jump(5))
		         output(i) = state(5);
                         else
		           output(i) = state(6);
           return_val = py_output;
           """
    return weave.inline(
        code, ['jump','state','output','size'],
        type_converters=weave.converters.blitz)

def cstep7(input, jump, state):
    size   = input.size
    output = input.copy()
    code = """
           for (int i=0; i<size; ++i)
             if (output(i) < jump(0))
	       output(i) = state(0);
	       else if (output(i) < jump(1))
	         output(i) = state(1);
	         else if (output(i) < jump(2))
		   output(i) = state(2);
                   else if (output(i) < jump(3))
		     output(i) = state(3);
                     else if (output(i) < jump(4))
		       output(i) = state(4);
                       else if (output(i) < jump(5))
		         output(i) = state(5);
                         else if (output(i) < jump(6))
		           output(i) = state(6);
                           else
		             output(i) = state(7);
           return_val = py_output;
           """
    return weave.inline(
        code, ['jump','state','output','size'],
        type_converters=weave.converters.blitz)

def cstep8(input, jump, state):
    size   = input.size
    output = input.copy()
    code = """
           for (int i=0; i<size; ++i)
             if (output(i) < jump(0))
	       output(i) = state(0);
	       else if (output(i) < jump(1))
	         output(i) = state(1);
	         else if (output(i) < jump(2))
		   output(i) = state(2);
                   else if (output(i) < jump(3))
		     output(i) = state(3);
                     else if (output(i) < jump(4))
		       output(i) = state(4);
                       else if (output(i) < jump(5))
		         output(i) = state(5);
                         else if (output(i) < jump(6))
		           output(i) = state(6);
                           else if (output(i) < jump(7))
		             output(i) = state(7);
                             else
		               output(i) = state(8);
           return_val = py_output;
           """
    return weave.inline(
        code, ['jump','state','output','size'],
        type_converters=weave.converters.blitz)
