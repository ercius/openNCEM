'''
Module containing definitions of basic math functions, which can be used for fitting.

If you add functions here, do not forget to update the lookup table at the bottom as well.
'''

import numpy as np
import scipy.special


def const( x, param):
    '''
    Constant function.
    '''
    
    return param[0]*np.ones(x.shape)

def linear( x, param):
    '''
    Linear function.
    '''
    
    return param[0]*x + param[1]

def powlaw( x, param):
    '''
    Power law.
    '''
    # A = param[0]
    # n = param[1]
    ## A*x^n
    return param[0]*np.power(x, param[1])
    

def voigt( x, param):
    '''
    Voigt peak function.
    '''
    # A = param[0]
    # mu = param[1]
    # sigma = param[2]
    # gamma = param[3]
    return param[0]*np.real(scipy.special.wofz((x-param[1] + 1.0j*param[3])/(param[2]*np.sqrt(2.0))))/(param[2]*np.sqrt(2.0*np.pi))


def sum_functions( x, funcs, param ):
    '''
    Sum functions in funcs over range x.
    
    input:
    - x             axis to evaluate functions on
    - funcs         list of functions to include with entrys (function_handle, number_of_arguments)
    - param         parameters for functions in funcs
    '''
    
    # estimator
    est = np.zeros(x.shape)
    
    n = 0
    # evaluate given functions
    for i in range(len(funcs)):
        est += lkp_funcs[funcs[i]][0]( x, param[n:n+lkp_funcs[funcs[i]][1]] )
        n += lkp_funcs[funcs[i]][1]
        
    return est
    

# lookup table for functions
lkp_funcs = { 'const': (const, 1),
              'linear': (linear, 2),
              'powlaw': (powlaw, 2),
              'voigt': (voigt, 4)
            }
