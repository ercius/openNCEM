"""
Module containing definitions of basic math functions, which can be used for fitting.

If you add functions here, do not forget to update the lookup table at the bottom as well.
"""

import numpy as np
import scipy.special


def const( x, param):
    """Constant function.

    f = param[0] .

    Parameters:
        x (np.ndarray):    Positions at which to evaluate function.
        param (np.ndarray):    Necessary parameters.

    Returns:
        (np.ndarray):    Values of function at x.

    """
    
    return param[0]*np.ones(x.shape)


def linear( x, param):
    """Linear function.

    f = param[0] * x + param[1] .

    Parameters:
        x (np.ndarray):    Positions at which to evaluate function.
        param (np.ndarray):    Necessary parameters.

    Returns:
        (np.ndarray):    Values of function at x.

    """
    
    return param[0]*x + param[1]


def powlaw( x, param):
    """Power law.

    f = param[0] * x^param[1] .

    Parameters:
        x (np.ndarray):    Positions at which to evaluate function.
        param (np.ndarray):    Necessary parameters.

    Returns:
        (np.ndarray):    Values of function at x.

    """
    return param[0]*np.power(x, param[1])
    

def voigt(x, param):
    """
    Voigt peak function.

    f = param[0]/(param[2] * sqrt(2*pi)) * Re( Faddeeva( ((x - param[1]) + i*param[3])/(param[2] * sqrt(2)) ) )

    Parameters:
        x (np.ndarray):    Positions at which to evaluate function.
        param (np.ndarray):    Necessary parameters.

    Returns:
        (np.ndarray):    Values of function at x.

    """
    # A = param[0]
    # mu = param[1]
    # sigma = param[2]
    # gamma = param[3]
    return param[0]*np.real(scipy.special.wofz((x-param[1] + 1.0j*param[3])/(param[2]*np.sqrt(2.0))))/(param[2]*np.sqrt(2.0*np.pi))


def sum_functions( x, funcs, param ):
    """
    Sum functions in funcs over range x.

    Parameters:
        x (np.ndarray):    Positions at which to evaluate the functions.
        funcs (list):    List of strings identifying function implemented in ncempy.algo.math.
        param (np.ndarray):    Concatenated parameters for functions in funcs.

    Returns:
        (np.ndarray):   Values of sum of functions at x.

    """
    
    # estimator
    est = np.zeros(x.shape)
    
    n = 0
    # evaluate given functions
    for i in range(len(funcs)):
        est += lkp_funcs[funcs[i]][0]( x, param[n:n+lkp_funcs[funcs[i]][1]] )
        n += lkp_funcs[funcs[i]][1]
        
    return est
    

# lookup table for functions
lkp_funcs = {'const': (const, 1),
              'linear': (linear, 2),
              'powlaw': (powlaw, 2),
              'voigt': (voigt, 4)
            }
'''(dict):    Look-up table for functions implemented in this module. Functions are identified by strings, the entries give handles to the functions as well as the number of arguments necessary.'''
