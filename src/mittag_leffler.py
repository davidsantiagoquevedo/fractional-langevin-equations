# -*- coding: utf-8 -*-
"""
Created on Tue Dec 7 2022
@author: davidsantiagoquevedo
"""
import numpy as np
from scipy.special import gamma
from scipy.integrate import quad

def mittag_leffler_point(z, alpha, beta, inf = 100):
    """This function evaluates the Mittag-Leffler function (MFF) of at a real value z
   
    Parameters
    ----------
    z : float64
        Evaluation point
    alpha : float64
        alpha coefficient of the MFF function. Now limited to real values but should admit complex values
    beta : float64
        beta coefficient of the MFF function. Now limited to real values but should admit complex values
    inf : int64
        Sufficiently large value for the  summation on the definition (in general 100 is enough to vanis 1/gamma dependency)
    """
    s_k = np.zeros(inf)
    k = np.arange(0, inf, 1, dtype=int)
    return (z**k/gamma(alpha*k + beta)).sum()

def mittag_leffler_vector(vect, alpha, beta, inf = 100):
    """This function evaluates the Mittag-Leffler function (MFF) in a 1xn numpy array
    In principle, this is the most computational efficient way of evaluation 
    (See: https://stackoverflow.com/questions/35215161/most-efficient-way-to-map-function-over-numpy-array)
    
    Parameters
    ----------
    vect : 1xn numpy array 
        Evaluation vector
    alpha : float64
        alpha coefficient of the MFF function. Now limited to real values but should admit complex values
    beta : float64
        beta coefficient of the MFF function. Now limited to real values but should admit complex values
    inf : int64
        Sufficiently large value for the  summation on the definition (in general 100 is enough to vanis 1/gamma dependency)
    """
    mtlf = lambda x: mittag_leffler_point(x, alpha, beta, inf)
    v_mtlf = np.vectorize(mtlf)
    return v_mtlf(vect)

def mittag_leffler_stable(z, a):
    """Stable Piecewise implementation for the Mittag-Leffler function
    (See: https://stackoverflow.com/questions/48645381/instability-in-mittag-leffler-function-using-numpy)
    
    Parameters
    ----------
    z : 1xn numpy array 
        Evaluation vector
    a : float64
        alpha coefficient of the MFF function. Now limited to real values but should admit complex values
    """
    z = np.atleast_1d(z)
    if a == 0:
        return 1/(1 - z)
    elif a == 1:
        return np.exp(z)
    elif a > 1 or all(z > 0):
        k = np.arange(100)
        return np.polynomial.polynomial.polyval(z, 1/gamma(a*k + 1))

    # a helper for tricky case, from Gorenflo, Loutchko & Luchko
    def _MLf(z, a):
        if z < 0:
            f = lambda x: (np.exp(-x*(-z)**(1/a)) * x**(a-1)*np.sin(np.pi*a)
                          / (x**(2*a) + 2*x**a*np.cos(np.pi*a) + 1))
            return 1/np.pi * quad(f, 0, np.inf)[0]
        elif z == 0:
            return 1
        else:
            return mittag_leffler_stable(z, a)
    return np.vectorize(_MLf)(z, a)