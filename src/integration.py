import types
import numpy as np

def convolution(f, g, t):
    """Computes the convolution between an analytical function f(x) and a secondary function
        that can be: another analytical function or a time series (e.a. noise series) 

    Args:
        f (function): main function to convolve
        g (function or ndarray): secondary function to convolve
        t (ndarray): vector with the domain of integration

    Returns:
        ndarray: result of the convolution along the domain t
    """
    N = len(t)
    if isinstance(g, np.ndarray):
        g_eval = g
    elif isinstance(f, types.FunctionType):
        g_eval = g(t)
    dt = t[1] - t[0]
    conv = np.zeros(N)
    for i in range(N):
        dt = t[1] - t[0]
        g_eval_i = g_eval[:i+1] # subset of g(t)
        t_i = t[:i+1] # subset of the time vector
        integrand = f(t[i] - t_i) * g_eval_i
        conv[i] = np.trapz(integrand, dx = dt)
    return conv