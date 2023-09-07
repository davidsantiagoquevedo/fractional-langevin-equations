# -*- coding: utf-8 -*-
"""
@author: davidsantiagoquevedo
"""
import pandas as pd
import numpy as np

def msd(trajectories, normalize = True):
    """Numeric estimation of the mean square displacement (MSD)
	
	Args:
	    trajectories (panda.DataFrame): data table with trajectories
	    normalize (boolean). Defaults to True: If True, normalizes the MSD with respect to its maximum value.
    
    Returns:
        pandas.Series: mean squared displacement
	"""
    msd = (trajectories**2).mean(axis = 1)
    if normalize:
        msd /= max(msd)
    return msd

def cov(series, t_ref, normalize = True, limit = False):
    """Numeric estimation of the (auto)covariance with respect to a reference point in time.
	
	Args:
	    series (pandas.DataFrame): data table with series
	    t_ref (int): reference time point for the covariance
	    normalize (boolean). Defaults to True: If True, normalizes the covariance with respect to its maximum value.
        limit (boolean). Defaults to False: calculate the autocovariance assuming the mean vanishes for every realization of the series,
    
    Returns:
        pandas.Series: covariance with respect to time point t_ref
	"""
    mean_t = series.mean(axis = 1)
    mean_tref = mean_t.loc[t_ref]

    X_tref = series.loc[t_ref]
    X_t_minus_mean = series.apply(lambda x: x - np.array(mean_t))
    X_tref_minus_mean = X_tref - mean_tref

    autocov = (X_t_minus_mean*X_tref_minus_mean).mean(axis = 1)
    
    if limit:
        autocov = (series*X_tref).mean(axis = 1)
    if normalize:
        autocov /= max(autocov)
    return autocov