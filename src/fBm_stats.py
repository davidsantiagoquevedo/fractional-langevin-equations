# -*- coding: utf-8 -*-
"""
@author: davidsantiagoquevedo
"""
import pandas as pd

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

def cov(series, t_ref, normalize = True):
    """Numeric estimation of the covariance with respect to a reference point in time.
	
	Args:
	    trajectories (pandas.DataFrame): data table with trajectories
	    t_ref (int): reference time point for the covariance
	    normalize (boolean). Defaults to True: If True, normalizes the covariance with respect to its maximum value.
    
    Returns:
        pandas.Series: covariance with respect to time point t_ref
	"""
    #TODO: check this implementation. Currently assuming that E[series(t)] = 0 for every t
    mean = series.mean(axis = 1)
    x_ref = series.iloc[t_ref]
    x_ref_mean = mean.iloc[t_ref]
    df_cov = pd.DataFrame()
    for i in range(len(x_ref)):
        df_cov[series.columns[i]] = series[series.columns[i]]*x_ref[i]
    cov = df_cov.mean(axis = 1)
    if normalize:
        cov /= max(cov)
    return cov