# -*- coding: utf-8 -*-
"""
@author: davidsantiagoquevedo
"""
from fbm import FBM
import pandas as pd
from tqdm import tqdm

def frac_brown_davies_harte(h, n, T, trajectories = 100):
    """An accurate fractional Brownian motion generator (1994)
    Simulation of benchmark trajectories using fbm library.
    (See: https://doi.org/10.1016/0378-4371(94)90531-2, https://github.com/crflynn/fbm)

    Parameters
    ----------
    h : float64
        Hurst exponente
    n: int64
        Number of data points
    T: int64
        lenght of the trajectory (max time)
    n: int64
        Number of data points
    trajectories: int64
        Number of trajectories (ensembles) generated
    """
    df_noise = {}
    df_traj = {}
    for tr in tqdm(range(trajectories)):
        f = FBM(n = n, hurst = h, length = T, method = 'daviesharte')
        if tr == 0:
            df_traj.update({"t" : f.times()})
        # Generate a fBm realization
        df_traj.update({f"traj_{tr}" : f.fbm()})
        df_noise.update({f"ns_{tr}" : f.fgn()})
    return pd.DataFrame(df_traj), pd.DataFrame(df_noise)