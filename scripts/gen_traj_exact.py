# -*- coding: utf-8 -*-
"""
@author: davidsantiagoquevedo
"""
import sys
h = float(sys.argv[1])
n = int(sys.argv[2])
T = int(sys.argv[3])
trajectories = int(sys.argv[4])

DATA_PATH = "data/"

sys.path.append("src/")
import fBm_exact as fbml
import utils as ut

df_fbm, df_fgn = fbml.frac_brown_davies_harte(h, n, T, trajectories)

ut.write_hdf5(DATA_PATH + f"fBm_dh-h-{h}-{trajectories}.hdf5", df_fbm)
ut.write_hdf5(DATA_PATH + f"fgn_dh-h-{h}-{trajectories}.hdf5", df_fbm)