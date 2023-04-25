# -*- coding: utf-8 -*-
"""
@author: davidsantiagoquevedo
"""
import h5py
import os
from tqdm import tqdm
import pandas as pd

def read_hdf5_data(path):
    """Read hdf5 structered in columns
    Adapted from:
        https://stackoverflow.com/questions/28170623/how-to-read-hdf5-files-in-python    
        https://stackoverflow.com/questions/65865756/how-extract-data-from-hdf5-in-python
  
    Args:
        path (str): path to file

    Returns:
        pandas.DataFrame: column-structured data    
    """
    with h5py.File(path, "r") as f:
        group_key = list(f.keys())[0]
        data = list(f[group_key])
        df = {}
        for col in data:
            df.update({col: f[group_key][col][:]})
    return pd.DataFrame(df)

def read_hdf5_all(h, DATA_PATH, trajectory_i = 1, trajectory_f = 100):
    """Read several hdf5 files that contain trajectories generated with a hurst number
    
    Args:
        h (float64): Hurst number
        DATA_PATH (str): path to folder where the data is located
        trajectory_i (int). Defaults to 1: first trajectory to read
        trajectory_i (int). Defaults to 100: last trajectory to read
    
    Returns:
        pandas.DataFrame: concat table with multiple trajectories
    """
    dictionary_data = {}
    for tr in tqdm(range(trajectory_i, trajectory_f + 1)):
        fbm_sample = read_hdf5_data(DATA_PATH + f"fBM-h-{h}-{tr}.hdf5")
        if len(dictionary_data) == 0:
            dictionary_data.update({"t" : fbm_sample["deets_t"]})
        dictionary_data.update({"traj_" + str(tr) : fbm_sample["deets_v"]})
    return pd.DataFrame(dictionary_data)

def write_hdf5(path, data):
    """Write hdf5 data from pandas dataframes
    
    Args:
        path (str): path to new file
        data (pandas.DataFrame): data table to write. The names of all the columns must be strings!
    """
    try:
        os.remove(path)
    except:
        print("File does not exist. Creating new file...")
    with h5py.File(path, "a") as data_file:
        data_file.create_group("values")
        group = data_file["values"]
        for col in data.columns:
            group[col] = data[col]
    data_file.close()    