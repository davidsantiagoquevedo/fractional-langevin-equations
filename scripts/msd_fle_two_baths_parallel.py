import multiprocessing as mp
import sys
import numpy as np
import ctypes as c
from tqdm import tqdm
import pandas as pd

sys.path.append("src/")
from fLe_twobath import fle_twobath
import fBm_stats as fbs
import utils as ut
import datetime

import warnings
warnings.filterwarnings("ignore") 
#UserWarning: Combination of increments n and Hurst value H invalid for Davies-Harte method. Reverting to Hosking method. Occurs when n is small and Hurst is close to 1. 
#warnings.warn

#DATA_PATH = "data/two_baths/"
DATA_PATH = "dolab/"

# Input
H = float(sys.argv[1])
A = float(sys.argv[2])
eta = float(sys.argv[3])
C = float(sys.argv[4])
theta_12 = float(sys.argv[5])
theta_H = float(sys.argv[6])

#Fixed params
h = 0.1
T = 20
n = int(T/h)
realizations = 40

save_all = True

batch_size = 4
assert(realizations%batch_size == 0)


task = f"msd-h-{H}-A{A}-eta{eta}-C{C}-t12_{theta_12}-tH{theta_H}"

t = datetime.datetime.now()
print(f"Start task {task}: {t.year}/{t.month}/{t.day} {t.hour}:{t.minute}:{t.second}")

def to_shared_array(arr, ctype):
    shared_array = mp.Array(ctype, arr.size, lock=False)
    temp = np.frombuffer(shared_array, dtype=arr.dtype)
    temp[:] = arr.flatten(order = "C")
    return shared_array

def to_numpy_array(shared_array, shape):
    """Create a numpy array backed by a shared memory Array."""
    arr = np.ctypeslib.as_array(shared_array)
    return arr.reshape(shape)

def msd(shared_msd, shape_msd, lock, process, H, T, A, eta, C, theta_12, theta_H, h = 0.01):   
    temp_msd = to_numpy_array(shared_msd, shape_msd)
    
    eq = fle_twobath(H)
    eq.params(T, h, v0 = 0, 
              A = A, eta = eta, C = C, 
              theta_12 = theta_12, theta_H = theta_H)
    eq.make_B_H()
    eq.solve()
    n = int(T/h)
    msd_i = eq.numerical**2
    with lock:
        if save_all:
            temp_msd[process] += eq.numerical
        else:
            temp_msd += msd_i
        
    
    

# protect the entry point
if __name__ == '__main__':
    if save_all:
        init_msd = np.zeros((realizations, n), dtype = np.float32)
    else:
        init_msd = np.zeros(n, dtype = np.float32)
    shared_msd = to_shared_array(init_msd, c.c_float)
    stored_msd = to_numpy_array(shared_msd, init_msd.shape)
    lock = mp.Lock()
    # execute in batches
    for i in tqdm(range(0, realizations, batch_size)):
        #print(f"Start thread {i}: {t.year}/{t.month}/{t.day} {t.hour}:{t.minute}:{t.second}")
        # execute all tasks in a batch
        processes = [mp.Process(target = msd, 
                                args=(shared_msd, init_msd.shape, lock, p,
                                      H, T, A, eta, C, theta_12, theta_H, h)) 
                     for p in range(i, i + batch_size)]
        # start all processes
        for process in processes:
            process.start()
        # wait for all processes to complete
        for process in processes:
            process.join()
        #print(f"End thread {i}: {t.year}/{t.month}/{t.day} {t.hour}:{t.minute}:{t.second}")
    # report that all tasks are completed
    print('Done', flush=True)

if save_all:
    t_ = np.arange(0, T, h)
    d = dict(enumerate(stored_msd, 1))
    df_msd = pd.DataFrame(d)
    df_msd.columns = "trj_"+ df_msd.columns.astype(str)
    df_msd["t"] = t_
else:
    t_ = np.arange(0, T, h)
    df_msd = pd.DataFrame(dict(zip(["t", "msd"], [t_, stored_msd/realizations])))

ut.write_hdf5(DATA_PATH + f"par_dt{h}-T{T}_{task}.hdf5", df_msd)
t = datetime.datetime.now()
print(f"End task {task}: {t.year}/{t.month}/{t.day} {t.hour}:{t.minute}:{t.second}") 
