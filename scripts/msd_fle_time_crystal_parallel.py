import multiprocessing as mp
import sys
import numpy as np
import ctypes as c
from tqdm import tqdm
import pandas as pd

sys.path.append("src/")
from fLe_timecrystal import fle
import utils as ut
import datetime

import warnings
warnings.filterwarnings("ignore") 
#UserWarning: Combination of increments n and Hurst value H invalid for Davies-Harte method. Reverting to Hosking method. Occurs when n is small and Hurst is close to 1. 
#warnings.warn

#DATA_PATH = "data/two_baths/"
DATA_PATH = "_raw/time_crystal/"

# Input
alpha = float(sys.argv[1])
linear = sys.argv[2]
M = float(sys.argv[3])
v0 = float(sys.argv[4])
eta_1 = float(sys.argv[5])
eta_2 = float(sys.argv[6])
T1 = float(sys.argv[7])
T2 = float(sys.argv[8])

h = float(sys.argv[9])
T = int(sys.argv[10])
n = int(T/h)
avg = int(sys.argv[11])
set = str(sys.argv[12])

save_all = True

batch_size = 4
assert(avg%batch_size == 0)

if save_all:
    prefix = "trj"
else:
    prefix = "msd"

task = f"{prefix}-set{set}-avg{avg}-dt{h}-T{T}-linear-{linear}-eta1_{eta_1}-eta2_{eta_2}-T1_{T1}-T2_{T2}-v0{v0}_M{M}-alpha{alpha}"

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

def msd(shared_msd, shape_msd, lock, process,
        T, h,
        alpha, linear,
        v0, M,
        eta_1, eta_2, T1, T2):   
    temp_msd = to_numpy_array(shared_msd, shape_msd)
    
    eq = fle(alpha, linear)
    eq.params(T = T, h = h,
              v0 = v0, M = M,
              eta_1 = eta_1, eta_2 = eta_2,
              T1 = T1, T2 = T2)
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
        init_msd = np.zeros((avg, n), dtype = np.float32)
    else:
        init_msd = np.zeros(n, dtype = np.float32)
    shared_msd = to_shared_array(init_msd, c.c_float)
    stored_msd = to_numpy_array(shared_msd, init_msd.shape)
    lock = mp.Lock()
    # execute in batches
    for i in tqdm(range(0, avg, batch_size)):
        # execute all tasks in a batch
        processes = [mp.Process(target = msd, 
                                args=(shared_msd, init_msd.shape, lock, p,
                                      T, h, alpha, linear,
                                      v0, M, eta_1, eta_2, T1, T2)) 
                     for p in range(i, i + batch_size)]
        # start all processes
        for process in processes:
            process.start()
        # wait for all processes to complete
        for process in processes:
            process.join()
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
    df_msd = pd.DataFrame(dict(zip(["t", "msd"], [t_, stored_msd/avg])))

ut.write_hdf5(DATA_PATH + f"{task}.hdf5", df_msd)
t = datetime.datetime.now()
print(f"End task {task}: {t.year}/{t.month}/{t.day} {t.hour}:{t.minute}:{t.second}") 
