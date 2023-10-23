import sys
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
from tqdm import tqdm

sys.path.append("src/")
from fLe_twobath import fle_twobath
import fBm_stats as fbs
import utils as ut
import datetime

DATA_PATH = "data/two_baths/"

H = float(sys.argv[1])
A = float(sys.argv[2])
eta = float(sys.argv[3])
C = float(sys.argv[4])
theta_12 = float(sys.argv[5])
theta_H = float(sys.argv[6])

t = datetime.datetime.now()
print(f"Start task msd-h-{H}-A{A}-eta{eta}-C{C}-t12_{theta_12}-tH{theta_H}: {t.year}/{t.month}/{t.day} {t.hour}:{t.minute}:{t.second}")

T = 20
h = 0.005
realizations = 20000

def msd(H, T, A, eta, C, theta_12, theta_H,
        realizations = 100, h = 0.01):
    for r in tqdm(range(realizations)):
        eq = fle_twobath(H)
        eq.params(T, h, v0 = 0, 
                  A = A, eta = eta, C = C, 
                  theta_12 = theta_12, theta_H = theta_H)
        eq.make_B_H()
        eq.solve()        
        if r == 0:
            df_msd = pd.DataFrame({"t": eq.t})
        df_msd["x_"+str(r)] = eq.numerical        
    df_msd.set_index("t", inplace = True)
    msd = fbs.msd(df_msd, False)
    return msd

df_msd = pd.DataFrame(msd(H = H, T = T, 
                          A = A, eta = eta, C = C, 
                          theta_12 = theta_12, theta_H = theta_H,
                          realizations = realizations, h = h), columns=["msd"])

ut.write_hdf5(DATA_PATH + f"msd-h-{H}-A{A}-eta{eta}-C{C}-t12_{theta_12}-tH{theta_H}.hdf5", df_msd)
t = datetime.datetime.now()
print(f"End task msd-h-{H}-A{A}-eta{eta}-C{C}-t12_{theta_12}-tH{theta_H}: {t.year}/{t.month}/{t.day} {t.hour}:{t.minute}:{t.second}")