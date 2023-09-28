import sys
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
from tqdm import tqdm

sys.path.append("src/")
from fLe import fle
import fBm_stats as fbs
import utils as ut

DATA_PATH = "data/fle/"

H = float(sys.argv[1])
T = int(sys.argv[2])
h = float(sys.argv[3])
realizations = int(sys.argv[4])

def msd(H, T, realizations = 100, h = 0.01):
    for r in tqdm(range(realizations)):
        eq = fle(H)
        #zeta = np.sqrt(3 - 2*H)
        #zeta = np.sqrt(2-(H*(2*H-1)))
        zeta = np.sqrt((2*H-1)/(H*(2*H-1)))
        eq.params(T, h, zeta = zeta)
        eq.make_B_H()
        eq.solve()        
        if r == 0:
            df_msd = pd.DataFrame({"t": eq.t})
        df_msd["x_"+str(r)] = eq.numerical        
    df_msd.set_index("t", inplace = True)
    msd = fbs.msd(df_msd, False)
    return msd

df_msd = pd.DataFrame(msd(H = H, T = T, realizations = realizations, h = h), columns=["msd"])
ut.write_hdf5(DATA_PATH + f"msd-h-{H}-{realizations}.hdf5", df_msd)
