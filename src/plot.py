import sys
sys.path.append("/src/")
import utils as ut
import fBm_stats as fbs
from fLe_timecrystal import fle

import numpy as np
import pandas as pd
from scipy import fftpack
from scipy.signal import find_peaks
from tqdm import tqdm
import matplotlib.ticker as tck

def plot_position(eq, 
                  ax,
                  color_fd, color_a = "black",
                  legend_main = False, legend_second = False):
    eq.make_B_H()
    eq.solve()
    eq.analytical_linear() 
    
    if legend_main:
        ax.plot(eq.t, eq.numerical, ls = "-", color = color_fd, label = r"$\alpha$ ="+str(eq.alpha))
    else:
        ax.plot(eq.t, eq.numerical, ls = "-", color = color_fd, label = "")
    if legend_second:
        ax.plot(eq.t, eq.analytical, ls = ":", color = color_a, label = "Analytical")
    else:
        ax.plot(eq.t, eq.analytical, ls = ":", color = color_a, label = "")
                
    ax.legend()
    
def numeric_msd(eq, avg, task_set, data_path, bootstrap = False, percentiles=[0.01,0.99], num_bootstrap_samples = 1000):
    T = eq.T
    h = eq.h
    v0 = eq.v0_in
    M = eq.M
    eta_1 = eq.eta_1
    eta_2 = eq.eta_2
    T1 = eq.T1
    T2 = eq.T2

    alpha = eq.alpha
    linear = eq.linear
    
    for i, tk in enumerate(task_set):
        f = f"trj-set{tk}-avg{avg}-dt{h}-T{T}-linear{linear}-eta1{eta_1}-eta2{eta_2}-T1{T1}-T2{T2}-v0{v0}-M{M}-alpha{alpha}"
        if i == 0:
            df_trj = ut.read_hdf5_data(data_path + f + ".hdf5")
            df_trj = df_trj.set_index("t")
        else:
            df_temp = ut.read_hdf5_data(data_path + f + ".hdf5")
            df_temp = df_temp.set_index("t")
            df_trj = pd.concat([df_trj, df_temp], axis = 1)
    
    msd = fbs.msd(df_trj, False)
    if bootstrap:
        for i in tqdm(range(num_bootstrap_samples)):
            # Randomly sample rows from the DataFrame with replacement
            columns_sampled = np.random.choice(df_trj.columns, size=int(len(df_trj.columns)*0.6), replace=True)
            df_sampled = df_trj[columns_sampled]
            # Compute mean squared displacement for the bootstrap sample
            msd = fbs.msd(df_sampled, normalize = False)
            if i == 0:
                df_msd_bootstrap = pd.DataFrame({f"sample_{i}" : np.array(msd)})   
            else:
                df_temp = pd.DataFrame({f"sample_{i}" : np.array(msd)})
                df_msd_bootstrap = pd.concat([df_msd_bootstrap, df_temp], axis = 1)

        lower = df_msd_bootstrap.quantile(percentiles[0],axis=1)
        upper = df_msd_bootstrap.quantile(percentiles[1],axis=1)
    else:
        lower = [0]*len(np.array(msd))
        upper = [0]*len(np.array(msd))
        
    numeric_msd = pd.DataFrame({"t" : np.array(df_trj.index),
                                "msd" : np.array(msd),
                                "lower" : np.array(lower),
                                "upper" : np.array(upper)})
    
    return numeric_msd

def plot_msd(eq, avg, task_set, data_path,
             ax, color_fd, color_a = "black", bootstrap = False,
             legend_main = False, legend_second = False, **kwargs):
    numeric = numeric_msd(eq, avg, task_set, data_path, bootstrap = bootstrap)
    if "truncate" in kwargs:
        numeric = numeric[numeric.t <= kwargs["truncate"]]
    if legend_main:
        ax.plot(numeric.t, numeric.msd, ls = "-", color = color_fd, label = r"$\alpha$ ="+str(eq.alpha))
        if bootstrap:
            ax.fill_between(numeric.t, numeric.lower, numeric.upper, alpha = .4, color = color_fd)            
    else:
        ax.plot(numeric.t, numeric.msd, ls = "-", color = color_fd, label = "")
    ax.set_xlim(xmin = eq.h)
    #try:
    print("Computing analytical....")
    if kwargs["analytical"]:
        eq_ = fle(eq.alpha, eq.linear)
        eq_.params(T = kwargs["T"], h = kwargs["h"],
        v0 = eq.v0, M = eq.M,
        eta_1 = eq.eta_1, eta_2 = eq.eta_2,
        T1 = eq.T1, T2 = eq.T2)
        eq_.make_B_H()
        if eq_.linear == 1:
            eq_.msd_linear()
        else:
            eq_.msd_non_linear()
        if legend_second:
            ax.plot(eq_.t, eq_.msd, ls = ":", color = color_a, label = "Analytical")
        else:
            ax.plot(eq_.t, eq_.msd, ls = ":", color = color_a, label = "")
        ax.set_xlim(xmin = kwargs["h"])
    #except:
    #    print("No analytical solution")
    ax.legend()
    
def add_trend(ax, x0, xf, func, text = "", xtext = False, dx = 0, dy = 0, **kwargs):
    t = np.arange(x0,xf,0.01)
    if xtext:
        xloc = xtext
    else:
        xloc = xf
    if kwargs:
        ax.plot(t,func(t, alpha = kwargs["alpha"]), color = "black", alpha = 0.5, ls = ":")
        ax.text(xloc + dx, func(xloc, alpha = kwargs["alpha"]) + dy, text)
    else:
        ax.plot(t, func(t), color = "black", alpha = 0.5, ls = ":")
        ax.text(xloc + dx, func(xloc) + dy, text)

def t(t):
    return t

def t2(t):
    return t**2

def t3(t):
    return t**3

def constant(t):
    return 10 -t+t

def talpha(t, **kwargs):
    alpha = kwargs["alpha"]
    return t**alpha

def add_perio_grid(ax, freq, T, times = 8):
    t = 1/freq
    while t <= T:
        ax.axvline(x = t, color = "black", alpha = 0.3, ls = ":")
        t += 1/freq
        
    def numfmt(x, pos): # your custom formatter function: divide by 100.0
        s = '{} $\pi$'.format(int(x / np.pi))
        return s

    yfmt = tck.FuncFormatter(numfmt)
    ax.xaxis.set_major_formatter(yfmt)
    ax.xaxis.set_major_locator(tck.MultipleLocator(base = times*np.pi))
    
    
def get_ft(eq, avg, task_set, data_path):
    numeric = numeric_msd(eq, avg, task_set, data_path)
    h = eq.h
    T = eq.T
    Fs = 1/h #sampling rate
    n = int(T/h) #number of observations
    
    df_fft = pd.DataFrame()
    df_fft["fft"] = fftpack.fft(np.array(numeric.msd))
    df_fft["A"] = df_fft["fft"].abs()
    df_fft["fr"] = Fs/n * np.linspace(0,n,int(n))
    
    return df_fft

def get_freq(ft, q_thr = 0.975, half = True):
    if half:
        N = int(len(ft)/2)
    else:
        N = len(ft)
    x = ft[:N].A
    threshold = x.quantile(q_thr)
    peaks, _ = find_peaks(x, height = threshold)
    
    return ft.iloc[peaks].sort_values(by = "A", ascending = False)

def plot_fft(ax, ft, color, half = True):
    if half:
        N = int(len(ft)/2)
    else:
        N = len(ft)
    ft[:N].set_index("fr")["A"].plot(ax = ax, 
        ls = "-", marker = "", label = "", color = color)
    ax.legend()
    ax.set_ylabel("")
    ax.set_xlabel("")
    
def add_1_npi(ax, n = 1):
    if n == 1:
        label = r"$\frac{1}{\pi}$"
        pos = 0.45
    if n == 2:
        label = r"$\frac{1}{2\pi}$"
        pos = 0.17
    if n == 3:
        n = np.sqrt(2)
        label = r"$\frac{\sqrt{2}}{2\pi}$"
        pos = 0.25
    ax.axvline(x = 1/(n*np.pi), color = "black", alpha = 0.3, ls = ":")
    ax.text(pos, 0.1, label, transform=ax.transAxes)

#### #### #### ##### #### #### ####
#### FUNCTION TO CHECK RESULTS ####
#### #### #### ##### #### #### ####

def plot_check(ax, eq, avg, task_set, data_path, trunc, color):
    numeric = numeric_msd(eq, avg, task_set, data_path, bootstrap = True)
    numeric = numeric[numeric.t <= trunc]
    eq_ = fle(eq.alpha, eq.linear)
    eq_.params(T = trunc, h = 0.1,
            v0 = eq.v0, M = eq.M,
            eta_1 = eq.eta_1, eta_2 = eq.eta_2,
            T1 = eq.T1, T2 = eq.T2)
    eq_.make_B_H()
    if eq.linear == 0:
        eq_.msd_non_linear()
    else:
        eq_.msd_linear()
    if all(np.array((numeric.lower == 0))):
        ax.plot(numeric.t, numeric.msd, label = r"$\alpha$ = "+str(eq.alpha))
    else:
        ax.plot(numeric.t, numeric.msd, label = r"$\alpha$ = "+str(eq.alpha), color = color)
        ax.fill_between(numeric.t, numeric.lower, numeric.upper, alpha = .4, color = color)
    
    ax.plot(eq_.t, eq_.msd, label = "", color = "black", ls = ":")
    ax.set_xlim(xmin = eq_.h)
    ax.legend()
