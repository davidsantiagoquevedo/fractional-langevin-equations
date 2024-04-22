import sys
sys.path.append("/src/")
import utils as ut
import fBm_stats as fbs
from fLe_timecrystal import fle

import numpy as np
import pandas as pd
from scipy import fftpack
from scipy.signal import find_peaks
import matplotlib.ticker as tck

def plot_position(eq, 
                  ax,
                  color_fd, color_a = "black",
                  legend_main = False, legend_second = False):
    eq.make_B_H()
    eq.solve()
    eq.get_analytical() 
    
    if legend_main:
        ax.plot(eq.t, eq.numerical, ls = "-", color = color_fd, label = r"$\alpha$ ="+str(eq.alpha))
    else:
        ax.plot(eq.t, eq.numerical, ls = "-", color = color_fd, label = "")
    if legend_second:
        ax.plot(eq.t, eq.analytical, ls = ":", color = color_a, label = "Analytical")
    else:
        ax.plot(eq.t, eq.analytical, ls = ":", color = color_a, label = "")
                
    ax.legend()
    
def numeric_msd(eq, avg, task_set, data_path):
    T = eq.T
    h = eq.h
    v0 = eq.v0
    M = eq.M
    eta_1 = eq.eta_1
    eta_2 = eq.eta_2
    T1 = eq.T1
    T2 = eq.T2

    alpha = eq.alpha
    linear = eq.linear
    
    for i, tk in enumerate(task_set):
        f = f"trj-set{tk}-avg{avg}-dt{h}-T{T}-linear{linear}-eta1{eta_1}-eta2{eta_2}-T1{T1}-T2{T2}-v0{v0}_M{M}-alpha{alpha}"
        if i == 0:
            df_trj = ut.read_hdf5_data(data_path + f + ".hdf5")
            df_trj = df_trj.set_index("t")
        else:
            df_temp = ut.read_hdf5_data(data_path + f + ".hdf5")
            df_temp = df_temp.set_index("t")
            df_trj = pd.concat([df_trj, df_temp], axis = 1)
    msd = fbs.msd(df_trj, False).reset_index()
    msd.columns = ["t", "msd"]
    return msd

def plot_msd(eq, avg, task_set, data_path,
             ax, color_fd, color_a = "black",
             legend_main = False, legend_second = False, **kwargs):
    numeric = numeric_msd(eq, avg, task_set, data_path)
    if "truncate" in kwargs:
        numeric = numeric[numeric.t <= kwargs["truncate"]]
    if legend_main:
        numeric.set_index("t").msd.plot(ax = ax, ls = "-", color = color_fd, label = r"$\alpha$ ="+str(eq.alpha))
    else:
        numeric.set_index("t").msd.plot(ax = ax, ls = "-", color = color_fd, label = "")
    try:
        print("Computing analytical....")
        if kwargs["analytical"]:
            eq_ = fle(eq.alpha, eq.linear)
            eq_.params(T = kwargs["T"], h = kwargs["h"],
            v0 = 0, M = eq.M,
            eta_1 = eq.eta_1, eta_2 = eq.eta_2,
            T1 = eq.T1, T2 = eq.T2)
            eq_.make_B_H()
            eq_.get_analytical_msd()
            if legend_second:
                ax.plot(eq_.t, eq_.analytical_msd, ls = ":", color = color_a, label = "Analytical")
            else:
                ax.plot(eq_.t, eq_.analytical_msd, ls = ":", color = color_a, label = "")
    except:
        print("No analytical solution")
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

def plot_fft(ax, ft, half = True):
    if half:
        N = int(len(ft)/2)
    else:
        N = len(ft)
    ft[:N].set_index("fr")["A"].plot(ax = ax, 
        ls = "-", marker = "", label = "")
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
        label = r"$\frac{\sqrt(2)}{2\pi}$"
        pos = 0.25
    ax.axvline(x = 1/(n*np.pi), color = "black", alpha = 0.3, ls = ":")
    ax.text(pos, 0.1, label, transform=ax.transAxes)
    