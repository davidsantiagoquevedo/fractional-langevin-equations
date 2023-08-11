from fbm import FBM
import matplotlib.pyplot as plt
plt.style.use("src/plot_style.mplstyle")
import pandas as pd
import numpy as np
import scipy.special as scp
from tqdm import tqdm
import sys
sys.path.append("src/")
import integration as itg
import mittag_leffler as ml
import fBm_stats as fbs

T = 10
h = 0.01
n = int(T/h)
x_n = [0]*n

v0 = 0
gamma = 1
eta = 1
m = 1
KBT = 1


# Analytical solution 
def solution_fle(t, noise, order, v0):
    t__ = np.array(t)
    noise__ = np.array(noise)
    def nonlinear_term(t):
        z = -gamma*t**(2-order)
        return t * ml.mittag_leffler_vector(z, 2-order, 2)
    conv = itg.convolution(nonlinear_term, noise__, t__)
    nonlinear = v0 * nonlinear_term(t__)
    return nonlinear + conv

def a_j(j):
    return j**(2*H-1) - (j-1)**(2*H-1)

def x_k(k, B_H):
    a_jj = 0
    if k-1 >= 1:    
        for j in range(1, k):
            #print(j, x_n[k-j], x_n[k-j-1])
            a_jj += a_j(j)*(x_n[k-j]+x_n[k-j-1])
    h2Hm = h**(2*H)/(2*(2*H-1)*m)/scp.gamma(2*H-1)
    return x_n[k-1] + h*v0 - h2Hm*a_jj + (eta*h/m)*B_H[k]

def solve(H):
    order = 2 - 2*H
    f = FBM(n = n, hurst = H, length = T, method = 'daviesharte')
    B_H = f.fbm()
    dB_H = n/T*np.diff(B_H)
    t = f.times()
    anl = solution_fle(t[:n], dB_H, order, v0)
    for k in range(1,n):
        x_n[k] = x_k(k,B_H)
    return t, x_n, anl

def plot_results(H, axi, panel, xlabel = False, ylabel = True):
    t, x_n, anl = solve(H)
    axi.plot(t[:n],x_n, ls="", marker = "^", label = "Numerical")
    axi.plot(t[:n],anl, ls = "-", label = "Analytical")
    if ylabel:
        axi.set_ylabel(r"$q(t)$")
    if xlabel:
        axi.set_xlabel(r"$t$")
    axi.set_title(panel)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.3)
    axi.text(0.80, 0.15, f"H = {H}", transform=axi.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
    

def calculate_msd(H, realizations = 100):
    for r in tqdm(range(realizations)):
        t, x_n, anl = solve(H)
        if r == 0:
            df_msd = pd.DataFrame({"t": t[:len(t)-1]})
        df_msd["x_"+str(r)] = x_n        
    df_msd.set_index("t", inplace = True)
    msd = fbs.msd(df_msd, False)
    return msd

def analytical_msd_limit(H, t):
    order = 2 - 2*H
    gm = scp.gamma(1 + order)
    return 2*KBT/(m*eta) * (t**order)/gm

def analytical_msd(H, t):
    order = 2 - 2*H
    z = -eta*t**(2-order)
    return 2*KBT/(m*eta) * (t**2) * ml.mittag_leffler_vector(z, 2-order, 3)
    
    
H = 0.51
T = 20
t = np.arange(0,T,h)
msd = calculate_msd(H, 20)

H = 0.8
anl = analytical_msd(H, t)
anl_lim = analytical_msd_limit(H, t)
fig, ax = plt.subplots(figsize=(10,10))
ax.plot(msd.index, msd, label = "numerical")
ax.plot(t, anl, label = "analytical")
ax.plot(t, anl_lim, ls = "--", label = "analytical limit")
ax.legend()
plt.show()

fig, ax = plt.subplots(2,3, figsize=(15,10), sharex = True)
H = 0.51
axi=ax[0][0]
print(H)
plot_results(H, axi, "a.", xlabel = False, ylabel = True)

H = 0.6
axi=ax[0][1]
print(H)
plot_results(H, axi, "b.", xlabel = False, ylabel = False)

H = 0.7
axi=ax[0][2]
print(H)
plot_results(H, axi, "c.", xlabel = False, ylabel = False)

H = 0.8
axi=ax[1][0]
print(H)
plot_results(H, axi, "d.", xlabel = True, ylabel = True)

H = 0.9
axi=ax[1][1]
print(H)
plot_results(H, axi, "e.", xlabel = True, ylabel = False)

H = 0.99
axi=ax[1][2]
print(H)
plot_results(H, axi, "f.", xlabel = True, ylabel = False)

handles, labels = ax[0][0].get_legend_handles_labels()
lgd = fig.legend(handles, labels, bbox_to_anchor = (0.87, 0.0), ncol = 3)
fig.tight_layout()
fig.savefig("outs/FE.png", dpi=100)

