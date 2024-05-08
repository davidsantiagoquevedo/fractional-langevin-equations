import sys
import matplotlib.pyplot as plt
plt.style.use("analysis/plot_style.mplstyle")
import numpy as np

sys.path.append("src/")
from fLe_timecrystal import *
import plot as pt


def plot(ax, eq, avg, task_set, data_path, trunc):
    numeric = pt.numeric_msd(eq, avg, task_set, data_path, bootstrap = True)
    numeric = numeric[numeric.t <= trunc]
    eq_ = fle(eq.alpha, eq.linear)
    eq_.params(T = trunc, h = 0.1,
            v0 = eq.v0, M = eq.M,
            eta_1 = eq.eta_1, eta_2 = eq.eta_2,
            T1 = eq.T1, T2 = eq.T2)
    eq_.make_B_H()
    eq_.msd_non_linear()
    if all(np.array((numeric.lower == 0))) :
        ax.plot(numeric.t, numeric.msd, label = r"$\alpha$ = "+str(eq.alpha))
    else:
        ax.fill_between(numeric.t, numeric.lower, numeric.upper, label = r"$\alpha$ = "+str(eq.alpha), alpha=.4)
    
    ax.plot(eq_.t, eq_.msd, label = "", color = "black", ls = ":")
    ax.set_xlim(xmin = eq.h)
    ax.legend()
    
fig, ax = plt.subplots(1,3, figsize = (15,5))
alpha = [0.05, 0.08, 0.4, 0.8]#, 0.3, 0.5, 0.7]

avg = 4000
task_set = ["001"]
data_path = "_raw/time_crystal/"

#Time crystal
T = 15
h = 0.005
v0 = 1.0
M = 1.0
eta_1 = 0.0
eta_2 = 1.0
T1 = 0.0
T2 = 1.0
linear = 0
axi = ax[0]
for a in alpha:
    print(a)
    eq = fle(a, linear)
    eq.params(T = T, h = h,
            v0 = v0, M = M,
            eta_1 = eta_1, eta_2 = eta_2,
            T1 = T1, T2 = T2)
    plot(axi, eq, avg, task_set, data_path, trunc = 15)
    
#Time glass
T = 15
h = 0.005
v0 = 0.0
M = 1.0
eta_1 = 1.0
eta_2 = 0.0
T1 = 1.0
T2 = 0.0
linear = 0
axi = ax[1]
for a in alpha:
    print(a)
    eq = fle(a, linear)
    eq.params(T = T, h = h,
            v0 = v0, M = M,
            eta_1 = eta_1, eta_2 = eta_2,
            T1 = T1, T2 = T2)
    plot(axi, eq, avg, task_set, data_path, trunc = 15)

alpha = [0.05, 0.08, 0.4, 0.8]
    
#Mixed phase
T = 15
h = 0.005
v0 = 1.0
M = 1.0
eta_1 = 1.0
eta_2 = 1.0
T1 = 1.0
T2 = 1.0
linear = 0
axi = ax[2]
for a in alpha:
    print(a)
    eq = fle(a, linear)
    eq.params(T = T, h = h,
            v0 = v0, M = M,
            eta_1 = eta_1, eta_2 = eta_2,
            T1 = T1, T2 = T2)
    plot(axi, eq, avg, task_set, data_path, trunc = 15)
    
fig.tight_layout()
fig.savefig("outs/fig2_check.png", dpi = 200)
