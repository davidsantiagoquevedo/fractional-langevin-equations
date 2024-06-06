import sys
import matplotlib.pyplot as plt
plt.style.use("analysis/plot_style.mplstyle")
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
import numpy as np

sys.path.append("src/")
from fLe_timecrystal import *
import plot as pt
 
fig, ax = plt.subplots(1,1, figsize = (5,5))

avg = 4000
task_set = ["002"]
data_path = "_raw/time_crystal/"
analytical = False

#Time crystal
# alpha = [0.01, 0.03, 0.05, 0.07, 0.09, 0.1]
# T = 100
# h = 0.005
# v0 = 1.0
# M = 1.0
# eta_1 = 0.0
# eta_2 = 1.0
# T1 = 0.0
# T2 = 1.0
# linear = 0
# title = "Fractional-order and colored noise"

# Time glass
# alpha = [0.01, 0.03, 0.05, 0.07, 0.09, 0.1]
# T = 100
# h = 0.005
# v0 = 1.0
# M = 1.0
# eta_1 = 1.0
# eta_2 = 0.0
# T1 = 1.0
# T2 = 0.0
# linear = 0
# title = "Fractional-order and white noise"

# Mixed phase
# alpha = [0.01, 0.03, 0.05, 0.07, 0.09, 0.1]
# T = 100
# h = 0.005
# v0 = 1.0
# M = 1.0
# eta_1 = 1.0
# eta_2 = 1.0
# T1 = 1.0
# T2 = 1.0
# linear = 0
# title = "Fractional-order, colored and white noise"

# Linear=1
alpha = [0.05, 0.1, 0.3, 0.5, 0.7, 0.9]
T = 100
h = 0.01
v0 = 1.0
M = 1.0
eta_1 = 1.0
eta_2 = 1.0
T1 = 1.0
T2 = 1.0
linear = 1
title = "Linear- and fractional-order dissipation"

prefix = "msd"

task = f"{prefix}-avg{avg}-dt{h}-T{T}-linear{linear}-eta1{eta_1}-eta2{eta_2}-T1{T1}-T2{T2}-v0{v0}-M{M}"

for i, a in enumerate(alpha):
    print(a)
    eq = fle(a, linear)
    eq.params(T = T, h = h,
            v0 = v0, M = M,
            eta_1 = eta_1, eta_2 = eta_2,
            T1 = T1, T2 = T2)
    pt.plot_check(ax, eq, avg, task_set, data_path, trunc = 100, color = colors[i], analytical = analytical)

ax.set_title(title)
fig.tight_layout()
fig.savefig(f"outs/{task}.png", dpi = 100)

ax.set_yscale("log")
ax.set_xscale("log")
fig.savefig(f"outs/log{task}.png", dpi = 100)
