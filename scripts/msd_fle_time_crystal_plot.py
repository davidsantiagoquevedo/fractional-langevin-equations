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
task_set = ["001"]
data_path = "_raw/time_crystal/"

#Time crystal
#alpha = [0.05, 0.08, 0.4, 0.8]
#T = 15
#h = 0.005
#v0 = 1.0
#M = 1.0
#eta_1 = 0.0
#eta_2 = 1.0
#T1 = 0.0
#T2 = 1.0
#linear = 0

#Time glass
#alpha = [0.01, 0.03, 0.05, 0.08, 0.4, 0.8]
#T = 15
#h = 0.005
#v0 = 1.0
#M = 1.0
#eta_1 = 1.0
#eta_2 = 0.0
#T1 = 1.0
#T2 = 0.0
#linear = 0

#Mixed phase
alpha = [0.05, 0.08, 0.4, 0.8]
T = 15
h = 0.005
v0 = 1.0
M = 1.0
eta_1 = 1.0
eta_2 = 1.0
T1 = 1.0
T2 = 1.0
linear = 0

#Linear=1
#alpha = [0.05, 0.1, 0.3, 0.5, 0.7, 0.9]
#T = 100
#h = 0.05
#v0 = 1.0
#M = 1.0
#eta_1 = 1.0
#eta_2 = 1.0
#T1 = 1.0
#T2 = 1.0
#linear = 1

#task_set = ["001"]

prefix = "msd"

task = f"{prefix}-avg{avg}-dt{h}-T{T}-linear{linear}-eta1{eta_1}-eta2{eta_2}-T1{T1}-T2{T2}-v0{v0}-M{M}"

for i, a in enumerate(alpha):
    print(a)
    eq = fle(a, linear)
    eq.params(T = T, h = h,
            v0 = v0, M = M,
            eta_1 = eta_1, eta_2 = eta_2,
            T1 = T1, T2 = T2)
    pt.plot_check(ax, eq, avg, task_set, data_path, trunc = 15, color = colors[i])

fig.tight_layout()
fig.savefig(f"outs/{task}.png", dpi = 100)
ax.set_yscale("log")
ax.set_xscale("log")
fig.savefig(f"outs/log{task}.png", dpi = 100)
