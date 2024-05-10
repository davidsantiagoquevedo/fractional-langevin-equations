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
#T = 15
#h = 0.005
#v0 = 1.0
#M = 1.0
#eta_1 = 0.0
#eta_2 = 1.0
#T1 = 0.0
#T2 = 1.0
#linear = 0
#alpha = [0.05, 0.08, 0.4, 0.8]

#Time glass
alpha = [0.01, 0.03, 0.05, 0.08, 0.4, 0.8]
T = 15
h = 0.005
v0 = 1.0
M = 1.0
eta_1 = 1.0
eta_2 = 0.0
T1 = 1.0
T2 = 0.0
linear = 0

#Mixed phase
#T = 15
#h = 0.005
#v0 = 1.0
#M = 1.0
#eta_1 = 1.0
#eta_2 = 1.0
#T1 = 1.0
#T2 = 1.0
#linear = 0

prefix = "msd"

task = f"{prefix}-avg{avg}-dt{h}-T{T}-linear{linear}-eta1{eta_1}-eta2{eta_2}-T1{T1}-T2{T2}-v0{v0}-M{M}"

axi = ax
for i, a in enumerate(alpha):
    print(a)
    eq = fle(a, linear)
    eq.params(T = T, h = h,
            v0 = v0, M = M,
            eta_1 = eta_1, eta_2 = eta_2,
            T1 = T1, T2 = T2)
    pt.plot_check(axi, eq, avg, task_set, data_path, trunc = 15, color = colors[i])

fig.tight_layout()
fig.savefig(f"outs/{task}.png", dpi = 100)
