import sys
import matplotlib.pyplot as plt
from tqdm import tqdm

sys.path.append("src/")

from fLe_timecrystal import fle
import plot_utils as pu
from plot import *

plt.style.use("analysis/plot_style.mplstyle")
data_path = "_raw/time_crystal/"
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
fig, ax = plt.subplots(2,2, figsize=(10, 10), sharex = False);
avg = 4000
task_set = ["002", "003", "004", "005"]

alpha = [0.05, 0.1, 0.3, 0.5, 0.7, 0.9]

##########################################################
print("Plot (a)")
axi = ax[0,0]

T = 15
h = 0.05
v0 = 1.0
M = 1.0
eta_1 = 1.0
eta_2 = 1.0
T1 = 1.0
T2 = 1.0

linear = 1

for i, a in enumerate(alpha):
    print(a)
    eq = fle(a, linear)
    eq.params(T = T, h = h,
            v0 = v0, M = M,
            eta_1 = eta_1, eta_2 = eta_2,
            T1 = T1, T2 = T2)
    plot_position(eq, axi, 
                  color_fd = colors[9 - i])

axi.set_xlabel("Time "+r"$t$")
axi.set_ylabel("Position "+r"$x(t)$")
pu.add_caption_letter(axi, "(a)")

##########################################################
print("Plot (b)")
axi = ax[0,1]

T = 15
h = 0.01
v0 = 1.0
M = 1.0
eta_1 = 1.0
eta_2 = 1.0
T1 = 1.0
T2 = 1.0

linear = 1
for i, a in enumerate(alpha):
    print(a)
    eq = fle(a, linear)
    eq.params(T = T, h = h,
             v0 = v0, M = M,
             eta_1 = eta_1, eta_2 = eta_2,
             T1 = T1, T2 = T2)
    if a == alpha[-1]:
        legend_anl = True
    else: 
        legend_anl = False
    # Main plot: anlytical + numeric + bootstrap
    plot_msd(eq = eq, avg = avg, task_set = task_set, data_path = data_path, 
             ax = axi, color_fd = colors[9 - i], analytical = True, T = 15, h = 0.05,
             legend_main = True, legend_second = legend_anl, bootstrap = True)

axi.set_ylabel("MSD " r"$\langle x^2 (t) \rangle$")
axi.set_xlabel("Time "+"$t$")
axi.get_legend().remove()

pu.add_caption_letter(axi, "(b)")

##########################################################
print("Plot (c)")
axi = ax[1,0]

T = 15
h = 0.01
v0 = 1.0
M = 1.0
eta_1 = 1.0
eta_2 = 1.0
T1 = 1.0
T2 = 1.0

linear = 0
for i, a in enumerate(alpha):
    print(a)
    eq = fle(a, linear)
    eq.params(T = T, h = h,
              v0 = v0, M = M,
              eta_1 = eta_1, eta_2 = eta_2,
              T1 = T1, T2 = T2)
    if a == alpha[-1]:
        legend_anl = True
    else: 
        legend_anl = False
    # Main plot: anlytical + numeric + bootstrap
    plot_msd(eq = eq, avg = avg, task_set = task_set, data_path = data_path, 
             ax = axi, color_fd = colors[9 - i], analytical = True, T = 15, h = 0.05,
             legend_main = True, legend_second = legend_anl, bootstrap = True)


axi.set_ylabel("MSD " r"$\langle x^2 (t) \rangle$")
axi.set_xlabel("Time "+"$t$")

handles, labels = axi.get_legend_handles_labels()
axi.get_legend().remove()
pu.add_caption_letter(axi, "(c)")

##########################################################
print("Plot (d)")
axi = ax[1,1]

T = 15
h = 0.01
v0 = 1.0
M = 1.0
eta_1 = 1.0
eta_2 = 1.0
T1 = 1.0
T2 = 1.0

linear = 0
for i, a in enumerate(alpha):
    print(a)
    eq = fle(a, linear)
    eq.params(T = T, h = h,
              v0 = v0, M = M,
              eta_1 = eta_1, eta_2 = eta_2,
              T1 = T1, T2 = T2)
    plot_msd(eq = eq, avg = avg, task_set = task_set, data_path = data_path, 
             ax = axi, color_fd = colors[9 - i], bootstrap = True, analytical = False)
    if i == 0:
        add_trend(axi, x0 = 1, xf = T, func = talpha, text = r"$~t^{\alpha}$", xtext = 3, dy = -0.7, alpha = a)
    else:
        add_trend(axi, x0 = 1, xf = T, func = talpha, alpha = a)

add_trend(axi, x0 = h, xf = 15, func = t, text = "~t", xtext = 0.1)
add_trend(axi, x0 = h, xf = 1, func = t2, text = "~t²", xtext = 0.1)
add_trend(axi, x0 = h, xf = 1, func = t3, text = "~t³", xtext = 0.1)

axi.set_xscale("log")
axi.set_yscale("log")
axi.set_ylabel("MSD " r"$\langle x^2 (t) \rangle$")
axi.set_xlabel("Time "+"$t$")

axi.get_legend().remove()
pu.add_caption_letter(axi, "(d)")

##########################################################
fig.legend(handles, labels, bbox_to_anchor = (0.85, 1.06), ncol = 4)
fig.tight_layout()

fig.savefig("outs/fig1_.png", dpi = 100)