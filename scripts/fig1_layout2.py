import sys
import matplotlib.pyplot as plt

sys.path.append("src/")

from fLe_timecrystal import fle
import plot_utils as pu
from plot import *

plt.style.use("analysis/plot_style.mplstyle")
data_path = "_raw/time_crystal/"
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
fig, ax = plt.subplots(1,3, figsize=(15, 5), sharex = False);
avg = 4000
h_anl = 0.1
task_set = ["002"]

alpha = [0.05, 0.1, 0.3, 0.5, 0.7, 0.9]

##########################################################
print("Plot (a)")
axi = ax[0]

T = 15
h = 0.1
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
axi = ax[1]

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
    # Main plot: anlytical + numeric + ci_normal
    plot_msd(eq = eq, avg = avg, task_set = task_set, data_path = data_path, 
             ax = axi, color_fd = colors[9 - i], analytical = True, T = 15, h = h_anl,
             legend_main = True, legend_second = legend_anl, ci_normal = True)

axi.set_ylabel("MSD " r"$\langle x^2 (t) \rangle$")
axi.set_xlabel("Time "+"$t$")
axi.get_legend().remove()

pu.add_caption_letter(axi, "(b)")

##########################################################
print("Plot (c)")
axi = ax[2]

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
    # Main plot: anlytical + numeric + ci_normal
    plot_msd(eq = eq, avg = avg, task_set = task_set, data_path = data_path, 
             ax = axi, color_fd = colors[9 - i], analytical = True, T = 15, h = h_anl,
             legend_main = True, legend_second = legend_anl, ci_normal = True)


axi.set_ylabel("MSD " r"$\langle x^2 (t) \rangle$")
axi.set_xlabel("Time "+"$t$")

handles, labels = axi.get_legend_handles_labels()
axi.get_legend().remove()
pu.add_caption_letter(axi, "(c)")

##########################################################
print("Plot (d)")
axins = axi.inset_axes([0.10, 0.53, 0.37, 0.37])

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
             ax = axins, color_fd = colors[9 - i], ci_normal = True, analytical = False)
    if i == 0:
        add_trend(axins, x0 = 1, xf = T, func = talpha, text = r"$~t^{\alpha}$", xtext = 3, dy = -0.93, alpha = a)
    else:
        add_trend(axins, x0 = 1, xf = T, func = talpha, alpha = a)

add_trend(axins, x0 = h, xf = 15, func = t, text = "~t", xtext = 0.1)
add_trend(axins, x0 = h, xf = 1, func = t2, text = "~t²", xtext = 0.1)
add_trend(axins, x0 = h, xf = 1, func = t3, text = "~t³", xtext = 0.1)

axins.set_xscale("log")
axins.set_yscale("log")

axins.get_legend().remove()
pu.resize_names(axins)

##########################################################
fig.legend(handles, labels, bbox_to_anchor = (0.95, 1.03), ncol = 8)
fig.tight_layout()

fig.savefig("outs/fig1_layout2.png", dpi = 100)
