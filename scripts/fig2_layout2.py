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
fig, ax = plt.subplots(2,3, figsize=(15, 10), sharex = False, sharey = False);
avg = 4000
task_set = ["002"]
h_anl =0.1
T_trunc =  20
alpha = [0.01, 0.03, 0.05, 0.07, 0.09, 0.1]

   
def add_axislabel_ins(axins, axi, x_text = "Freq. "+ r"$(1/t)$", y_text = "Amp."):
    props = dict(facecolor='white', edgecolor='white')
    axins.text(0.13, 0.65, x_text, transform=axi.transAxes, bbox=props, fontsize = 10)
    axins.text(0.04, 0.76, y_text, transform=axi.transAxes, bbox=props, fontsize = 10, rotation = 90)
  
#################################################
 ################ TIME CRYSTAL #################
#################################################
axi = ax[0][0]
axins = ax[1][0]
axins_ = axi.inset_axes([0.59, 0.15, 0.37, 0.37])
#axins_ = axi.inset_axes([0.0, 1.17, 1.0, 0.4])

T = 100
h = 0.005
v0 = 1.0
M = 1.0
eta_1 = 0.0
eta_2 = 1.0
T1 = 0.0
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
    plot_msd(eq = eq, avg = avg, task_set = task_set, data_path = data_path, 
             ax = axins, color_fd = colors[9 - i], analytical = True, T = T_trunc, h = h_anl, 
             legend_main = True, legend_second = legend_anl, ci_normal = True, truncate = T_trunc)
    plot_msd(eq = eq, avg = avg, task_set = task_set, data_path = data_path, 
             ax = axi, color_fd = colors[9 - i], analytical = False,
             legend_main = False)
    if i == 0:
        add_trend(axi, x0 = 1, xf = T, func = talpha, text = r"$~t^{\alpha}$", xtext = 3, dy = -0.7, alpha = a)
    else:
        add_trend(axi, x0 = 1, xf = T, func = talpha, alpha = a)
    
    df_ft = get_ft(eq, avg, task_set, data_path)
    plot_fft(axins_, df_ft, color = colors[9 - i], half = True)
    df_freq = get_freq(df_ft)
    
add_1_npi(axins_, n = 2)
pu.resize_names(axins_)

add_trend(axi, x0 = h, xf = T, func = t, text = "~t", xtext = 0.01)
add_trend(axi, x0 = h, xf = 1, func = t2, text = "~t²", xtext = 0.01)
add_trend(axi, x0 = h, xf = 1, func = t3, text = "~t³", xtext = 0.01)
add_perio_grid(axins, 1/np.pi, T_trunc, times = 1)

axi.set_xscale("log")
axi.set_yscale("log")
axins_.set_xscale("log")
axins_.set_yscale("log")

axi.set_ylabel("MSD " r"$\langle x^2 (t) \rangle$")
axi.set_xlabel("Time "+"$t$")

axins.set_ylabel("MSD " r"$\langle x^2 (t) \rangle$")
axins.set_xlabel("Time "+"$t$")

axins_.set_ylabel("FFT - Amplitude ")
axins_.set_xlabel("Freq. "+ r"$(1/t)$")


handles, labels = axins.get_legend_handles_labels()
axins.get_legend().remove()

pu.add_caption_letter(axi, "(a)")
pu.add_caption_letter(axins, "(b)")

#################################################
 ################ TIME GLASS #################
#################################################

axi = ax[0][1]
axins = ax[1][1]
axins_ = axi.inset_axes([0.59, 0.15, 0.37, 0.37])
#axins_ = axi.inset_axes([0.0, 1.17, 1.0, 0.4])

T = 15
h = 0.005
v0 = 1.0
M = 1.0
eta_1 = 1.0
eta_2 = 0.0
T1 = 1.0
T2 = 0.0

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
    plot_msd(eq = eq, avg = avg, task_set = task_set, data_path = data_path, 
             ax = axins, color_fd = colors[9 - i], analytical = True, T = T_trunc, h = h_anl, 
             legend_main = True, legend_second = legend_anl, ci_normal = True, truncate = T_trunc)
    plot_msd(eq = eq, avg = avg, task_set = task_set, data_path = data_path, 
             ax = axi, color_fd = colors[9 - i], analytical = False,
             legend_main = False)
    if i == 0:
        add_trend(axi, x0 = 1, xf = T, func = talpha, text = r"$~t^{\alpha}$", xtext = 3, dy = -0.7, alpha = a)
    else:
        add_trend(axi, x0 = 1, xf = T, func = talpha, alpha = a)
    
    df_ft = get_ft(eq, avg, task_set, data_path)
    plot_fft(axins_, df_ft, color = colors[9 - i], half = True)
    df_freq = get_freq(df_ft)
    
add_1_npi(axins_, n = 1)
pu.resize_names(axins_)

add_trend(axi, x0 = h, xf = T, func = t, text = "~t", xtext = 0.01)
add_trend(axi, x0 = h, xf = 1, func = t2, text = "~t²", xtext = 0.01)
add_trend(axi, x0 = h, xf = 1, func = t3, text = "~t³", xtext = 0.01)
add_perio_grid(axins, 1/np.pi, T_trunc, times = 1)

axi.set_xscale("log")
axi.set_yscale("log")
axins_.set_xscale("log")
axins_.set_yscale("log")

axi.set_ylabel("MSD " r"$\langle x^2 (t) \rangle$")
axi.set_xlabel("Time "+"$t$")

axins.set_ylabel("MSD " r"$\langle x^2 (t) \rangle$")
axins.set_xlabel("Time "+"$t$")

axins_.set_ylabel("FFT - Amplitude ")
axins_.set_xlabel("Freq. "+ r"$(1/t)$")


handles, labels = axins.get_legend_handles_labels()
axins.get_legend().remove()

pu.add_caption_letter(axi, "(c)")
pu.add_caption_letter(axins, "(d)")

#################################################
 ################ MIXED PHASE #################
#################################################

axi = ax[0][2]
axins = ax[1][2]
axins_ = axi.inset_axes([0.59, 0.15, 0.37, 0.37])
#axins_ = axi.inset_axes([0.0, 1.17, 1.0, 0.4])

T = 15
h = 0.005
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
    plot_msd(eq = eq, avg = avg, task_set = task_set, data_path = data_path, 
             ax = axins, color_fd = colors[9 - i], analytical = True, T = T_trunc, h = h_anl, 
             legend_main = True, legend_second = legend_anl, ci_normal = True, truncate = T_trunc)
    plot_msd(eq = eq, avg = avg, task_set = task_set, data_path = data_path, 
             ax = axi, color_fd = colors[9 - i], analytical = False,
             legend_main = False)
    if i == 0:
        add_trend(axi, x0 = 1, xf = T, func = talpha, text = r"$~t^{\alpha}$", xtext = 3, dy = -0.7, alpha = a)
    else:
        add_trend(axi, x0 = 1, xf = T, func = talpha, alpha = a)
    
    df_ft = get_ft(eq, avg, task_set, data_path)
    plot_fft(axins_, df_ft, color = colors[9 - i], half = True)
    df_freq = get_freq(df_ft)
    
add_1_npi(axins_, n = 3)
pu.resize_names(axins_)

add_trend(axi, x0 = h, xf = T, func = t, text = "~t", xtext = 0.01)
add_trend(axi, x0 = h, xf = 1, func = t2, text = "~t²", xtext = 0.01)
add_trend(axi, x0 = h, xf = 1, func = t3, text = "~t³", xtext = 0.01)
add_perio_grid(axins, 1/np.pi, T_trunc, times = 1)

axi.set_xscale("log")
axi.set_yscale("log")
axins_.set_xscale("log")
axins_.set_yscale("log")

axi.set_ylabel("MSD " r"$\langle x^2 (t) \rangle$")
axi.set_xlabel("Time "+"$t$")

axins.set_ylabel("MSD " r"$\langle x^2 (t) \rangle$")
axins.set_xlabel("Time "+"$t$")

axins_.set_ylabel("FFT - Amplitude ")
axins_.set_xlabel("Freq. "+ r"$(1/t)$")


handles, labels = axins.get_legend_handles_labels()
axins.get_legend().remove()

pu.add_caption_letter(axi, "(e)")
pu.add_caption_letter(axins, "(f)")

###### ###### ######
 ###### SAVE ######
###### ###### ######
fig.legend(handles, labels, bbox_to_anchor = (0.95, 1.07), ncol = 8)
fig.tight_layout()

fig.savefig("outs/fig2_layout2_lt.png", dpi = 100)
