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
fig, ax = plt.subplots(1,3, figsize=(15, 5), sharex = False, sharey = False);
avg = 4000
task_set = ["002", "003", "004", "005"]
h_anl =0.1

alpha = [0.01, 0.03, 0.05, 0.07, 0.09, 0.1]

   
def add_axislabel_ins(axins, axi, x_text = "Freq. "+ r"$(1/t)$", y_text = "Amp."):
    props = dict(facecolor='white', edgecolor='white')
    axins.text(0.13, 0.65, x_text, transform=axi.transAxes, bbox=props, fontsize = 10)
    axins.text(0.04, 0.76, y_text, transform=axi.transAxes, bbox=props, fontsize = 10, rotation = 90)
  
#################################################
 ################ TIME CRYSTAL #################
#################################################
axi = ax[0]
axins = axi.inset_axes([0.6, 1.17, 0.4, 0.4])
axins_ = axi.inset_axes([0.0, 1.17, 0.4, 0.4])

T = 15
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
             ax = axi, color_fd = colors[9 - i], analytical = True, T = 15, h = h_anl, 
             legend_main = True, legend_second = legend_anl, ci_normal = True)
    plot_msd(eq = eq, avg = avg, task_set = task_set, data_path = data_path, 
             ax = axins, color_fd = colors[9 - i], analytical = False,
             legend_main = False)
    
    df_ft = get_ft(eq, avg, task_set, data_path)
    plot_fft(axins_, df_ft, color = colors[9 - i], half = True)
    df_freq = get_freq(df_ft)
    
add_1_npi(axins_, n = 2)
pu.resize_names(axins)
pu.resize_names(axins_)

add_trend(axins, x0 = h, xf = 15, func = t, text = "~t", xtext = 0.01, fontsize = 11)
add_trend(axins, x0 = h, xf = 1, func = t2, text = "~t²", xtext = 0.01, fontsize = 11)
add_trend(axins, x0 = h, xf = 1, func = t3, text = "~t³", xtext = 0.01, fontsize = 11)
add_perio_grid(axi, 1/np.pi, 15, times = 1)

axins.set_xscale("log")
axins.set_yscale("log")
axins_.set_xscale("log")
axins_.set_yscale("log")

axi.set_ylabel("MSD " r"$\langle x^2 (t) \rangle$")
axi.set_xlabel("Time "+"$t$")
axins.set_ylabel("MSD " r"$\langle x^2 (t) \rangle$")
axins.set_xlabel("Time "+"$t$")
axins_.set_ylabel("FFT - Amplitude ")
axins_.set_xlabel("Freq. "+ r"$(1/t)$")


handles, labels = axi.get_legend_handles_labels()
axi.get_legend().remove()

pu.add_caption_letter(axi, "(a)")
pu.add_caption_letter(axins_, "(a-i)", fontsize = 12)
pu.add_caption_letter(axins, "(a-ii)", fontsize = 12)

#################################################
 ################ TIME GLASS #################
#################################################

axi = ax[1]
axins = axi.inset_axes([0.6, 1.17, 0.4, 0.4])
axins_ = axi.inset_axes([0.0, 1.17, 0.4, 0.4])

T = 15
h = 0.005
v0 = 0.0
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
             ax = axi, color_fd = colors[9 - i], analytical = True, T = 15, h = h_anl, 
             legend_main = True, legend_second = legend_anl, ci_normal = True)
    plot_msd(eq = eq, avg = avg, task_set = task_set, data_path = data_path, 
             ax = axins, color_fd = colors[9 - i], analytical = False,
             legend_main = False)
    
    df_ft = get_ft(eq, avg, task_set, data_path)
    plot_fft(axins_, df_ft, color = colors[9 - i], half = True)
    df_freq = get_freq(df_ft)
    
add_1_npi(axins_, n = 1)
pu.resize_names(axins)
pu.resize_names(axins_)

add_trend(axins, x0 = h, xf = 15, func = t, text = "~t", xtext = 0.01, fontsize = 11)
add_trend(axins, x0 = h, xf = 1, func = t2, text = "~t²", xtext = 0.01, fontsize = 11)
add_trend(axins, x0 = h, xf = 1, func = t3, text = "~t³", xtext = 0.01, fontsize = 11)
add_perio_grid(axi, 1/np.pi, 15, times = 1)

axins.set_xscale("log")
axins.set_yscale("log")
axins_.set_xscale("log")
axins_.set_yscale("log")

axi.set_ylabel("MSD " r"$\langle x^2 (t) \rangle$")
axi.set_xlabel("Time "+"$t$")
axins.set_ylabel("MSD " r"$\langle x^2 (t) \rangle$")
axins.set_xlabel("Time "+"$t$")
axins_.set_ylabel("FFT - Amplitude ")
axins_.set_xlabel("Freq. "+ r"$(1/t)$")


handles, labels = axi.get_legend_handles_labels()
axi.get_legend().remove()

pu.add_caption_letter(axi, "(b)")
pu.add_caption_letter(axins_, "(b-i)", fontsize = 12)
pu.add_caption_letter(axins, "(b-ii)", fontsize = 12)

#################################################
 ################ MIXED PHASE #################
#################################################

axi = ax[2]
axins = axi.inset_axes([0.6, 1.17, 0.4, 0.4])
axins_ = axi.inset_axes([0.0, 1.17, 0.4, 0.4])

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
             ax = axi, color_fd = colors[9 - i], analytical = True, T = 15, h = h_anl, 
             legend_main = True, legend_second = legend_anl, ci_normal = True)
    plot_msd(eq = eq, avg = avg, task_set = task_set, data_path = data_path, 
             ax = axins, color_fd = colors[9 - i], analytical = False,
             legend_main = False)
    
    df_ft = get_ft(eq, avg, task_set, data_path)
    plot_fft(axins_, df_ft, color = colors[9 - i], half = True)
    df_freq = get_freq(df_ft)
    
add_1_npi(axins_, n = 3)
pu.resize_names(axins)
pu.resize_names(axins_)

add_trend(axins, x0 = h, xf = 15, func = t, text = "~t", xtext = 0.01, fontsize = 11)
add_trend(axins, x0 = h, xf = 1, func = t2, text = "~t²", xtext = 0.01, fontsize = 11)
add_trend(axins, x0 = h, xf = 1, func = t3, text = "~t³", xtext = 0.01, fontsize = 11)
add_perio_grid(axi, 1/np.pi, 15, times = 1)

axins.set_xscale("log")
axins.set_yscale("log")
axins_.set_xscale("log")
axins_.set_yscale("log")

axi.set_ylabel("MSD " r"$\langle x^2 (t) \rangle$")
axi.set_xlabel("Time "+"$t$")
axins.set_ylabel("MSD " r"$\langle x^2 (t) \rangle$")
axins.set_xlabel("Time "+"$t$")
axins_.set_ylabel("FFT - Amplitude ")
axins_.set_xlabel("Freq. "+ r"$(1/t)$")


handles, labels = axi.get_legend_handles_labels()
axi.get_legend().remove()

pu.add_caption_letter(axi, "(c)")
pu.add_caption_letter(axins_, "(c-i)", fontsize = 12)
pu.add_caption_letter(axins, "(c-ii)", fontsize = 12)

###### ###### ######
 ###### SAVE ######
###### ###### ######
fig.legend(handles, labels, bbox_to_anchor = (0.95, 1.2), ncol = 8)
fig.tight_layout()

fig.savefig("outs/fig2_.png", dpi = 100)
