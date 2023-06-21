from fbm import FBM
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as scp
import tqdm as tqdm
import sys
sys.path.append("src/")
import integration as itg
import mittag_leffler as ml

T = 10
H = 0.7
h = 0.01
n = int(T/h)
x_n = [0]*n

order = 2 - 2*H
v0 = 0
gamma = 1
eta = 1
m = 1


f = FBM(n = n, hurst = H, length = T, method = 'daviesharte')
B_H = f.fbm()
dB_H = n/T*np.diff(B_H)
t = f.times()

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

def x_k(k):
    a_jj = 0
    if k-1 >= 1:    
        for j in range(1, k):
            print(j, x_n[k-j], x_n[k-j-1])
            a_jj += a_j(j)*(x_n[k-j]+x_n[k-j-1])
    h2Hm = h**(2*H)/(2*(2*H-1)*m)/scp.gamma(2*H-1)
    return x_n[k-1] + h*v0 - h2Hm*a_jj + (eta*h/m)*B_H[k]

anl = solution_fle(t[:n], dB_H, order, v0)
for k in range(1,n):
    x_n[k] = x_k(k)

fig, ax = plt.subplots()
ax.plot(t[:n],x_n, label = "Numeric")
ax.plot(t[:n],anl, ls = ":", label = "Analytical")
ax.legend()
plt.show()


