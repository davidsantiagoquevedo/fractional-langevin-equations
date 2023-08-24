from fbm import FBM
import numpy as np
from scipy.special import gamma
import integration as itg
import mittag_leffler as ml


class fle():
    def __init__(self, H):
        self.H = H
    
    def params(self, T, h = 0.01, v0 = 1, eta = 1, m = 1, KBT = 1):
        self.T = T
        self.h = h
        self.n = int(self.T/self.h)
        self.x_n = [0]*self.n

        self.v0 = v0
        self.eta = eta
        self.m = m
        self.KBT = KBT
        
        H = self.H
        self.zeta = np.sqrt((2*m*KBT*eta) * gamma(1.5 - H) * gamma(0.5 + H) / (gamma(2*H - 1) * gamma(2 - 2*H)))
        #self.zeta = 1
        
    def external_B_H(self, B_H, t):
        n = self.n
        T = self.T
          
        self.B_H = B_H
        self.dB_H = [T*(B_H[i+1] - B_H[i])/n for i in range(0,n)]
        
        self.t_BH = t
        self.t = self.t_BH[:n]
        
    def make_B_H(self, method = 'daviesharte'):
        n = self.n
        H = self.H
        T = self.T
          
        f = FBM(n = n, hurst = H, length = T, method = method)
        self.B_H = f.fbm()
        self.dB_H = n/T*np.diff(self.B_H)
        
        self.t_BH = f.times()
        self.t = self.t_BH[:self.n]
         
    # Analytical solution 
    def get_analytical(self):
        H = self.H
        
        eta = self.eta
        v0 = self.v0
        m = self.m
        zeta = self.zeta
        
        t = self.t
        noise = self.dB_H
        
        order = 2 - 2*H
        t__ = np.array(t)
        noise__ = np.array(noise)
        
        def nonlinear_term(t):
            z = -eta*t**(2-order)
            return t * ml.mittag_leffler_vector(z, 2-order, 2)
        
        conv = itg.convolution(nonlinear_term, noise__, t__)
        
        nonlinear = v0 * nonlinear_term(t__)
        
        self.analytical = nonlinear + (zeta/m) * conv
    
    # Numerical solution
    def a_j(self, j):
        H = self.H
        return j**(2*H-1) - (j-1)**(2*H-1)
    
    def x_k(self, k):
        H = self.H
        B_H = self.B_H
        
        h = self.h
        x_n = self.x_n
        
        v0 = self.v0
        eta = self.eta
        m = self.m
        zeta = self.zeta
         
        a_jj = 0
        if k-1 >= 1:    
            for j in range(1, k):
                #print(j, x_n[k-j], x_n[k-j-1])
                a_jj += self.a_j(j)*(x_n[k-j]+x_n[k-j-1])
        h2Hm = h**(2*H)/(2*(2*H-1)*m)/gamma(2*H-1)
        #h2Hm = h**(2*H)/((2*H-1)*m)
        return x_n[k-1] + h*v0 - eta*h2Hm*a_jj + (zeta*h/m)*B_H[k]

    def solve(self):
        n = self.n
        x_n = self.x_n
        
        for k in range(1,n):
            x_n[k] = self.x_k(k)
        
        self.numerical = x_n
        
    # Mean squarred displacement
    def get_msd_analytical(self):
        H = self.H
        
        eta = self.eta
        m = self.m
        KBT = self.KBT
        t = self.t
        
        order = 2 - 2*H
        #*gamma(1-order)
        z = -eta*t**(2-order)
        self.msd_analytical = 2*KBT/(m) * (t**2) * ml.mittag_leffler_vector(z, 2-order, 3)
        
    def get_msd_analytical_limit(self):
        H = self.H
        
        eta = self.eta
        m = self.m
        KBT = self.KBT
        t = self.t
        
        order = 2 - 2*H
        #*gamma(1-order)
        self.msd_analytical_limit = 2*KBT/(m*eta) * (t**order)/(gamma(1 + order))