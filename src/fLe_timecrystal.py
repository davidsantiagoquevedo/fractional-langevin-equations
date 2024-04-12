from fbm import FBM
import numpy as np
from scipy.special import gamma
import pandas as pd
import integration as itg
import mittag_leffler as ml

class fle():
    def __init__(self, alpha, linear = True):
        self.alpha = alpha
        self.H = (2-alpha)/2
        self.linear = linear
    
    def params(self, T, h = 0.01, v0 = 0, M = 1, eta_1 = 1, eta_2 = 1, T1 = 1, T2 = 1):
        """
        Function to set the parameters of the system to solve
        
        Args:
            T (int): Final time. Upper bound of evaluation.
            h (float, optional): Size of time step. Defaults to 0.01.
            v0 (int, optional): Initial velocity of the system. Defaults to 0.
            eta_1 (float, optional): Damping of the first bath (white). Defaults to 1.
            eta_2 (float, optional): Damping of the second bath (colored). Defaults to 1.
            T1 (float, optional): Damping of the first bath (white). Defaults to 1.
            T2 (float, optional): Damping of the second bath (colored). Defaults to 1.
        """
        self.kB = 1
        self.T = T
        self.h = h
        self.n = int(self.T/self.h)
        self.x_n = np.zeros(self.n)
                
        self.M = M
        self.eta_1 = eta_1
        self.eta_2 = eta_2
        self.T1 = T1
        self.T2 = T2
        
        if self.linear:
            kB = self.kB
            self.gamma = eta_1
            self.zeta = eta_2
                        
            self.theta_1 = np.sqrt(2*kB*T2*eta_1)
            self.theta_2 = np.sqrt(kB*T2*eta_2)
        else:
            kB = self.kB
            self.gamma = 0
            self.zeta = eta_1 + eta_2
            alpha = self.alpha
            if eta_1 != 0:
                self.t_alpha = (M/eta_1)^(1/(2-alpha))
            else:
                self.t_alpha = 1
            t_alpha = self.t_alpha
            self.theta_1 = np.sqrt(2*kB*T1*eta_1*t_alpha)
            self.theta_2 = np.sqrt(kB*T2*eta_2)
            
        if v0 == 0: #static
            self.v0 = v0
            self.v02 = 0
        
        elif v0 == 1: #thermal
            self.v0 = 0
            if eta_1 == 0:
                self.v02 = self.kB*T2/M
            if eta_2 == 0:
                self.v02 = self.kB*T1/M
            else:
                self.v02 = self.kB/M*(T1+T2)/2
            
        # Amplitudes for stochastic process
        H = self.H
        self.A = self.theta_1
        self.A_H = self.theta_2/np.sqrt(H*(2*H-1)*gamma(2*H-1))
        
    def external_B(self, B, t):
        n = self.n
        T = self.T
          
        self.B = B
        self.dB = [(B[i+1] - B[i])/(T/n) for i in range(0,n)]
        
        self.t_BH = t
        self.t = self.t_BH[:n]
        
    def make_B_H(self, method = 'daviesharte'):
        n = self.n
        H = self.H
        T = self.T
        A = self.A
        A_H = self.A_H
        np.random.seed()
          
        f = FBM(n = n, hurst = H, length = T, method = method)
        self.B_H = f.fbm()
        self.dB_H = n/T*np.diff(self.B_H)
               
        f12 = FBM(n = n, hurst = 0.5, length = T, method = method)
        self.B_12 = f12.fbm()
        self.dB_12 = n/T*np.diff(self.B_12)
        
        self.t_BH = f.times()
        self.t = self.t_BH[:self.n]
        self.B = A*self.B_12 + A_H*self.B_H
        self.dB = A*self.dB_12 + A_H*self.dB_H
            
    # Numerical solution
    def a_j(self, j):
        H = self.H
        return j**(2*H-1) - (j-1)**(2*H-1)
    
    def x_k(self, k):
        H = self.H
        B = self.B
        
        h = self.h
        x_n = self.x_n
        
        v0 = self.v0
        M = self.M
        zeta = self.zeta
        gmma = self.gamma
         
        a_jj = 0
        if k-1 >= 1:    
            for j in range(1, k):
                a_jj += self.a_j(j)*(x_n[k-j]+x_n[k-j-1])
        h2Hm = h**(2*H)/(2*(2*H-1)*gamma(2*H-1))
        return (M/(M+h*gmma))*x_n[k-1] + (h*M/(M+h*gmma))*v0 + (h/(M+h*gmma))*B[k] - (zeta/(M+h*gmma))*h2Hm*a_jj

    def solve(self):
        n = self.n
        x_n = self.x_n
        
        for k in range(1,n):
            x_n[k] = self.x_k(k)
        
        self.numerical = x_n
        
        
    def get_analytical(self):
        alpha = self.alpha
        
        M = self.M
        gmma = self.gamma
        zeta = self.zeta
        
        t = self.t
        noise = self.dB
        v0 =self.v0
                
        t__ = np.array(t)
        noise__ = np.array(noise)
                
        def relaxation_function(t):
            inf = 40
            z = -(gmma/M)*t
            G = pd.DataFrame()
            for n in range(inf):
                t_ = ((-zeta/M)**n)*t**((2-alpha)*n+1)
                G[f"n{n}"] = ml.prabhakar_mittag_leffler(z, 1 , (2-alpha)* n + 2, n+1) * t_ 
            
            return np.array(G.sum(axis = 1))
        
        conv = itg.convolution(relaxation_function, noise__, t__)
        linear_v0 = v0 * relaxation_function(t__)
        self.analytical = linear_v0 + (1/M) * conv
        
    def get_analytical_msd(self):
        alpha = self.alpha
        t = self.t
        
        zeta = self.zeta
        gmma = self.gamma
        M = self.M
        T1 = self.T1
        T2 = self.T2
        assert(T1==T2)
        kBT = self.kB*T2
                             
        def msd(t):
            z = -(gmma/M)*t
            G = pd.DataFrame()
            inf = 70
            for n in range(inf):
                t_ = (t**((2-alpha)*n+2)) * ((-zeta/M)**n)
                G[f"n{n}"] = ml.prabhakar_mittag_leffler(z, 1, (2-alpha)*n + 3, n+1) * t_
            return np.array(G.sum(axis = 1))
        
        self.analytical_msd = 2*kBT/M*msd(t)
        
    def relaxation_non_linear(self):
        alpha = self.alpha
        zeta = self.zeta
        M = self.M
        
        t = self.t
        dB = self.dB
        
        def relaxation_function(t):
            z = (-zeta / M) * t**(2 - alpha)
            return t * ml.mittag_leffler(z, 2 - alpha, 2)
        
        self.G = relaxation_function(t)
        self.G_conv_dB = itg.convolution(relaxation_function, dB, t)
    
    def get_analytical_colored(self):
        assert(self.eta_1 == 0 and self.linear == False)
        
        M = self.M
        v0 = self.v0
        
        self.relaxation_non_linear()
        self.analytical = v0*self.G +  (1/M) * self.G_conv_dB
        
    def get_analytical_colored_msd(self):
        assert(self.eta_1 == 0 and self.linear == False)
        
        alpha = self.alpha
        
        zeta = self.zeta
        M = self.M
        T2 = self.T2
        v02 = self.v02
        kBT = self.kB*T2
        
        t = self.t
        self.relaxation_non_linear()
        G2 = self.G*self.G
        z = -(zeta/M)*t**(2-alpha)
        self.analytical_colored_msd = v02 * G2 + 2*kBT/(M) * ((t**2) * ml.mittag_leffler(z, 2-alpha, 3) - 1/2*G2)