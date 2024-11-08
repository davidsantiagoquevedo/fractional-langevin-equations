from fbm import FBM
import numpy as np
from scipy.special import gamma
import pandas as pd
import integration as itg
import mittag_leffler as ml

# This script uses the old definition of constants
# Now for a more general description we use
## alpha instead of H
## A -> M
## eta -> zeta
## C -> gamma
## theta_H -> A_H (depends on the physical constants and the model)
## theta_12 -> A (depends on the physical constants and the model)
class fle_twobath():
    def __init__(self, H):
        self.H = H
    
    def params(self, T, h = 0.01, v0 = 1, M = 1, eta_H = 1, eta_12 = 1, theta_H = 1, theta_12 = 1, fluctuation_dissipation = True):
        """
        Function to set the parameters of the system to solve
        


        Args:
            T (int): Final time. Upper bound of evaluation.
            h (float, optional): Size of time step. Defaults to 0.01.
            v0 (int, optional): Initial velocity of the system. Defaults to 1.ç
            A (float, optional): Amplitud of second order derivative. Defaults to 1.
            eta_H (float, optional): Amplitud of fractional derivative. Defaults to 1.
            eta_12 (float, optional): Amplitud of first order derivative. Defaults to 1.
            theta_H (float, optional): Amplitud of first order derivative. Defaults to 1.
            theta_12 (float, optional): Amplitud of first order derivative. Defaults to 1.
            fluctuation_dissipation (bool, optional): If the system preserves fluctuation-dissipation, the model constructs the proper
            coefficientes. Defaults to True
        """
        self.kBT = 1
        self.T = T
        self.h = h
        self.n = int(self.T/self.h)
        self.x_n = np.zeros(self.n)

        self.v0 = v0
        self.M = M
        self.eta_H = eta_H
        self.eta_12 = eta_12
        self.theta_H = theta_H
        self.theta_12 = theta_12
        if fluctuation_dissipation:
            kBT = self.kBT
            H = self.H
            
            self.theta_12 = np.sqrt(2*kBT*eta_12)*theta_12
            self.theta_H = np.sqrt(kBT*eta_H/(H*(2*H-1)*gamma(2*H-1)))*theta_H
        
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
        theta_H = self.theta_H
        theta_12 = self.theta_12
        np.random.seed()
          
        f = FBM(n = n, hurst = H, length = T, method = method)
        self.B_H = f.fbm()
        self.dB_H = n/T*np.diff(self.B_H)
               
        f12 = FBM(n = n, hurst = 0.5, length = T, method = method)
        self.B_12 = f12.fbm()
        self.dB_12 = n/T*np.diff(self.B_12)
        
        self.t_BH = f.times()
        self.t = self.t_BH[:self.n]
        self.B = theta_12*self.B_12 + theta_H*self.B_H
        self.dB = theta_12*self.dB_12 + theta_H*self.dB_H
        
    def get_analytical(self, relaxation_type = "sub"):
        H = self.H
        eta_H = self.eta_H
        M = self.M
        eta_12 = self.eta_12 
        t = self.t
        noise = self.dB
        v0 =self.v0
                
        t__ = np.array(t)
        noise__ = np.array(noise)
        
        def relaxation_function_sup(t):            
            order = 2 - 2*H
            z = -(eta_H/M)*t**(2-order)
            G = pd.DataFrame()
            inf = 40
            for n in range(inf):
                t_ = (t**(1+n)) * ((-eta_12/M)**n)
                G[f"n{n}"] = ml.prabhakar_mittag_leffler(z, 2-order, 2 + n, n+1) * t_
            
            return np.array(G.sum(axis = 1))
        
        def relaxation_function_sub(t):            
            order = 2 - 2*H
            z = -(eta_12/M)*t
            G = pd.DataFrame()
            inf = 40
            for n in range(inf):
                t_ = ((-eta_H/M)**n)*t**((2-order)*n+1)
                G[f"n{n}"] = ml.prabhakar_mittag_leffler(z, 1 , (2-order)* n + 2, n+1) * t_ 
            
            return np.array(G.sum(axis = 1))
        
        if relaxation_type == "sub":
            relaxation_function = relaxation_function_sub
        if relaxation_type == "sup":
            relaxation_function = relaxation_function_sup
        
        conv = itg.convolution(relaxation_function, noise__, t__)
        
        nonlinear = v0 * relaxation_function(t__)
        
        self.analytical = nonlinear + (1/M) * conv
    
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
        eta_H = self.eta_H
        eta_12 = self.eta_12
         
        a_jj = 0
        if k-1 >= 1:    
            for j in range(1, k):
                a_jj += self.a_j(j)*(x_n[k-j]+x_n[k-j-1])
        h2Hm = h**(2*H)/((M+h*eta_12)*2*(2*H-1)*gamma(2*H-1))
        return (M/(M+h*eta_12))*x_n[k-1] + (h*M/(M+h*eta_12))*v0 - eta_H*h2Hm*a_jj + (h/(M+h*eta_12))*B[k]
    #TODO: se a_j for negative powers

    def solve(self):
        n = self.n
        x_n = self.x_n
        
        for k in range(1,n):
            x_n[k] = self.x_k(k)
        
        self.numerical = x_n
