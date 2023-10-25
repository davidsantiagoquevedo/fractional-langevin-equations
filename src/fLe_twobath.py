from fbm import FBM
import numpy as np
from scipy.special import gamma

class fle_twobath():
    def __init__(self, H):
        self.H = H
    
    def params(self, T, h = 0.01, v0 = 1, A = 1, eta = 1, C = 1, theta_H = 1, theta_12 = 1,):
        """
        Function to set the parameters of the system to solve
        


        Args:
            T (int): Final time. Upper bound of evaluation.
            h (float, optional): Size of time step. Defaults to 0.01.
            v0 (int, optional): Initial velocity of the system. Defaults to 1.ç
            A (float, optional): Amplitud of second order derivative. Defaults to 1.
            eta (float, optional): Amplitud of fractional derivative. Defaults to 1.
            C (float, optional): Amplitud of first order derivative. Defaults to 1.
            theta_H (float, optional): Amplitud of first order derivative. Defaults to 1.
            theta_12 (float, optional): Amplitud of first order derivative. Defaults to 1.
        """
        self.T = T
        self.h = h
        self.n = int(self.T/self.h)
        self.x_n = np.zeros(self.n)

        self.v0 = v0
        self.A = A
        self.eta = eta
        self.C = C
        self.theta_H = theta_H
        self.theta_12 = theta_12
        
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
        A = self.A
        eta = self.eta
        C = self.C
         
        a_jj = 0
        if k-1 >= 1:    
            for j in range(1, k):
                a_jj += self.a_j(j)*(x_n[k-j]+x_n[k-j-1])
        h2Hm = h**(2*H)/((A+h*C)*2*(2*H-1)*gamma(2*H-1))
        return (A/(A+h*C))*x_n[k-1] + (h*A/(A+h*C))*v0 - eta*h2Hm*a_jj + (h/(A+h*C))*B[k]
    #TODO: se a_j for negative powers

    def solve(self):
        n = self.n
        x_n = self.x_n
        
        for k in range(1,n):
            x_n[k] = self.x_k(k)
        
        self.numerical = x_n
