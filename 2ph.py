'''
Two-photon emission from H-nebula.

Call two_ph(T,N) to use as a blackbox

Inputs:
T: temperature [K]
N: is number density [#/cm^3] (Ne~Np for low rho, hence only one N is required as input) 

Output:
lam: wavelength [micron]
ynu: emission coefficient multiplied by frequency [erg cm^3 s^-1]

-Platon Karpov
'''
from __future__ import print_function
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def main():
    #Ne = Np, low density limit
    T = 1e4
    N = 1e0       

    nu = np.logspace(14,16,1000)
    lam, ynu = two_ph(T,N)

    plt.figure(figsize=(7,5))
    plt.loglog(lam,ynu, basex=10, basey=10)
    plt.xlim(0.1, 10)
    plt.ylim(1e-26, 2e-23)
    plt.xlabel('Wavelength [$\mu m$]')
    plt.ylabel(r'$\nu\gamma_{\nu} \; [erg \; cm^3 \; s^{-1}]$')
    plt.title("Two-photon emission (Brant's 4.1)")
    
    plt.figure(figsize=(5,8))
    plt.semilogy(nu,ynu/nu, basey=10)
    plt.xlim(3e14,1e15)
    plt.ylim(1e-40,1e-38)
    plt.xlabel('Frequency')
    plt.ylabel(r'$\gamma_{\nu} \; [erg \; cm^3 \; s^{-1} \; Hz^{-1}]$')
    plt.title("Two-photon emission (Zack's 4.1)")    
    

def two_ph(T,N):
    #compute two-photon emission
    #import this function to use as a blackbox
    nu = np.logspace(14,16,1000)
    A = 8.23
    ynu = alpha(T)*g2(nu)/(1+N*q(T)/A)*nu
    
    c = 2.99792e8
    lam = [c/i/1e-6 for i in nu]    #convert to wavelength in microns
    return lam,ynu
    
def alpha(T):
    #absorption coefficients
    Temps = [5e3, 1e4, 1.5e4, 2e4]
    a = [1.38, 0.838, 0.625, 0.506]
    return np.interp(T, Temps, a)*1e-13

def g(nu):
    #from the Table 4.12 of Osterbrock & Ferland
    #rmk: it does not have the range necessary
    nu_vals = [0.0,1.23,2.47,3.70,4.93,6.17,7.40,8.64,9.87,11.10,12.34]
    nu_vals = [i*1e14 for i in nu_vals]
    g = [0.0,0.303,0.978,1.836,2.78,3.78,4.8,5.8,6.78,7.74,8.62]
    g = [i*1e-27 for i in g]

    return np.interp(nu, nu_vals, g)

def g2(nu):
    #fitting formula for density dependence from Schirmer 2016
    h = 6.626e-27
    nu12 = 1.23e15*2
    A2q=8.23
    beta=1.53
    gamma=0.8
    delta=0.88
    C=202
    y = nu/nu12*(1-nu/nu12)
    A = C*(y*(1-(4*y)**gamma)+delta*y**beta*(4*y)**gamma)
    return h*nu/nu12*A/A2q

def q(T):
    #collisional transition rates
    Temps = [1e4, 2e4]
    q = [5.37e-4, 4.75e-4]
    return np.interp(T, Temps, q)
    
if __name__ == '__main__':
    main()

