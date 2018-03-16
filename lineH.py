'''
Line emission spectrum of H at T=1e4K

-Platon Karpov
'''
from __future__ import print_function
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def main():
    #Ne = Np, low density limit  
    N = 1e0
    h = 6.626e-27
    alpha = 3.03e-14    #for H-beta from Table 4.2
    c = 2.99792e10
    
    #Find emission coefficient of H-beta
    lam = Rydberg(2,4)
    nu = c/lam
    jb = alpha*h*N**2*nu/(4*np.pi)  #Eq.4.14 of O&F
    
    #Table 4.2 of O&F
    j = [[2.87, 1, 0.466,0.256,0.158,0.105,0.0730,0.0529,0.0154,0.0064],[0.352,0.354,0.354,0.350,0.350,0.344,0.344]]    
    n2 = [[3,4,5,6,7,8,9,10,15,20],[4,5,6,8,10,15,20]]
    n1 = [2,3]
    balmer = []
    plt.figure(figsize=(7,5))
    for k in range(len(n2)):
        for i in range(len(j[k])):                    
            lam = Rydberg(n1[k], n2[k][i])
            if n1[k]==2:
                if n2[k][i] in n2[1]:
                    balmer.append(j[k][i]*jb)
                plt.vlines(x=lam/1e-8, ymin=0, ymax = j[k][i]*jb)
            elif n1[k]==3:
                plt.vlines(x=lam/1e-8, ymin=0, ymax = balmer[i]*j[k][i])                

    plt.axhline(y=0, color='k')
    plt.xlim(1e3, 2e4)
    plt.xlabel('Wavelength [$\AA$]')
    plt.ylabel('j [$erg \; cm^{-3} \; s^{-1}$]')
    plt.title('Hydrogen line spectrum')

def Rydberg(n1,n2):
    R = 1.097e5
    return (R*(1/n1**2-1/n2**2))**(-1)

if __name__ == '__main__':
    main()