'''
Line emission spectrum of H at T=1e4K

The lines have Gaussian form.

-Platon Karpov
'''
from __future__ import print_function
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

def main():
    h = 6.626e-27
    c_ph = 2.99792e10
    rate = 1.47234e+53  #Ly photon rate
    
    #Find emission coefficient of H-beta
    lam = Rydberg(2,4)
    nu = c_ph/lam
    d = 7e28 #cm
    jb = h*nu*rate/(4*np.pi*(d)**2)*1e-29 #Flux of H-beta in Janskies
    
    #Table 4.2 of O&F
    j = [[2.87, 1, 0.466,0.256,0.158,0.105,0.0730,0.0529,0.0154,0.0064],
         [0.352,0.354,0.354,0.350,0.350,0.344,0.344]]    
    n2 = [[3,4,5,6,7,8,9,10,15,20],[4,5,6,8,10,15,20]]
    n1 = [2,3]
    balmer = []
    plt.figure(figsize=(7,5))
    f = open('line_flux.txt', 'w+')
    f.write('lambda[A]\tFlux[microJy]\tError\n')
    print('lambda[A]\tFlux[microJy]\tError')
    
    for k in range(len(n2)):
        for i in range(len(j[k])):                    
            lam = Rydberg(n1[k], n2[k][i])
            
            #parameters for the Gaussian
            b = lam/1e-8
            c = 1e-4*b
            x = np.linspace(b-5*c, b+5*c, 100)
            if n1[k]==2:
                if n2[k][i] in n2[1]:
                    balmer.append(j[k][i]*jb)
                a = j[k][i]*jb
            elif n1[k]==3:
                a = balmer[i]*j[k][i]
            
            #integrate each line
            inter, err = integrate.quad(lambda x: gauss(x,lam,a,b,c), b-5*c, b+5*c)
            f.write('%.3e\t%.3e\t%.3e\n'%(b,inter,err))
            print('%.3e\t%.3e\t%.3e'%(b,inter,err))
            
            plt.plot(x, gauss(x,lam,a,b,c), color = 'k')


    plt.axhline(y=0, color='k')
    plt.xlim(1e3, 2e4)
    plt.xlabel('Wavelength [$\AA$]')
    plt.ylabel('Flux [$erg \; cm^{-2} \; s^{-1} \; Hz^{-1}$]')
    plt.title('Hydrogen line spectrum')

def Rydberg(n1,n2):
    R = 1.097e5
    return (R*(1/n1**2-1/n2**2))**(-1)

def gauss(x,lam,a,b,c):
    return a*np.exp(-(x-b)**2/(2*c**2))

if __name__ == '__main__':
    main()