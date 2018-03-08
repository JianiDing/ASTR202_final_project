# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 11:46:47 2018

@author: Brent
"""

# Start with standard imports 
import matplotlib.pylab as plt
import numpy as np
from scipy import interpolate

### Interpolation of free-free gaunt factor ###

# read in the data table for free free gaunt factor
# this is basically our 'z' data
g_table = 'gauntff_noerr.dat'
g_data = np.loadtxt(g_table, skiprows=42) # skipping the header

# follow the header to get log(u) and log(gam2) values for
# interpolation   
g_logu = np.arange(-16, 13.2, 0.2)
g_loggam2 = np.arange(-6, 10.2, 0.2)

# create the interpolation
g_interp = interpolate.interp2d(g_loggam2, g_logu, g_data)

### Evaluation of formula ###

# fundamental constants
Ry = 2.17987e-11 # erg
h = 6.62606885e-27 # erg * s
k = 1.3807e-16 # erg * k^-1
c = 2.998e10 # cm / s

# Test cases. The numerical values of these are in the MNRA paper. 
# The values match relatively closely
#print(np.log10(Ry/(k * 1e10)))
#print(np.log10(Ry* 30/(k * 3)))
#print(np.log10(2.17e-19/(k * 1e10)))
#print(np.log10(1.603e-4 * 30/(k*3)))

# get gam2 and u in terms of T and nu
gam2 = lambda T, nu: (h * (nu))/(k * T)
u = lambda T: Ry / (k * T)

# main function to be evaluated 

def j (T, lamb):
    # convert a wavelength (lamb) in micrometers to a frequency
    nu = (c / (lamb * 10**-4))
    
    # interpolate the gaunt factor
    gff = g_interp(np.log10(gam2(T, nu)),np.log10(u(T)))
    print(-h * nu/(k * T))
    print('nu is: ', nu)
    print('gff is: ', gff)
    print('logu is: ', np.log10(u(T)))
    print('loggam2 is: ', np.log10(gam2(T, nu)))
    
    # evaluate the function 
    return 5.44 * 10**(-39) * (T**(-1/2)) * gff * np.exp(-h * nu/(k * T)) * nu

x_values = np.logspace(-1, 1, 1000)
y_values = j(10000, x_values)

plt.loglog(x_values, y_values)
plt.show()