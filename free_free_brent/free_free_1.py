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
#print(np.log10(Ry/(k * 1e10))) # gam2 low
#print(np.log10(Ry* 30**2/(k * 3))) # gam2 high
#print(np.log10(2.17e-19/(k * 1e10))) # u low
#print(np.log10(1.603e-4/(k*3)))# u high

# get gam2 and u in terms of T and nu
gam2 = lambda T: Ry/(k * T)
u = lambda T, nu: (h * nu) / (k * T)


# main function to be evaluated 

def j (T, lamb):
    # convert a wavelength (lamb) in micrometers to a frequency
    nu = (c / (lamb * 10**-4))
    
    # interpolate the gaunt factor
    gff = g_interp(np.log10(gam2(T)),np.log10(u(T, nu)))
    
    # transpose the matrix. This is just for consistensy sake with the 
    # horizontal inputs 
    gff = np.resize(gff, (1000,))
    
    # evaluate the function 
    output = 5.44 * 10**(-39) * (4*np.pi) * (T**(-1/2)) * gff * np.exp(-h * nu/(k * T)) * nu
    return output

x_values = np.logspace(-1, 1, 1000)
y_values = j(10000, x_values)
#print(x_values)

plt.loglog(x_values, y_values)
plt.show()