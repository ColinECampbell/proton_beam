"""
Implementation of the 1D Bragg Kleeman equation.

Supplementary Material for "Optimal Targeting of a Tumor through Proton Beam 
Therapy", Journal of Young Investigators.
Authors: Kiran Pant and Colin Campbell
Python 3.6.9
"""

import matplotlib.pyplot as plt
import numpy as np

# constants for water (Newhauser table 2)
rho = 1             #g cm^(-3)
alpha = 2.633E-3
p = 1.735


dx=1E-3                     # differential step (depth into material)
Estart=[150,200,250]        # starting energies (MeV)

# Below iteratively calculates energy loss as fn of depth until all energy is 
# lost, then notes that 0 energy remains until x = 50
plot_xvals = []
plot_yvals = []
for E in Estart:
    x=0 
    x_vals = []
    dEdx_vals = []
    while E>0: 
        x=x+dx  
        dEdx=E**(1-p)/(rho*alpha*p) 
        dE=dEdx*dx 
        E=E-dE 
        if E < 0: break
        x_vals.append(x)
        dEdx_vals.append(dEdx)
    while x < 50:
        x = x + dx
        x_vals.append(x)
        dEdx_vals.append(0)
    
    plot_xvals.append(x_vals)
    plot_yvals.append(dEdx_vals)

# make the plots  
plt.figure()
plt.plot(plot_xvals[0],plot_yvals[0],'b--',lw=2,label='150 MeV')
plt.plot(plot_xvals[1],plot_yvals[1],'g-.',lw=2,label='200 MeV')
plt.plot(plot_xvals[2],plot_yvals[2],'r-',lw=2,label='250 MeV')
plt.legend(loc='best')
plt.xlabel('depth (cm)',fontsize=15)
plt.ylabel(r'$\mathrm{d}E/\mathrm{d}x\;(\mathrm{MeV}\;\mathrm{cm}^2\;\mathrm{g}^{-1})$',fontsize=15)
plt.show()


# or on a log scale
plt.figure()
plt.semilogy(plot_xvals[0],plot_yvals[0],'b--',lw=2,label='150 MeV')
plt.semilogy(plot_xvals[1],plot_yvals[1],'g-.',lw=2,label='200 MeV')
plt.semilogy(plot_xvals[2],plot_yvals[2],'r-',lw=2,label='250 MeV')
plt.legend(loc='best')
plt.xlabel('depth (cm)',fontsize=15)
plt.ylabel(r'$\mathrm{d}E/\mathrm{d}x\;(\mathrm{MeV}\;\mathrm{cm}^2\;\mathrm{g}^{-1})$',fontsize=15)
plt.show()