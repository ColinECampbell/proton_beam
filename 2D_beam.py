'''
Creates energy dispersion figures in the case of 2D beams, either directly or
by loading an output file from 2D_beam_iteration.py.

Supplementary Material for "Optimal Targeting of a Tumor through Proton Beam 
Therapy", Journal of Young Investigators.
Authors: Kiran Pant and Colin Campbell
Python 3.6.9
'''

import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.stats import norm
from scipy.integrate import simps
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm
from matplotlib.patches import Circle, PathPatch
import mpl_toolkits.mplot3d.art3d as art3d
import time

# FUNCTION DEFINITION =========================================================
def determine_dE(Estart):
    """
    Computes the B-K curve. Returns x coordinates, dEdx values as two lists.
    Estart is the initial energy.
    """
    # constants for water (Newhauser table 2)
    rho = 1     #g cm^(-3)
    alpha = 2.633E-3
    p = 1.735
    
    # Below iteratively calculates energy loss as fn of depth until all energy is lost
    x=0 
    E=Estart 
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

    # Append a few "0" entries at the end (ease of plotting)    
    for i in range(10):
        x_vals.append(x+dx)
        dEdx_vals.append(0)
    
    return x_vals, dEdx_vals        

def calculate_row(Estart):
    
    # Set up container arrays in polar coordinates
    rvals, thetavals = np.linspace(dx,L,int(N-1)), np.linspace(-np.pi/2,np.pi/2,int(N-1))

    # Output array of energy deposition
    a = np.zeros([len(rvals),len(thetavals)],float)
    
    # parameters for rms spread (Abril Table 1 & eq 5)
    a1 = -0.058
    a2 = -1.868
    b1 = 9.39E-3
    b2 = 1.56E-3
    C1 = b1*(Estart**a1)
    C2 = b2*(Estart**a2)
    
    # populate entries of array via a double loop
    for i,r in enumerate(rvals):
        # determine dEdx for this depth
        dEdx_index = np.searchsorted(linear_x,r)
        if dEdx_index == len(linear_dEdx):
            dEdx_val = 0
        else:
            dEdx_val = linear_dEdx[dEdx_index]

        # determine spread of Gaussian at this depth    (Abril eqn 4)     
        sigma = 0.1+(C1*(r*1E4)+C2*(r*1E4)**2)/(r*1E4)  # in radians
        
        
        for j,theta in enumerate(thetavals): 
            # store energy deposited at this theta, r
            # We want dE/dA * r, but the expression simplifies:
            # dE = (dE/dr)*dr*(norm()*dtheta)
            # dA = r*dr*dtheta
            # Thus we end up with just dE/dr * norm()
            a[i,j] += norm.pdf(theta,0,sigma) * dEdx_val
    
    num_integral = simps(simps(a, thetavals), np.array(rvals))
    print(num_integral)
            
    return a
    
# BEGIN MAIN ALGORITHM ========================================================
TIME0 = time.time()    

L = 10                                          # How far into the material will we consider?
dx = 0.05                                       # Differential step of energy deposition as fn of depth
N = L/dx                                        # Number of points to consider as fn of depth

# Use below to run trial and save data ...
#E0 = 80                                         # Initial energy 
#linear_x,linear_dEdx = determine_dE(E0)         # Energy deposition as fn of depth
#master_array = calculate_row(E0)
#np.savetxt('proton_beam_sample_figure.txt',master_array)

# ... or below to load previous results for plotting
master_array = np.loadtxt(r'../Python Output/2D_proton_beam_array_99_0-1.txt')

# MAKE FIGURES ================================================================
rvals, thetavals = np.linspace(dx,L,int(N-1)), np.linspace(-np.pi/2,np.pi/2,int(N-1))
Rvals, THETAvals = np.meshgrid(rvals,thetavals)

# make a plot of the 2D dispersion ---------------------
X, Y = Rvals*np.cos(THETAvals), Rvals*np.sin(THETAvals)

fig = plt.figure(figsize=(9,4))
ax0 = fig.add_subplot(1, 2, 1, projection='3d')
surf = ax0.plot_surface(X, Y, master_array.T, rstride=1, cstride=1, cmap="plasma",\
                       linewidth=0, alpha = 0.8, antialiased=False)

ax0.set_yticks([-10,-5,0,5,10])
ax0.set_zticks([])
ax0.set_zticklabels([],{'fontsize':'x-small'})
ax0.set_xlabel('\nx (cm)', fontsize=10)
ax0.set_ylabel('\ny (cm)',fontsize=10)

ax0.azim = 270
ax0.elev = 90
p0 = Circle((6, 0.0),2, fc='k',alpha=0.5)
ax0.add_patch(p0)
art3d.pathpatch_2d_to_3d(p0, z=0, zdir="z")

ax0.annotate("target",
            xy=(0.63, 0.52), xycoords='axes fraction',
            xytext=(0.75, 0.85), textcoords='axes fraction',
            arrowprops=dict(arrowstyle="fancy",
                            connectionstyle="arc3,rad=-0.3", ec = 'g',fc='k'),
            fontsize=10 
            )
            
ax0.annotate("beam entry point",
            xy=(0.28, 0.52), xycoords='axes fraction',
            xytext=(0.15, 0.1), textcoords='axes fraction',
            arrowprops=dict(arrowstyle="fancy",
                            connectionstyle="arc3,rad=-0.3", ec = 'g',fc='k'),
            fontsize=10 
            )

#----

ax = fig.add_subplot(1, 2, 2, projection='3d')
X, Y = Rvals*np.cos(THETAvals), Rvals*np.sin(THETAvals)
surf = ax.plot_surface(X, Y, master_array.T, rstride=1, cstride=1, cmap="plasma",\
                       linewidth=0, alpha = 0.8, antialiased=False)

ax.set_xlabel('x (cm)', fontsize=10)
ax.set_ylabel('y (cm)',fontsize=10)
ax.set_zlabel('\n\n\n\nMass Stopping Power\n(MeV '+'$\mathrm{g}^{-1}\ \mathrm{cm}^{-2}$'+')\n\n\n',fontsize=10)
ax.azim = 300

# This bit to add a circular patch representing the tumor
p = Circle((6, 0.0), 2, fc='k',alpha=0.5)
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=0, zdir="z")

# Add a colorbar
plt.colorbar(surf,shrink=0.5,pad=.15)

plt.subplots_adjust(top=0.954,bottom=0.046,left=0.023,right=0.9,hspace=0.29,wspace=0.01)
#plt.savefig('Appendix C.png',dpi=800)

plt.show()

TIME1 = time.time()
print('Algorithm completed in {0:.2f} minutes.'.format((TIME1-TIME0)/60))