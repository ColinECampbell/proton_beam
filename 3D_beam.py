'''
Generates a 3D beam and saves the data to file. Plotted in 3D_beam_plot.py.

Supplementary Material for "Optimal Targeting of a Tumor through Proton Beam 
Therapy", Journal of Young Investigators.
Authors: Kiran Pant and Colin Campbell
Python 3.6.9
'''

import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.stats import norm, uniform
from scipy.integrate import simps
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.patches import Circle, PathPatch
import mpl_toolkits.mplot3d.art3d as art3d
import matplotlib.cm as cm
import pickle
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
    dE_vals = []
    while E>0: 
        x=x+dx  
        dEdx=E**(1-p)/(rho*alpha*p) 
        dE=dEdx*dx 
        E=E-dE 
        if E < 0: break
        x_vals.append(x)
        dEdx_vals.append(dEdx)
        dE_vals.append(dE)
        
    # Append a few "0" entries at the end (ease of plotting)    
    for i in range(10):
        x_vals.append(x+dx)
        dEdx_vals.append(0)
        dE_vals.append(0)
    
    return x_vals, dEdx_vals, dE_vals      

def calculate_row(Estart):
    
    # Set up container arrays in polar coordinates
    rvals, thetavals = np.linspace(0.01,L,N), np.linspace(-.99*np.pi/2,.99*np.pi/2,N)
    phivals = np.linspace(0,2*np.pi,N)
    
    # Output array of energy deposition
    a = np.zeros([N,N,N],float)
    
    x,y,z,dEdV,total_energy = [],[],[],[],[]


    # parameters for rms spread (Abril Table 1 & eq 5)
    a1 = -0.058
    a2 = -1.868
    b1 = 9.39E-3
    b2 = 1.56E-3
    C1 = b1*(Estart**a1)
    C2 = b2*(Estart**a2)
        
    # populate entries of array via a double loop
    for i,r in enumerate(rvals):
        if i%10 == 0: print("starting i={0}".format(i))
        # determine dEdx for this depth
        dEdx_index = np.searchsorted(linear_x,r)
        if dEdx_index == len(linear_dEdx):
            dEdx_val = 0
            dE_val = 0
        else:
            dEdx_val = linear_dEdx[dEdx_index]
            dE_val = linear_dE[dEdx_index]

        # determine spread of Gaussian at this depth    (Abril eqn 4)     
        sigma = 0.1+(C1*(r*1E4)+C2*(r*1E4)**2)/(r*1E4)  # in radians
       
        for j,theta in enumerate(thetavals):
            for k,phi in enumerate(phivals): 
                # store energy deposited at this theta, r, phi
                # We want dE/dV * (r**2 * sin(theta)), but the expression simplifies:
                # dE = (dE/dr)*dr*(norm()*dtheta)*(uniform()*dphi)
                # dV = r**2 *sin(theta) * dr * dtheta * dphi
                # Thus we end up with just (dE/dr)*(norm())*(uniform())
                a[i,j,k] += norm.pdf(theta,0,sigma) * dEdx_val * uniform.pdf(phi,0,2*np.pi)
                # and store Cartesian data for output
                x.append(r * np.sin(theta) * np.cos(phi))
                y.append(r * np.sin(theta) * np.sin(phi))
                z.append(r * np.cos(theta))
                total_energy.append(a[i,j,k])
    
    t = simps(simps(simps(a, phivals),thetavals),rvals)
    
    print(t) 
       
    return x,y,z,total_energy
    
# BEGIN MAIN ALGORITHM ========================================================
TIME0 = time.time()    

L = 10                                                      # How far into the material will we consider?
dx = 0.05
N = int(L/dx)
dphi = 2*np.pi/N
dtheta = np.pi/N
E0 = 99                                                     # Initial energy
s = 0.1                                                     # Beam spread
linear_x,linear_dEdx,linear_dE = determine_dE(E0)           # Energy deposition as fn of depth


x,y,z,total_energy = calculate_row(E0)

energy_min,energy_max = min(total_energy),max(total_energy)
energy_range = abs(energy_max-energy_min)
energy_normed = [(x+energy_min)/energy_range for x in total_energy]
#
# Trim data
x_plt,y_plt,z_plt,s_plt,c_plt = [],[],[],[],[]
for i in range(len(x)):
    if energy_normed[i] > 0.01:
        x_plt.append(x[i])
        y_plt.append(y[i])
        z_plt.append(z[i])
        s_plt.append(10)
        c_plt.append(total_energy[i])

# Save data?
data = {'E0':E0,'xvals':x_plt,'yvals':y_plt,'zvals':z_plt,'energy_vals':c_plt}
E_str = str(E0).replace('.','-')
s_str = str(s).replace('.','-')
f = open(r'../Python Output/3D_plot_'+E_str+'_'+s_str+'.data','wb')
pickle.dump(data,f)
f.close()


# MAKE FIGURES ================================================================

# (Don't normally generate this figure -- just use 3D_beam_plot.py)
#fig = plt.figure(figsize=(8,4))
#ax = fig.add_subplot(1, 1, 1, projection='3d')
#ax.scatter(x_plt,y_plt,z_plt,alpha=0.2,s=s_plt,c='.5')
#plt.tight_layout()
#
#plt.show()

#print("Numerical energy = {0:.2f}, or {1:.2f} % error.".format\
#      (float(sum(total_energy)),float(100*abs(sum(total_energy)-E0)/E0)))

TIME1 = time.time()
print('Algorithm completed in {0:.2f} minutes.'.format((TIME1-TIME0)/60))

