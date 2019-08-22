"""
Loads data from 3D_beam.py and plots a surface containing 80% of energy
deposited as depth varies.

Supplementary Material for "Optimal Targeting of a Tumor through Proton Beam 
Therapy", Journal of Young Investigators.
Authors: Kiran Pant and Colin Campbell
Python 3.6.9
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib import cm
from mpl_toolkits.mplot3d.axes3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
import pickle
import matplotlib.cm as cm
from scipy.signal import savgol_filter as sgf

f = open(r'../Python Output/3D_plot_80_0-1.data','rb')
data = pickle.load(f)
f.close()

rmax = (2**0.5)*max(data['xvals'])
zmax = max(data['zvals'])

N = 50                                     # how many values to have in our z,r array
array_zvals = np.linspace(0,zmax,N)
array_rvals = np.linspace(0,rmax,N)
A = np.zeros((N,N))
tumor_energy = 0

target_center = 6
target_radius = 2

# Populate array that stores energy indexed by depth z and radial distance r
for i in range(len(data['xvals'])):
    # For array
    r = (data['xvals'][i]**2 + data['yvals'][i]**2)**0.5
    z_index = np.searchsorted(array_zvals,data['zvals'][i])
    r_index = np.searchsorted(array_rvals,r)
    A[z_index,r_index] += data['energy_vals'][i]
    
# Now for each z, what radius corresponds to encompassing 50% of total?
radius_list = []
energy_list = []
for z_index in range(N):
    row_sum = sum(A[z_index,:])
    energy_list.append(row_sum)
    tot = 0
    for r_index in range(N):
        tot += A[z_index,r_index]
        if tot >= 0.5*row_sum:
            radius_list.append(array_rvals[r_index])
            break
        
#radius_list can be noisy; smooth it
radius_list_smooth = sgf(radius_list,27,3)            
        
# Now we need to plot this surface (use cylindrical coordinates)
theta = np.linspace(0,2*np.pi,N)
Z,THETA = np.meshgrid(array_zvals,theta)
R = np.outer(np.ones(N),radius_list_smooth)

energy_min = min(energy_list)
energy_max = max(energy_list)
energy_list_colors = [(energy_min+x)/energy_max for x in energy_list]
C = np.outer(np.ones(N),energy_list_colors)
X = R * np.cos(THETA)
Y = R * np.sin(THETA)

# Separately, we plot a sphere representing the target (centered at 6 cm, radius 2 cm)
# spherical coordinates
theta, phi = np.linspace(0,np.pi,100),np.linspace(0,2*np.pi,100)
THETA, PHI = np.meshgrid(theta,phi)
R = np.outer(np.ones(100),target_radius)

X_target = R*np.sin(THETA)*np.cos(PHI)
Y_target = R*np.sin(THETA)*np.sin(PHI)
Z_target = target_center+R*np.cos(THETA)

norm = matplotlib.colors.Normalize(vmin=0, vmax=1)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
surf = ax.plot_surface(X,Y,Z,facecolors=plt.cm.plasma(norm(C)),alpha=1,linewidth=0,rcount=100,ccount=100)
ax.plot_wireframe(X_target,Y_target,Z_target,color='k',alpha=0.15,linewidth=1,rcount=10,ccount=10)

#ax.set_aspect('equal')
# label axes to match orientation in 2D projection (beam points along X on entry)
ax.set_xlabel('Y (cm)')
ax.set_ylabel('Z (cm)')
ax.set_zlabel('X (cm)')
#ax.grid(False)

#colorbar
m = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
m.set_array([])
plt.colorbar(m,shrink=0.5,pad=.15)

plt.figtext(0.75, 0.495, "Relative Energy Deposition",rotation='vertical',horizontalalignment='center',verticalalignment='center')

#plt.savefig('3D Optimal Beam.png',dpi=800)

plt.show()

