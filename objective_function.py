"""
Loads output of 2D_beam_iteration.py and makes a plot of the objective 
function.
   
Supplementary Material for "Optimal Targeting of a Tumor through Proton Beam 
Therapy", Journal of Young Investigators.
Authors: Kiran Pant and Colin Campbell
Python 3.6.9
"""

import numpy as np
import matplotlib.pyplot as plt        
import pickle

E_total = np.loadtxt(r'../Python Output/2D_Energy_total_array.txt')
E_tumor = np.loadtxt(r'../Python Output/2D_Energy_tumor_array.txt')
f = open(r'../Python Output/array_stats.data','rb')
array_stats = pickle.load(f)
f.close()

f = 2*E_tumor - E_total

mlist = ['o','s','^','*','d']
clist = ['k','r','g','b','m']

plt.figure()
for i in range(np.shape(f)[0]):
    plt.plot(array_stats['Evals'],f[i,:],marker=mlist[i],c=clist[i],label='s = {0:.2f}'.format(array_stats['svals'][i]))
plt.legend(loc='best',ncol=5,fontsize=8)
plt.xlabel('Initial Beam Energy (MeV)')
plt.ylabel('Objective Function (MeV)')
plt.xticks(np.arange(83,107,2))
plt.ylim(-50,25)
plt.show()

