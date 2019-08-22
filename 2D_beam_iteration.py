"""
Calculates 2D beam properties for a range of energies & widths. Stores each
beam to an array and the total energy, tumor energy arrays for objective 
function calculations.

Supplementary Material for "Optimal Targeting of a Tumor through Proton Beam 
Therapy", Journal of Young Investigators.
Authors: Kiran Pant and Colin Campbell
Python 3.6.9
"""
import numpy as np
from scipy.stats import norm
from scipy.integrate import simps
import time
import pickle

def master(E0,s):    
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
        a = np.zeros([int(N-1),int(N-1)],float)
        
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
            sigma = s+(C1*(r*1E4)+C2*(r*1E4)**2)/(r*1E4)  # in radians            
            
            for j,theta in enumerate(thetavals): 
                # store energy deposited at this theta, r
                # We want dE/dA * r, but the expression simplifies:
                # dE = (dE/dr)*dr*(norm()*dtheta)
                # dA = r*dr*dtheta
                # Thus we end up with just dE/dr * norm()
                a[i,j] += norm.pdf(theta,0,sigma) * dEdx_val
        
        num_integral = simps(simps(a, thetavals), np.array(rvals))
        pct_error = 100*abs(E0-num_integral)/(E0)
        if pct_error >= 5:
            # this rarely occurs; just re-run the trial if so
            return 'bad run','bad run'
        else:
            print("Initial energy = {0:.3g} vs. Numerical Integral = {1:.3g}".format(E0,num_integral))
            print("Percent Error = {0:.2f}".format(pct_error))
            return a, num_integral
    
    # BEGIN MAIN ALGORITHM ========================================================
    L = 10                                          # How far into the material will we consider?
    dx = 0.05                                       # Differential step of energy deposition as fn of depth
    N = L/dx                                        # Number of points to consider as fn of depth
         
    while True:
        linear_x,linear_dEdx = determine_dE(E0)         # Energy deposition as fn of depth
        master_array, num_integral = calculate_row(E0)
        
        if type(master_array) == str:
            dx /= 2
            N = L/dx
            print('High error; restarting.')
        else:
            break
        
    # Save data?
    E_str = str(E0).replace('.','-')
    s_str = str(s).replace('.','-')
    np.savetxt(r'../Python Output/2D_proton_beam_array_'+E_str+'_'+s_str+'.txt',master_array)
        
    
    # EVAL TUMOR ==============================================================
    rvals, thetavals = np.linspace(dx,L,int(N-1)), np.linspace(-np.pi/2,np.pi/2,int(N-1))
    Rvals, THETAvals = np.meshgrid(rvals,thetavals)
    
    rho = 2
    d = 6
    tumor = np.zeros([int(N-1),int(N-1)],float)
    for i,r in enumerate(rvals):
        for j,theta in enumerate(thetavals):  
            h = (r**2 + d**2 -2*d*r*np.cos(theta))**0.5
            if h >= rho:
                tumor[i,j] = -0.05

    # Calculate energy deposited inside tumor
    reduced_master = master_array.copy()
    for i,r in enumerate(rvals):
        for j,theta in enumerate(thetavals):
            if tumor[i,j] < 0:
                reduced_master[i,j] = 0
    num_integral_tumor = simps(simps(reduced_master, thetavals), rvals)

    return num_integral, num_integral_tumor
        
# CALL ABOVE FUNCTION IN A LOOP
    
svals = np.arange(0.1,0.35,0.05)
Evals = np.arange(85,107,2)

E_total = np.zeros([len(svals),len(Evals)])
E_tumor = np.zeros([len(svals),len(Evals)])
array_stats = {'svals':svals,'Evals':Evals}
f = open(r'../Python Output/array_stats.data','wb')
pickle.dump(array_stats,f)
f.close()


TIME0 = time.time()

for i,s in enumerate(svals):
    for j,E_val in enumerate(Evals):
        energy_total, energy_tumor = master(E_val,s)
        
        E_total[i,j] = energy_total
        E_tumor[i,j] = energy_tumor
        array_stats[i,j] = (s,E_val)
        
        np.savetxt(r'../Python Output/2D_Energy_total_array.txt',E_total)
        np.savetxt(r'../Python Output/2D_Energy_tumor_array.txt',E_tumor)
        
        TIME_current = time.time()
        print("Combination {0} of {1} finished after {2:.2f} minutes elapsed."\
              .format(i*len(Evals)+j+1,len(svals)*len(Evals),(TIME_current-TIME0)/60))
        
        