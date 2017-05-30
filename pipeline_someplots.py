"   Project for PGF5339: Mini-course on structure formation at large-scales "
"   Carolina Queiroz & Maria Benito "
"   March - May, 2017   "

import numpy as np
import matplotlib.pyplot as plt
import sys
"""
First, we need to check that the rates are normally distributed.
We plot for each k, an histogramn for the rates.
"""
# Number of k's
number_k    = 49
# Number of simulations
number_simu = 400

kbin = [0.00254, 0.00661, 0.01068, 0.01475, 0.01882, 0.02289, 0.02696, 0.03104,
        0.03511, 0.03918, 0.04325, 0.04732, 0.05139, 0.05546, 0.05954, 0.06361,   
        0.06768, 0.07175, 0.07582, 0.07989, 0.08396, 0.08804, 0.09211,
        0.09618, 0.10025, 0.10432, 0.10839, 0.11246, 0.11654, 0.12061,
        0.12468, 0.12875, 0.13282, 0.13689, 0.14096, 0.14504, 0.14911,
        0.15318, 0.15725, 0.16132, 0.16539, 0.16946, 0.17354, 0.17761,
        0.18168, 0.18575, 0.18982, 0.19389, 0.19796]
# Convert into np array
kbin = np.asarray(kbin)

# Load CAMB spectrum
datafile = '/Users/maria/Desktop/structure formation/code/data/cosmo_project_matterpower.dat'
Pk  = np.loadtxt(datafile,unpack=True)
kth = Pk[0]
Pth = Pk[1]

# Rate P0/P4
# Theoretic rate
f       = 0.55
# Bias
b       = 1.
beta    = f/b

# Theoretic monopole
p0_th = Pth*(1+(2/3.)*beta+(1/5.)*beta**2)
# Theoretic quadrupole
p2_th = Pth*((4/3.)*beta+(4/7.)*beta**2)

with open('/Users/maria/Dropbox/camb/data/Proj_multipoles_11.dat', 'r') \
as f: 
    # Monopole
    p0 = np.zeros((number_k, number_simu))
    # Quadrupole
    p2 = np.zeros((number_k, number_simu))
    # Count number of simulations
    simu = 0
    # Define if line corresponds to a monopole or a quadrupole
    i    = 0 
    for line in f:
        if not line.startswith("#"):
            Data = line.split()
            try:
                if i%2 == 0:
                    for k in range(number_k):
                        p0[k][simu] = float(Data[k])
                else:
                    for k in range(number_k):
                        p2[k][simu] = float(Data[k])
                    simu += 1
            except IndexError:
                print i, "  ", simu
                #print Data
            i += 1

with open('/Users/maria/Dropbox/camb/data/Proj_multipoles_12.dat', 'r') \
as f: 
    # Define if line corresponds to a monopole or a quadrupole
    i    = 0 
    for line in f:
        if not line.startswith("#"):
            Data = line.split()
            try:
                if i%2 == 0:
                    for k in range(number_k):
                        p0[k][simu] = float(Data[k])
                else:
                    for k in range(number_k):
                        p2[k][simu] = float(Data[k])
                    simu += 1
            except IndexError:
                print i, "  ", simu
                #print Data
            i += 1

"""
Plot theoretic and simulated monopole and quadrupole. Just for one simulation.
"""
p0_est = np.zeros(number_k)
p2_est = np.zeros(number_k)
simu = 25
for i in range(number_k):
    p0_est[i] = p0[i][simu]
    p2_est[i] = p2[i][simu]
plt.figure()
plt.plot(kth, p0_th, lw=2., c='#0d47a1', label=r'$P_0$')
plt.plot(kbin, p0_est, 'o', c='#0d47a1')
plt.plot(kth, p2_th, lw=2., c='#e65100', label=r'$P_2$')
plt.plot(kbin, p2_est, 'o', c='#e65100')
plt.xlim(0.001, 0.3)
plt.ylim(1e3, 1e5)
plt.yscale('log')
plt.xscale('log')
plt.grid(True,which="both",alpha=0.2)
plt.legend(frameon=False, loc=3)
plt.ylabel(r'$P(k)\,\,{\rm [h^{-3}\,Mpc^{3}]}$')
plt.xlabel(r'$k\,\,{\rm [h\,Mpc^{-1}]}$')
plt.show()

"""
Plot theoretic and simulated monopole and quadrupole. For all simulations.
"""
fig = plt.figure()
# Plot monopole
ax  = fig.add_subplot(121)
for k in range(number_k):
    kbin_simu = np.ones(number_simu)*kbin[k]
    plt.plot(kbin_simu, p0[k], '.', c='#0d47a1')    
#plt.tick_params(axis='both', labelsize=14)
#plt.rc('font', family='helvetica', size=22)
plt.plot(kth, p0_th, lw=2., c='#0d47a1', label=r'$P_0$')
plt.xlim(0.001, 0.3)
plt.ylim(1e3, 1e5)
plt.yscale('log')
plt.xscale('log')
plt.grid(True,which="both",alpha=0.2)
plt.ylabel(r'$P(k)\,\,{\rm [h^{-3}\,Mpc^{3}]}$')
plt.xlabel(r'$k\,\,{\rm [h\,Mpc^{-1}]}$')
plt.legend(frameon=False, loc=3)
plt.text(0.85, 0.9, 'N = '+str(number_simu), 
         horizontalalignment='center', verticalalignment='center', 
         transform = ax.transAxes)
# Plot quadrupole
ax  = fig.add_subplot(122)
for k in range(number_k):
    kbin_simu = np.ones(number_simu)*kbin[k]
    plt.plot(kbin_simu, p2[k], '.', c='#e65100')    
plt.plot(kth, p2_th, lw=2., c='#e65100', label=r'$P_2$')
plt.xlim(0.001, 0.3)
plt.ylim(1e3, 1e5)
plt.yscale('log')
plt.xscale('log')
plt.grid(True,which="both",alpha=0.2)
plt.ylabel(r'$P(k)\,\,{\rm [h^{-3}\,Mpc^{3}]}$')
# Set yaxis right position
ax.yaxis.set_label_position('right')
ax.yaxis.tick_right()
plt.xlabel(r'$k\,\,{\rm [h\,Mpc^{-1}]}$')
plt.legend(frameon=False, loc=3)
plt.text(0.85, 0.9, 'N = '+str(number_simu), 
         horizontalalignment='center', verticalalignment='center', 
         transform = ax.transAxes)
plt.show()

"""
Plot simulated rate. For all simulations.
"""
# Simulated rate
rate = np.zeros((number_k, number_simu))
for k in range(number_k):
    for simu in range(number_simu):
        rate[k][simu] = p0[k][simu] / p2[k][simu]
        
fig = plt.figure()
ax  = fig.add_subplot(111)
for k in range(number_k):
    kbin_simu = np.ones(rate[k].shape)*kbin[k]
    plt.plot(kbin_simu, rate[k], '.', c='#2e7d32')
    plt.ylim([-2,6])
    plt.xlabel(r'$ k \,\,{\rm [h \,Mcp^{-1}]} $'); plt.ylabel(r'$P_0/P_2$')
    plt.tick_params(axis='both', labelsize=14)
    plt.rc('font', family='helvetica', size=22)
plt.grid()
plt.text(0.85, 0.9, 'N = '+str(number_simu), 
         horizontalalignment='center', verticalalignment='center', 
         transform = ax.transAxes)
#plt.savefig('N_'+str(number_simu)+'_rate_th_uncertainties_2.png', bbox_inches='tight')
plt.show()
