"   Project for PGF5339: Mini-course on structure formation at large-scales "
"   Carolina Queiroz & Maria Benito "
"   March - May, 2017   "

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from scipy import interpolate
import sys

#-------------------------------------------------------------------------------
#               STATISTICAL UNCERTAINTIES IN RATE P0/P2
# (1) Check that rates(k) are normally distribute for number_simu simulations
#-------------------------------------------------------------------------------
# Number of k's
number_k    = 49
# Number of simulations
number_simu = 400
#  k values
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

# Load simulated spectrums
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
                print Data
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
                print Data
            i += 1
            
# Rate P0/P4
f       = 0.55
# Bias
b       = 1.
beta    = f/b
# Theoretic rate
rate_th = np.ones(kbin.shape)*(1+(2/3.)*beta+(1/5.)*beta**2)/((4/3.)*beta+(4/7.)*beta**2)
# Simulated rate
rate = np.zeros((number_k, number_simu))
for k in range(number_k):
    for simu in range(number_simu):
        rate[k][simu] = p0[k][simu] / p2[k][simu]

# Rate mean and SD for each k-bin
rate_mean = np.zeros(number_k)
rate_SD   = np.zeros(number_k)
for k in range(number_k):
    mean         = sum(rate[k])/len(rate[k])
    rate_mean[k] = mean
    SD           = np.sqrt(1./(len(rate[k])-1)*sum(np.power(rate[k] - mean, 2)))
    rate_SD[k]   = SD
    
# Gaussian function
def gaussian(x, a, mean, sigma):
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

#-------------------------------------------------------------------------------
#               DO THE RATES FOLLOW A GAUSSIAN DISTRIBUTION?
#-------------------------------------------------------------------------------
name=raw_input("Do you want to visualize if the rates follow a Gaussian distribution? (y/n) ")
if name=='y':
    # Plot histogram and corresponding gaussian for each k-bin           
    for k in range(number_k):
        f = plt.figure(k)  
        ax = f.add_subplot(111) 
        count, bins, ignored = plt.hist(rate[k], 50, color='b')
        sigma = rate_SD[k]
        mu    = rate_mean[k]
        x = np.linspace(min(bins)-0.2, max(bins)+0.2, 50)
        plt.plot(x, gaussian(x, max(count), mu, sigma), lw=3., c='r')
        plt.xlabel("rate"); plt.ylabel("Counts")
        plt.text(0.85, 0.9, 'k = '+str(kbin[k]), horizontalalignment='center',
                verticalalignment='center', transform = ax.transAxes)
        plt.text(0.1, 0.9, 'N = '+str(number_simu), 
                horizontalalignment='center', verticalalignment='center', 
                transform = ax.transAxes)
        plt.tick_params(axis='both', labelsize=14)
        plt.rc('font', family='helvetica', size=16)
        #plt.savefig('N_'+str(number_simu)+'_k_'+str(k)+'.png', bbox_inches='tight')

#-------------------------------------------------------------------------------
#                              PLOT RATES P0/P2
#-------------------------------------------------------------------------------       
name=raw_input("Do you want to visualize the rates P0/P2? (y/n) ")
if name=='y':
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    name=raw_input("Do you want to visualize all the simulated rates P0/P2? (y/n) ")
    if name=='y':
        for k in range(number_k):
            kbin_simu = np.ones(rate[k].shape)*kbin[k]
            # plot all simulation rates
            plt.plot(kbin_simu, rate[k], 'x', c='#2e7d32')
            plt.ylim([-2,6])
            plt.xlabel(r'$ k \,\,{\rm [h \,Mcp^{-1}]} $'); plt.ylabel(r'$P_0/P_2$')
            plt.tick_params(axis='both', labelsize=14)
            plt.rc('font', family='helvetica', size=22)
    plt.plot(kbin, rate_th, lw=2.5, ls='--', c='#d35400', label='theory')
    plt.errorbar(kbin, rate_mean, yerr=rate_SD, fmt='o', lw=2., c='#6a1b9a',
             capsize=2., capthick=2., elinewidth=2.)
    plt.grid()
    plt.legend(frameon=False)
    #plt.savefig('N_'+str(number_simu)+'_rate_th_uncertainties_2.png', bbox_inches='tight')
    plt.show()
    
#-------------------------------------------------------------------------------
#                              WRITE RATES P0/P2
#-------------------------------------------------------------------------------       
name=raw_input("Do you want to write all the simulated rates P0/P2 into file? (y/n) ")
if name=='y':
    out = open('rate_f_055.dat','w')
    
    for simu in range(number_simu):
        for k in range(number_k):
            out.write('%.5f   '%(kbin[k]))
            out.write('%.5f   '%(rate_mean[k]))
            out.write('%.5f'%(rate_SD[k]))
            out.write('\n')

    #kfile.close()
    out.close()

#-------------------------------------------------------------------------------
#                              PLOT RESIDUALS
#-------------------------------------------------------------------------------
name=raw_input("Do you want to plot residuals for the rates P0/P2? (y/n) ")
if name=='y':
    residuals = rate_th - rate_mean
    fig = plt.figure()
    ax  = fig.add_subplot(111)

    plt.errorbar(kbin, residuals, yerr=rate_SD, fmt='.', c='#6a1b9a',
                 capsize=5., capthick=3., elinewidth=3.)
    plt.plot(kbin, np.zeros(kbin.shape), lw=2., ls='--', alpha=0.5, c='#17202a')
    plt.ylim([-2,2])         
    plt.xlabel(r'$ k \,\,{\rm [h \,Mcp^{-1}]} $'); plt.ylabel(r'${\rm residuals}$')
    plt.grid()
    plt.tick_params(axis='both', labelsize=14)
    plt.rc('font', family='helvetica', size=22)  
    #plt.savefig('N_'+str(number_simu)+'_residuals.png', bbox_inches='tight')
    plt.show()
    
#-------------------------------------------------------------------------------
#                       UNCERTAINTIES IN RATE P0/P2
#-------------------------------------------------------------------------------
name=raw_input("Do you want to write uncertainites in rate P0/P2 into file? (y/n) ")
if name=='y':
    # Open file
    out = open('/Users/maria/Dropbox/camb/data/rate_uncertainties.dat','w')
    
    # Write into file
    for i in range(number_k):
        out.write('%.5f   '%(kbin[i]) )
        out.write('%.5f'%(rate_SD[i]))
        out.write('\n')

    #for i in range(Nbin-1):
    #    out.write('%.5f '%(p2_kbin_1[i]))

    #out.write('\n')
    
    # Close file
    out.close()

    

#-------------------------------------------------------------------------------
#                              SCAN GRID IN f
#-------------------------------------------------------------------------------
# Compute uncertainties on f       
# 1 sigma chi2
sigma       = 1
ci          = scipy.stats.chi2.cdf(sigma**2, 1)
chi2_1sigma = scipy.stats.chi2.ppf(ci, 1)

# Grid in f
f_grid    = np.arange(0.1, 1.5, 0.001)
beta_grid = f_grid/b
rate_grid = (1+(2/3.)*beta_grid+(1/5.)*np.power(beta_grid,2))/((4/3.)*beta_grid+(4/7.)*np.power(beta_grid,2))

# chi2
chi       = np.zeros(f_grid.shape)
# minimum chi
chi_min = 100.
f_min   = 200000.
# uncertainty_chi
chi_sigma = []
f_sigma   = []

# Scan grid in f
for i in range(f_grid.shape[0]):
    rth    = rate_grid[i]
    chi[i] = sum(np.power(rate_mean - rth,2)/np.power(rate_SD, 2))
    if chi[i] < chi_min:
        f_min   = f_grid[i]
        chi_min = chi[i]
print "Best Fit f is ", f_min, ", with a chi2/dof of ", chi_min/(len(kbin)-1)
pos_chi = np.where(abs(chi-chi_min) <= chi2_1sigma)
print "+sigma_f = ", max(f_grid[pos_chi])-f_min
print "-sigma_f = ", f_min-min(f_grid[pos_chi]) 
               
#plt.figure()
#plt.plot(f_grid, chi, lw=3., c='#d4ac0d')
#plt.ylim([0,30])
#plt.show()

#-------------------------------------------------------------------------------
#                  PLOT THEORETIC AND ESTIMATED RATES
#-------------------------------------------------------------------------------
# Estimated rate
beta_min = f_min/b
rate_min = (1+(2/3.)*beta_min+(1/5.)*beta_min**2)/((4/3.)*beta_min+(4/7.)*beta_min**2)
rate_est = np.ones(kbin.shape)*rate_min

fig = plt.figure()
ax  = fig.add_subplot(111)
plt.plot(kbin, rate_th, lw=2.5, ls = '--', c='#d35400', label="theory")
plt.plot(kbin, rate_est, lw=2.5, c='#d35400', label = "best fit")
plt.errorbar(kbin, rate_mean, yerr=rate_SD, fmt='o', lw=2., c='#6a1b9a',
             capsize=2., capthick=2., elinewidth=2.)
plt.grid()
plt.ylim([-2,6])
plt.xlabel(r'$ k \,\,{\rm [h \,Mcp^{-1}]} $'); plt.ylabel(r'$P_0/P_2$')
plt.tick_params(axis='both', labelsize=14)
plt.rc('font', family='helvetica', size=18) 
plt.legend(frameon=False) 
#plt.savefig('N_'+str(number_simu)+'_rate_th_bf.png', bbox_inches='tight')       
plt.show()


#-------------------------------------------------------------------------------
#             COMPUTE BIAS A IN QUADRUPOLE, i.e P_2est = A*P_2th
#-------------------------------------------------------------------------------
# Plot quadrupole
p2_th  = Pth*((4/3.)*beta+(4/7.)*beta**2)
p2_est = Pth*((4/3.)*beta_min+(4/7.)*beta_min**2)

fig = plt.figure()
ax  = fig.add_subplot(111)
for k in range(number_k):
    kbin_simu = np.ones(rate[k].shape)*kbin[k]
    plt.plot(kbin_simu, p2[k], '.', c='g')
    #plt.ylim([-2,6])
plt.plot(kth, p2_th, lw=2., c='#17202a', label='theory')
plt.plot(kth, p2_est, lw=2., c='#6c3483', label='best fit')
plt.xlabel(r'$ k \,\,{\rm [h \,Mcp^{-1}]} $'); plt.ylabel(r'$P_2\,\,[h^{-3}\,Mpc^{3}]$')
plt.tick_params(axis='both', labelsize=14)
plt.rc('font', family='helvetica', size=22)
plt.plot(kth, p2_th, lw=2., c='#17202a')
plt.legend(frameon=False, loc = 3) 
plt.grid()
plt.text(0.85, 0.9, 'N = '+str(number_simu), 
         horizontalalignment='center', verticalalignment='center', 
         transform = ax.transAxes)
plt.xlim([0.001, 0.3])
plt.ylim([100, 100000])         
plt.yscale('log')
plt.xscale('log')    
#plt.savefig('N_'+str(number_simu)+'_quadrupole_th_est.png',bbox_inches='tight')

# Interpolate CAMB spectrum
Pm = interpolate.InterpolatedUnivariateSpline(Pk[0][:],Pk[1][:],k=3) 

# Redefine quadrupole at kbin
p2_th  = Pm(kbin)* ((4/3.)*beta+(4/7.)*beta**2) 

# P_2 mean and uncertainty
p2_mean = np.zeros(number_k)
p2_SD   = np.zeros(number_k)
for k in range(number_k):
    mean       = sum(p2[k])/len(p2[k])
    p2_mean[k] = mean
    SD         = np.sqrt(1./(len(p2[k])-1)*sum(np.power(p2[k] - mean, 2)))
    p2_SD[k]   = SD

A_grid  = np.arange(0.1, 1., 0.001)
chi     = np.zeros(A_grid.shape)
chi_min = 100.
A_min   = 20000.

for i in range(A_grid.shape[0]):
    chi[i] = np.sum(np.power(A_grid[i]*p2_th - p2_mean, 2)/np.power(p2_SD, 2))
    if chi[i] < chi_min:
        A_min   = A_grid[i]
        chi_min = chi[i]
print "Best Fit A is ", A_min, ", with a chi2/dof of ", chi_min/(len(kbin)-1)
pos_chi = np.where(abs(chi-chi_min) <= chi2_1sigma)
print "+sigma_A = ", max(A_grid[pos_chi])-A_min
print "-sigma_A = ", A_min-min(A_grid[pos_chi])


# Plot quadrupole
ax  = fig.add_subplot(111)
for k in range(number_k):
    kbin_simu = np.ones(number_simu)*kbin[k]
    plt.plot(kbin_simu, p2[k], '.', c='#e65100')    
plt.plot(kbin, p2_th, lw=2., c='#e65100', label=r'$P_2$')
plt.errorbar(kbin, p2_mean, yerr=p2_SD, fmt='o', c='#3e2723',
             capsize=2., capthick=2., elinewidth=2.)
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
#plt.savefig('quadrupole_bias.png', bbox_inches='tight')
plt.show()

#-------------------------------------------------------------------------------
#                     ESTIMATE RATE WITH BIAS
#-------------------------------------------------------------------------------
# Simulated rate
rate_bias = np.zeros((number_k, number_simu))
for k in range(number_k):
    for simu in range(number_simu):
        rate_bias[k][simu] = A_min*p0[k][simu] / p2[k][simu]

# Mean and SD for each k-bin
rate_mean_bias = np.zeros(number_k)
rate_SD_bias   = np.zeros(number_k)
for k in range(number_k):
    mean              = sum(rate_bias[k])/len(rate_bias[k])
    rate_mean_bias[k] = mean
    SD                = np.sqrt(1./(len(rate_bias[k])-1)*
                        sum(np.power(rate_bias[k] - mean, 2)))
    rate_SD_bias[k]   = SD
    
# chi2
chi       = np.zeros(f_grid.shape)
# minimum chi
chi_min = 100.
f_min   = 200000.
# uncertainty_chi
chi_sigma = []
f_sigma   = []

# Scan grid in f
for i in range(f_grid.shape[0]):
    rth    = rate_grid[i]
    chi[i] = sum(np.power(rate_mean_bias - rth,2)/np.power(rate_SD_bias, 2))
    if chi[i] < chi_min:
        f_min   = f_grid[i]
        chi_min = chi[i]
print "Best Fit f is ", f_min, ", with a chi2/dof of ", chi_min/(len(kbin)-1)
pos_chi = np.where(abs(chi-chi_min) <= chi2_1sigma)
print "+sigma_f = ", max(f_grid[pos_chi])-f_min
print "-sigma_f = ", f_min-min(f_grid[pos_chi]) 

# Estimated rate
beta_min = f_min/b
rate_min_bias = (1+(2/3.)*beta_min+(1/5.)*beta_min**2)/((4/3.)*beta_min+(4/7.)*beta_min**2)
rate_est = np.ones(kbin.shape)*rate_min
rate_est_bias = np.ones(kbin.shape)*rate_min_bias

fig = plt.figure()
ax  = fig.add_subplot(111)
plt.plot(kbin, rate_th, lw=2.5, ls = '--', c='#d35400', label="theory")
plt.plot(kbin, rate_est_bias, lw=2.5, c='#d35400', label = "best fit")
plt.errorbar(kbin, A_min*rate_mean, yerr=rate_SD, fmt='o', lw=2., c='#6a1b9a',
             capsize=2., capthick=2., elinewidth=2.)
plt.grid()
plt.ylim([-2,6])
plt.xlabel(r'$ k \,\,{\rm [h \,Mcp^{-1}]} $'); plt.ylabel(r'$P_0/P_2$')
plt.tick_params(axis='both', labelsize=14)
plt.rc('font', family='helvetica', size=18) 
plt.legend(frameon=False) 
#plt.savefig('N_'+str(number_simu)+'_rate_th_bf.png', bbox_inches='tight')       
plt.show()