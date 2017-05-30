"""
Pipeline for obtaining "oberved" uncertainties on rate P0/P4
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from scipy import interpolate
import sys
"""
First, we need to check that the rates are normally distributed.
We plot for each k, an histogramn for the rates.
"""

# Number of k's
number_k    = 49
# Number of simulations
number_simu = 1
#  Load k values and rate uncertainties
datafile  = 'rate_uncertainties.dat'
test      = np.loadtxt(datafile,unpack=True)
kbin      = test[0]
rate_SD   = test[1]

# Theoretic rate for MOG
f0      = 0.5
w       = 1./0.15
b       = 1.

# Load CAMB spectrum
datafile = 'cosmo_project_matterpower.dat'
Pk  = np.loadtxt(datafile,unpack=True)
kth = Pk[0]
Pth = Pk[1]

f_th    = f0*np.sqrt(np.power(np.sin(2*np.pi*kth*w), 2))
beta_th = f_th/b
rate_th = (1+(2/3.)*beta_th+(1/5.)*beta_th**2)/((4/3.)*beta_th+(4/7.)*beta_th**2)

p0_th = Pk[1] * (1+(2/3.)*beta_th+(1/5.)*beta_th**2)
p2_th = Pk[1] * ((4/3.)*beta_th+(4/7.)*beta_th**2)

# Load simulated spectrums
with open('/Users/maria/Dropbox/camb/data/Proj_multipoles_MOG.dat', 'r') \
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

rate_mean = np.zeros(number_k)

for k in range(number_k):
    rate_mean[k] = p0[k][0]/p2[k][0]


#-------------------------------------------------------------------------------
#                          PLOT RATES
#-------------------------------------------------------------------------------       
fig = plt.figure()
ax  = fig.add_subplot(111)
#plt.plot(kbin, rate_th, lw=2.5, ls='--', c='#d35400', label='theory')
plt.plot(kth, rate_th, lw=2.5, ls='--', c='g', label='theory MOG')
plt.plot(kth, f_th, lw=2.5, ls='-.', c='g', label='f MOG')
plt.errorbar(kbin, rate_mean, yerr=rate_SD, fmt='o', lw=2., c='#6a1b9a',
             capsize=2., capthick=2., elinewidth=2.)
plt.xlim([0.0, 0.25])
plt.ylim([-1., 25])
plt.xlabel(r'$ k \,\,{\rm [h \,Mcp^{-1}]} $'); plt.ylabel(r'$P_0/P_2$')
plt.tick_params(axis='both', labelsize=14)
plt.rc('font', family='helvetica', size=22)             
plt.grid()
plt.legend(frameon=False)
#plt.savefig('N_'+str(number_simu)+'_rate_MOG.png', bbox_inches='tight')
plt.show()

#-------------------------------------------------------------------------------
# Plot theoretic and simulated monopole and quadrupole. Just for one simulation.
#-------------------------------------------------------------------------------
p0_est = np.zeros(number_k)
p2_est = np.zeros(number_k)
for i in range(number_k):
    p0_est[i] = p0[i][0]
    p2_est[i] = p2[i][0]
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


#-------------------------------------------------------------------------------
#                              SCAN GRID IN f
#-------------------------------------------------------------------------------
# Compute uncertainties on f       
# 1 sigma chi2 for 1 dof
sigma       = 1
ci          = scipy.stats.chi2.cdf(sigma**2, 1)
chi2_1sigma = scipy.stats.chi2.ppf(ci, 1)

# Bias
b = 1.

# Grid in f
f_0_grid = np.arange(0.1, 1.5, 0.01)
w_grid   = np.arange(0.1, 10, 0.01)

F_0, W   = np.meshgrid(f_0_grid, w_grid)

F    = np.zeros((number_k, F_0.shape[0], F_0.shape[1]))
BETA = np.zeros((number_k, F_0.shape[0], F_0.shape[1]))
RATE = np.zeros((number_k, F_0.shape[0], F_0.shape[1]))
for k in range(number_k):
    F[k] = F_0*np.sqrt(np.power(np.sin(2*np.pi*kbin[k]*W), 2))
    BETA[k] = F[k]/b
    RATE[k] = (1+(2/3.)*BETA[k]+(1/5.)*np.power(BETA[k],2))/((4/3.)*BETA[k]+(4/7.)*np.power(BETA[k],2))

# chi2
chi       = np.zeros((F_0.shape[0], F_0.shape[1]))
# minimum chi
chi_min = 1e10
f_0_min = 200000.
w_min   = 2000000.
# uncertainty_chi
chi_sigma = []
f_sigma   = []

# Scan grid in f
for i in range(F_0.shape[0]):
    for j in range(F_0.shape[1]):
        rth = np.zeros(number_k)
        for k in range(number_k):
            rth[k]    = RATE[k][i][j]
        chi[i][j] = sum(np.power(rate_mean - rth,2)/np.power(rate_SD, 2))
        if chi[i][j] < chi_min:
            f_0_min = f_0_grid[j]
            w_min   = w_grid[i]
            chi_min = chi[i][j]
        
        
print "Best Fit f_0  is ", f_0_min, ", with a chi2/dof of ", chi_min/(len(kbin)-1)
print "Best Fit 1./w is ", w_min
pos_chi = np.where(abs(chi-chi_min) <= chi2_1sigma)

sys.exit(-1)
print "+sigma_f_0 = ", max(F_0[pos_chi])-f_0_min
print "-sigma_f_0 = ", f_0_min-min(F_0[pos_chi]) 


beta_min = f_0_min*np.sqrt(np.power(np.sin(2*np.pi*kth*w_min), 2))
rate_min =(1+(2/3.)*beta_min+(1/5.)*beta_min**2)/((4/3.)*beta_min+(4/7.)*beta_min**2)

#-------------------------------------------------------------------------------
#                          PLOT RATES
#-------------------------------------------------------------------------------       
fig = plt.figure()
ax  = fig.add_subplot(111)
#plt.plot(kbin, rate_th, lw=2.5, ls='--', c='#d35400', label='theory')
plt.plot(kth, rate_th, lw=2.5, ls='--', c='g', label='theory')
plt.plot(kth, rate_min, lw=2.5, ls='-', c='#d35400', label='best fit')
plt.errorbar(kbin, rate_mean, yerr=rate_SD, fmt='o', lw=2., c='#6a1b9a',
             capsize=2., capthick=2., elinewidth=2.)
plt.xlim([0.0, 0.25])
plt.ylim([-1., 25])
plt.xlabel(r'$ k \,\,{\rm [h \,Mcp^{-1}]} $'); plt.ylabel(r'$P_0/P_2$')
plt.tick_params(axis='both', labelsize=14)
plt.rc('font', family='helvetica', size=22)             
plt.grid()
plt.legend(frameon=False)
#plt.savefig('N_'+str(number_simu)+'_rate_MOG.png', bbox_inches='tight')
plt.show()

