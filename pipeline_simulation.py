"   Project for PGF5339: Mini-course on structure formation at large-scales "
"   Carolina Queiroz & Maria Benito "
"   March - May, 2017   "

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import sys

# Define theoretical variables 
f_th    = 0.55

# Number of k's
number_k    = 49
# Number of simulations
number_simu = 1
#  Load k values and rate uncertainties
datafile  = '/Users/maria/Dropbox/camb/data/rate_uncertainties.dat'
test      = np.loadtxt(datafile,unpack=True)
kbin      = test[0]
rate_SD   = test[1]

# Load simulated spectrum
with open('/Users/maria/Dropbox/camb/data/Proj_multipoles_f_055.dat', 'r') \
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
#                  WRITE ESTIMATED RATES
#-------------------------------------------------------------------------------
name=raw_input("Do you want to write estimated rated into file? (y/n) ")
if name=='y':
    out = open('rate_f_055.dat','w')
    
    for k in range(number_k):
        out.write('%.5f   '%(kbin[k]))
        out.write('%.5f   '%(rate_mean[k]))
        out.write('%.5f'%(rate_SD[k]))
        out.write('\n')

    #kfile.close()
    out.close()

#-------------------------------------------------------------------------------
#                              SCAN GRID IN f
#-------------------------------------------------------------------------------
# Compute uncertainties on f       
# 1 sigma chi2 for 1 dof
sigma       = 1
ci          = scipy.stats.chi2.cdf(sigma**2, 1)
chi2_1sigma = scipy.stats.chi2.ppf(ci, 1)

# Grid in f
f_grid    = np.arange(0.1, 1.5, 0.001)
# Bias
b = 1.
beta_grid = f_grid/b
rate_grid = (1+(2/3.)*beta_grid+(1/5.)*np.power(beta_grid,2))/((4/3.)*beta_grid+(4/7.)*np.power(beta_grid,2))

# chi2
chi       = np.zeros(f_grid.shape)
# minimum chi
chi_min = 1e10
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
               
#-------------------------------------------------------------------------------
#                  PLOT THEORETIC AND ESTIMATED RATES
#-------------------------------------------------------------------------------
beta_th = f_th/b
rate_th = (1+(2/3.)*beta_th+(1/5.)*beta_th**2)/((4/3.)*beta_th+(4/7.)*beta_th**2)
rate_th = np.ones(number_k)*rate_th

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
plt.legend(frameon=False, loc=4) 
#plt.savefig('N_'+str(number_simu)+'_rate_th_bf.png', bbox_inches='tight')       
plt.show()