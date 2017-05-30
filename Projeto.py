# -*- coding: utf-8 -*-
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import sys
import time


start=time.time()

### Size of the box
nx = 256
ny = 256
nz = 256
# Cell size in units of Mpc.h^-1
l  = 10.         
# Length each side box
Lx = nx*l
Ly = ny*l
Lz = nz*l

# Grid 3d    
cub = np.zeros((nx,ny,nz))
# Volume
box_vol = Lx*Ly*Lz

# Normalization   
name=raw_input("Do you want to use any normalization? (y/n) ")
if name=='y':
    Norm1=box_vol
    Norm2=l**3
elif name=='n':
    Norm1=1.
    Norm2=1.
Norm3=1.

# Load theoretical spectrum 
datafile = '/Users/maria/Desktop/structure formation/code/cosmo_project_matterpower.dat'
Pk = np.loadtxt(datafile,unpack=True)
# Interpolate theoretical spectrum
Pm = interpolate.InterpolatedUnivariateSpline(Pk[0][:],Pk[1][:],k=3)  

### Def frequency 
ni_x = np.fft.fftfreq(nx)
idx  = np.ones(nx)
ni_y = np.fft.fftfreq(ny)
idy  = np.ones(ny)
ni_z = np.fft.fftfreq(nz)
idz  = np.ones(nz)

# Physical k: k=2 * pi/cell_size * frequency
k_x = 2*np.pi*ni_x/l
k_y = 2*np.pi*ni_y/l
k_z = 2*np.pi*ni_z/l

# Einstein sum
q_ijk_x  = np.einsum( 'i,j,k', k_x , idy , idz )
q_ijk_y  = np.einsum( 'i,j,k', idx , k_y , idz )
q_ijk_z  = np.einsum( 'i,j,k', idx , idy , k_z )

q2_ijk_x = np.einsum( 'i,j,k', k_x*k_x , idy , idz )
q2_ijk_y = np.einsum( 'i,j,k', idx , k_y*k_y , idz )
q2_ijk_z = np.einsum( 'i,j,k', idx , idy , k_z*k_z )

k_abs=np.sqrt( (q2_ijk_x) + (q2_ijk_y) + (q2_ijk_z) )
#-------------------------------------------------------------------------------
#                         Steps (1), (2), (3)
#-------------------------------------------------------------------------------
### Log-normal map 
dk   = 0.0005
dr   = 0.5
veck = np.arange(np.min(Pk[0][:]),np.max(Pk[0][:]),dk)
vecr = np.arange(0.00000001,1000,dr)

kr   = np.outer(veck,vecr)
j0   = np.sin(kr)/(0.00000000001+kr)
integrand = j0.T*Pm(veck)*veck**2

print 'Calculating csi(r)...'

csi_ln = (1./(2*np.pi**2))*(dk)*np.sum(integrand,axis=1)

# Define RSD model (f, bias, etc.)
f = 0.5
b = 1.
beta = f/b

mono_norm = 1 + 2./3.*beta + 1./5.*beta**2

b_bar = b*np.sqrt(mono_norm)

### csi_gaussian 
csi_g = np.log(1 + b_bar**2 * csi_ln)

integrand2 = j0*csi_g*vecr**2

print 'Calculating pG(k)...'

pG = (4*np.pi)*dr*np.sum(integrand2,axis=1)

pG_interp = interpolate.InterpolatedUnivariateSpline(veck,pG,k=3)
#-------------------------------------------------------------------------------
#                           Def Gaussian modes  
#-------------------------------------------------------------------------------
# Setting the min_kabs equal to min_kcamb: Can I do this?!
k_abs[np.where(k_abs==0)] = np.min(Pk[0][:])  
        
# NO SERÍA AQUI DONDE INTRODUCIMOS LAS RSD???!
mu_k = np.ndarray.flatten(q_ijk_z/(k_abs+0.00000000001))


k_abs  = k_abs.flatten()
pG_kmu = (1 + beta*np.power(mu_k, 2))**2/mono_norm * pG_interp(k_abs)

#sigma_ijk=np.sqrt( 2*pG_interp(k_abs)*Norm1 )
sigma_ijk=np.sqrt(2*pG_kmu*Norm1)

k_abs=k_abs.reshape(nx,ny,nz)

# AQUI INTRODUZIR AS RSDs!
sigma_ijk=sigma_ijk.reshape(nx,ny,nz)

A   = np.random.normal(0.0,sigma_ijk)
phi = 2*np.pi*np.random.random(np.shape(sigma_ijk))
delta_ijk_til = A*np.exp(1j*phi)

###   Def Fourier transform   
delta_real = (1./Norm2)*(np.fft.ifftn(delta_ijk_til)).real

print 'ok!'
#-------------------------------------------------------------------------------
#                                Def delta_LN 
#-------------------------------------------------------------------------------
varg = np.var(delta_real)
delta_ln1 = np.exp(delta_real-(varg)/2.) - 1.0

delta_ijk_til2 = np.fft.fftn(delta_ln1)
delta_m = delta_ijk_til2.flatten()
#---------------------------------------------------------------------------
#                Introduce shot noise (Poisson counts)
#---------------------------------------------------------------------------    
Nmean1=1.0

Npoisson1 = np.random.poisson(Nmean1*(1+delta_ln1))

delta_ln1 = (Npoisson1-Nmean1)/Nmean1

delta_ijk_til3_1 = np.fft.fftn(delta_ln1)

delta_m_1=delta_ijk_til3_1.flatten()

k_m = k_abs.flatten()

Nbin=50
#kmin=np.min(k_m)
kmin = 5.*np.min(k_m)
kmax = 0.2 #np.max(k_m)
k_bin=np.linspace(kmin,kmax,Nbin)
k_bin_mean=[]

p_bin_1=[]

for i in range(Nbin-1):
    k_bin_mean.append( (k_bin[i]+k_bin[i+1])/2. )
    pos_k_m=np.where((k_m>=k_bin[i]) & (k_m<k_bin[i+1]))[0]
    p_bin_1.append(1./(Norm3*len(pos_k_m))*np.sum(np.abs(delta_m_1[pos_k_m])**2))

k_bin_mean=np.asarray(k_bin_mean)
p_bin_1=np.asarray(p_bin_1)

#Shot1 = np.power(k_bin_mean,0)*1./(Nmean1/Norm2)
Shot1 = 1./(Nmean1/Norm2)

#-------------------------------------------------------------------------------
#                            Define Multipoles
#-------------------------------------------------------------------------------
mu_k  = np.ndarray.flatten(q_ijk_z/(k_abs+0.00000000001))
abs_mu_k = np.abs(mu_k)

mu_bin = np.linspace(0,1,10)
mu_bin_mean = np.zeros(len(mu_bin)-1)
p_kmubin_1=np.zeros((Nbin-1,len(mu_bin)-1))

for i in range(Nbin-1):    
    for j in range(len(mu_bin)-1):
        mu_bin_mean[j] = mu_bin[j] + (mu_bin[j+1]-mu_bin[j])/2.
        pos_kmu_m = np.where( (k_m>=k_bin[i]) & (k_m<k_bin[i+1]) 
                            & (abs_mu_k >= mu_bin[j]) & (abs_mu_k < mu_bin[j+1]) )[0]
        p_kmubin_1[i,j] = np.mean(np.abs(delta_m_1[pos_kmu_m])**2)
        #print i,"  ", j, "  ", pos_kmu_m.shape

print("Llegamos hasta aqui")

# Normalize and subtract shot noise
p_est_kmu = p_kmubin_1/(Norm1/Norm2**2)

# Step in mu
step_mu = mu_bin[1]-mu_bin[0]

# Monopole
p0_kbin_1 = np.zeros(Nbin-1)

# Polynome de Legendre
Legendre2 = -0.5+1.5*np.power(mu_bin_mean,2)

# Quadrupole
p2_kbin_1 = np.zeros(Nbin-1)

for i in range(Nbin-1):
    p0_kbin_1[i] = np.sum((p_est_kmu[i,:] - Shot1)) * step_mu
    p2_kbin_1[i] = 5*np.sum( (p_est_kmu[i,:] - Shot1) *Legendre2) * step_mu
   

beta = f

plt.figure()    
plt.plot(k_bin_mean,Pm(k_bin_mean)*(1+(2/3.)*beta+(1/5.)*beta**2),'k-')
plt.plot(k_bin_mean,Pm(k_bin_mean)*((4/3.)*beta+(4/7.)*beta**2),'b-')
plt.plot(k_bin_mean,p0_kbin_1,'k.')
plt.plot(k_bin_mean,p2_kbin_1,'b.')
plt.xscale('log')
plt.yscale('log')
plt.show()
    
sys.exit(-1)

# Theoretical multipoles
def P0(beta):
    return 1+(2/3.)*beta+(1/5.)*beta**2

def P2(beta):
    return (4/3.)*beta+(4/7.)*beta**2

def B(u,beta):
    return (1+beta*u**2)**2

#####

#####   Redshift space distortions  #####

mu_k=q_ijk_z/(k_abs+0.00000000001)          

#####

end=time.time()

print start-end
