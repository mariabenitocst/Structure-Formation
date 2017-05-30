"   Project for PGF5339: Mini-course on structure formation at large-scales "
"   Carolina Queiroz & Maria Benito "
"   March - May, 2017   "
"   Last update: RSDs   "

# -*- coding: utf-8 -*-
import numpy as np
from scipy import interpolate
#from scipy.integrate import quad
#from scipy.integrate import simps
#import matplotlib.pyplot as plt
#import sys
import time

start=time.time()

Nsimu = 1

#kfile = open('Proj_kbins.dat','w')
out = open('Proj_multipoles.dat','w')

for i in range(Nsimu):

 print 'Simu: %d'%(i+1)

#---------------------------------------------------------------------------#

#####   Size of the box

 nx = 256
 ny = 256
 nz = 256

#####   Cell size in units of Mpc.h^-1

 l  = 10.

#####   Length each side box

 Lx = nx*l
 Ly = ny*l
 Lz = nz*l

#####   Grid 3D

 cub = np.zeros((nx,ny,nz))

#####   Volume

 box_vol = Lx*Ly*Lz
 cell_vol = l**3

#####   Normalization

 Norm1=box_vol
 Norm2=cell_vol

#####   Load theoretical spectrum (CAMB)

# [k]=h.Mpc^-1 ; [P]=h^-3.Mpc^3

 datafile = '/Users/maria/Desktop/structure formation/code/data/cosmo_project_matterpower.dat'
 Pk = np.loadtxt(datafile,unpack=True)

# Interpolate theoretical spectrum

 Pm = interpolate.InterpolatedUnivariateSpline(Pk[0][:],Pk[1][:],k=3)

#####   Def frequency

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

#---------------------------------------------------------------------------
#                         Steps (i), (ii), (iii)
#---------------------------------------------------------------------------

 f = 0.55         #matter growth rate
 b = 1.          #bias
 beta = f/b

 mono_norm = 1 + (2./3.)*beta + (1./5.)*beta**2
 b_bar = b*np.sqrt(mono_norm)
 
 dk   = 0.0005
 dr   = 0.5
 veck = np.arange(np.min(Pk[0][:]),np.max(Pk[0][:]),dk)
 vecr = np.arange(0.00000001,1000,dr)

 kr   = np.outer(veck,vecr)
 j0   = np.sin(kr)/(0.00000000001+kr)
 integrand = j0.T*Pm(veck)*veck**2

 csi_ln = (1./(2*np.pi**2))*(dk)*np.sum(integrand,axis=1)

 csi_g = np.log(1 + b_bar**2 * csi_ln)
 #csi_g = np.log(1 + csi_ln)

 integrand2 = j0*csi_g*vecr**2
 pG = (4*np.pi)*dr*np.sum(integrand2,axis=1)

 pG_interp = interpolate.InterpolatedUnivariateSpline(veck,pG,k=3)

#---------------------------------------------------------------------------
#                           Def Gaussian modes  
#---------------------------------------------------------------------------

# Setting the min_kabs equal to min_kcamb

 k_abs[np.where(k_abs==0)] = np.min(Pk[0][:]) 

#####   Introduce the RSDs

 mu_k = np.ndarray.flatten(q_ijk_z/(k_abs+0.00000000001))
 #mu2_k = np.ndarray.flatten(q2_ijk_z/(k_abs+0.00000000001)**2)

 k_abs  = k_abs.flatten()
 
 #f0      = 0.5
 #w       = 1./0.15
 #b       = 1.
 #f       = f0*np.sqrt(np.power(np.sin(2*np.pi*k_abs*w), 2))
 #beta    = f/b
 
 pG_kmu = ((1 + beta*mu_k*mu_k)**2)/mono_norm * pG_interp(k_abs)
 #pG_kmu = ((1 + beta*mu_k*mu_k)**2) * pG_interp(k_abs)

 sigma_ijk = np.sqrt( 2*pG_kmu*Norm1 )

 k_abs=k_abs.reshape(nx,ny,nz)
 sigma_ijk=sigma_ijk.reshape(nx,ny,nz)

 A = np.random.normal(0.0,sigma_ijk)
 phi = 2*np.pi*np.random.random(np.shape(sigma_ijk))
 delta_ijk_til = A*np.exp(1j*phi)

#####   Def Fourier transform

 delta_real = (1./Norm2)*(np.fft.ifftn(delta_ijk_til)).real

#-------------------------------------------------------------------------------
#                                Def delta_LN 
#-------------------------------------------------------------------------------

 varg = np.var(delta_real)
 delta_ln1 = np.exp(b*delta_real-(b**2 * varg)/2.) - 1.0

 delta_ijk_til2 = np.fft.fftn(delta_ln1)
 delta_m = delta_ijk_til2.flatten()

#---------------------------------------------------------------------------
#                Introduce shot noise (Poisson counts)
#---------------------------------------------------------------------------

 Nmean1=1.0
 Npoisson1 = np.random.poisson(Nmean1*(1+delta_ln1))
 delta_ln1 = (Npoisson1-Nmean1)/Nmean1
 delta_ijk_til_3 = np.fft.fftn(delta_ln1)
 delta_m_1=delta_ijk_til_3.flatten()
 k_m = k_abs.flatten()

 Nbin=50
 kmin = 5.0*np.min(k_m)
 kmax = 0.2
 k_bin=np.linspace(kmin,kmax,Nbin)
 k_bin_mean=[]
 p_bin_1=[]

 for i in range(Nbin-1):
    k_bin_mean.append( (k_bin[i]+k_bin[i+1])/2. )
    pos_k_m=np.where((k_m>=k_bin[i]) & (k_m<k_bin[i+1]))[0]
    p_bin_1.append(np.mean(np.abs(delta_m_1[pos_k_m])**2))

 k_bin_mean=np.asarray(k_bin_mean)
 p_bin_1=np.asarray(p_bin_1)

# for i in range(Nbin-1):
#    kfile.write('%.5f   '%(k_bin_mean[i]))
# kfile.write('\n')

 Shot1=np.power(k_bin_mean,0)*1./(Nmean1/Norm2)

#-------------------------------------------------------------------------------
#                Subtract shot noise (Poisson counts)
#-------------------------------------------------------------------------------

 p_bin_1 = p_bin_1/(Norm1/Norm2**2) - Shot1

#-------------------------------------------------------------------------------
#                            Define Multipoles
#-------------------------------------------------------------------------------

 Nmu=9
 abs_mu_k = np.abs(mu_k)

 mu_bin = np.linspace(0,1,Nmu)
 mu_bin_mean = np.zeros(Nmu-1)
 p_kmubin_1=np.zeros((Nbin-1,Nmu-1))

 for i in range(Nbin-1):
    for j in range(Nmu-1):
        mu_bin_mean[j] = mu_bin[j] + (mu_bin[j+1]-mu_bin[j])/2.
        pos_kmu_m = np.where( (k_m>=k_bin[i]) & (k_m<k_bin[i+1]) & (abs_mu_k >= mu_bin[j]) & (abs_mu_k < mu_bin[j+1]) )[0]
        p_kmubin_1[i,j] = np.mean(np.abs(delta_m_1[pos_kmu_m])**2)

# Convert nan's into zeros

 for j in range(Nbin-1):
    for i in range(Nmu-1):
        if np.isnan(p_kmubin_1[j,i]):
            p_kmubin_1[j,i]=0

 step_mu = mu_bin[1]-mu_bin[0]

# Monopole

 p0_kbin_1 = np.zeros(Nbin-1)

# Polynome de Legendre

 Legendre2 = -0.5+1.5*np.power(mu_bin_mean,2)

# Quadrupole

 p2_kbin_1 = np.zeros(Nbin-1)

 for i in range(Nbin-1):
    p0_kbin_1[i] = np.sum(p_kmubin_1[i,:]) * step_mu
    p2_kbin_1[i] = 5*np.sum(p_kmubin_1[i,:]* Legendre2) * step_mu
    
 # Bias in quadrupole
 A = 0.784

 p0_kbin_1 = p0_kbin_1/(Norm1/Norm2**2) - Shot1
 p2_kbin_1 = 1./A*p2_kbin_1/(Norm1/Norm2**2)

 for i in range(Nbin-1):
    out.write('%.5f '%(p0_kbin_1[i]))

 out.write('\n')

 for i in range(Nbin-1):
    out.write('%.5f '%(p2_kbin_1[i]))

 out.write('\n')

#kfile.close()
out.close()

end=time.time()

print start-end
