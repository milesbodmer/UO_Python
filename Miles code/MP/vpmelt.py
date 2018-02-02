#!/opt/local/bin/python2.7
# Test various effective medium estimates
# This script looks at the difference between Berryman's true symmetric
# self-consistent approach and Schmeling's asymmetric self-consistent approach
# The script plots elastic moduli, seismic velocities, attenuation and Vp/Vs
# Michele Paulatto - Imperial College London - January 2018
#
# Licenced under Creative Commons Attribution 4.0 International (CC BY 4.0)
# You are free to copy, use, modify and redistribute this work provided that you 
# give appropriate credit to the original author, provide a link to the licence 
# and indicate if changes were made.
# Full terms at: https://creativecommons.org/licenses/by/4.0/


# Import own libraries
import sca
import schmeling85
#import critical

# Import common libraries
import numpy as np
import matplotlib.pyplot as plt
from math import pi,e,log

###############################################################

# I have set up a couple of options for testing different fluids
fluid='basalt'
# Elastic moduli are in GPa, density is in g/cm^2
if fluid == 'magma':
# Silicic magma
    vpvs = 1.73
    Ki = 1.27e1
    Gi = 1.e-8  # I set this to a very small number > 0 to avoid numerical instabilities
    rof = 2.3
elif fluid == 'water':
# Water
    vpvs = 1.8
    Ki = 2.2
    Gi = 1.e-8  # I set this to a very small number > 0 to avoid numerical instabilities
    rof = 1.05
if fluid == 'basalt':
# Basaltic magma from Mainprice
    Ki = 1.491e1
    Gi = 1.e-8  # I set this to a very small number > 0 to avoid numerical instabilities
    rof = 2.7
    
# Bulk and shear modulus of solid matrix
if fluid in ['magma','water']:
# Intermediate silica content igneous rock
  vp0=6.5
  vs0=vp0/vpvs
  rom=2.7
  Gm = rom*vs0**2
  Km = rom*(vp0**2-4.0/3.0*vs0**2)
elif fluid == 'basalt':
# MOR gabbro from Mainprice
  rom=2.84
  c11=136.03
  c12=56.31
  c44=39.86
  Km=(c11+2*c12)/3.0
  Gm=c44

#print Km,Gm

# Porosity array
dc=0.005
nc=int(1./dc-1)
phi=np.linspace(dc,dc*nc,nc)

# Calculate density of the composite (weighted average)
roe=sca.Voight(rom,rof,phi)
#print roe

#---------------------------------------------
# Bounds
# Reuss
Kr = sca.Reuss(Km,Ki,phi)
Gr = sca.Reuss(Gm,Gi,phi)
Vpr = np.sqrt((Kr+4.0/3.0*Gr)/roe)
Vsr = np.sqrt(Gr/roe)
# Voight
Kv = sca.Voight(Km,Ki,phi)
Gv = sca.Voight(Gm,Gi,phi)
Vpv = np.sqrt((Kv+4.0/3.0*Gv)/roe)
Vsv = np.sqrt(Gv/roe)
# Hashin-Shtrikman
Khs1,Ghs1 = sca.hs_bounds(Km,Ki,Gm,Gi,1.-phi,phi)
Khs2,Ghs2 = sca.hs_bounds(Ki,Km,Gi,Gm,phi,1.-phi)
Vphs1 = np.sqrt((Khs1+4.0/3.0*Ghs1)/roe)
Vshs1 = np.sqrt(Ghs1/roe)
Vphs2 = np.sqrt((Khs2+4.0/3.0*Ghs2)/roe)
Vshs2 = np.sqrt(Ghs2/roe)

#---------------------------------------------
# Use Berryman's method

# Material 1
am=1.0     # Solid grains aspect ratio
ai=1.0     # fluid inclusions aspect ratio
# Unrelaxed, high frequency case
b_Ku1,b_Gu1 = sca.mod_b(Km,Ki,Gm,Gi,am,ai,1.-phi,phi)
#b_Ku1,b_Gu1 = sca.mod_a(Km,Ki,Gm,Gi,ai,phi)
# Relaxed, low frequency case - Gassman
K,G = sca.mod_b(Km,0.0,Gm,0.0,am,ai,1.-phi,phi)
b_Kr1,b_Gr1 = sca.gassman_2(K,G,Km,Gm,Ki,Gi,phi)
#
b_Vpu1=np.sqrt((b_Ku1+4.0/3.0*b_Gu1)/roe)
b_Vsu1=np.sqrt(b_Gu1/roe)
b_Vpr1=np.sqrt((b_Kr1+4.0/3.0*b_Gr1)/roe)
b_Vsr1=np.sqrt(b_Gr1/roe)

# Material 2
am=0.1     # Solid grains aspect ratio
ai=0.1     # fluid inclusions aspect ratio
# Unrelaxed, high frequency
b_Ku2,b_Gu2 = sca.mod_b(Km,Ki,Gm,Gi,am,ai,1.-phi,phi)
#b_Ku2,b_Gu2 = sca.mod_a(Km,Ki,Gm,Gi,ai,phi)
# Relaxed, low frequency - Gassman
K,G = sca.mod_b(Km,0.0,Gm,0.0,am,ai,1.-phi,phi)
b_Kr2,b_Gr2 = sca.gassman_2(K,G,Km,Gm,Ki,Gi,phi)
#
b_Vpu2=np.sqrt((b_Ku2+4.0/3.0*b_Gu2)/roe)
b_Vpr2=np.sqrt((b_Kr2+4.0/3.0*b_Gr2)/roe)
b_Vsu2=np.sqrt(b_Gu2/roe)
b_Vsr2=np.sqrt(b_Gr2/roe)


#---------------------------------------------
# Use Schmeling's method
# Here the solid phase is treated as a special phase that acts as a matrix
# and is not added explicitly. 

# Set up porosity matrix, start with all zeros
# The elements of c[i,j] are defined as
# i=0: films, i=1: tubes, i=2: spheroids 
# j=0: isolated pores, j=1: interconnected pores
cs=np.zeros([3,2])
# Add all porosity to phase 2,1 = connected spheroidal inclusions
cs[2,1]=1.0
nmax=100 # maximum number of iterations

# Material 1
alpha1=0.1    # aspect ratio of films
kappa=1.0     # shape parameter of tubes (0 = cuspate, 1=cylindrical)
alpha3=1.0    # aspect ratio of spheroids (1.0 =  spheres)
s_Ku1,s_Gu1,s_Kr1,s_Gr1=schmeling85.mod_s(Km,Gm,Ki,alpha1,kappa,alpha3,phi,nmax,cs)
s_Vpu1=np.sqrt((s_Ku1+4.0/3.0*s_Gu1)/roe)
s_Vpr1=np.sqrt((s_Kr1+4.0/3.0*s_Gr1)/roe)
s_Vsu1=np.sqrt(s_Gu1/roe)
s_Vsr1=np.sqrt(s_Gr1/roe)

# Material 2
alpha3=0.1
s_Ku2,s_Gu2,s_Kr2,s_Gr2=schmeling85.mod_s(Km,Gm,Ki,alpha1,kappa,alpha3,phi,nmax,cs)
s_Vpu2=np.sqrt((s_Ku2+4.0/3.0*s_Gu2)/roe)
s_Vpr2=np.sqrt((s_Kr2+4.0/3.0*s_Gr2)/roe)
s_Vsu2=np.sqrt(s_Gu2/roe)
s_Vsr2=np.sqrt(s_Gr2/roe)

# Calculate attenuation
b_qp1 = sca.qimax(b_Ku1,b_Kr1)
b_qp2 = sca.qimax(b_Ku2,b_Kr2)
b_qs1 = sca.qimax(b_Gu1,b_Gr1)
b_qs2 = sca.qimax(b_Gu2,b_Gr2)


# ==== Plot =============
# x limits (porosity)
pmin=0
pmax=1.0

fig = plt.figure()
ax = fig.add_subplot(221)
plt.xlabel('Phi')
plt.ylabel('M')
plt.xlim(pmin,pmax)
plt.ylim(0,90)
plt.plot(phi,Khs1,lw='2',c='r')
plt.plot(phi,Khs2,lw='2',c='r')

plt.plot(phi,b_Ku1,c='black',linewidth=1.5)
#plt.plot(phi,b_Kr1,c='black',linestyle='--',linewidth=2)
plt.plot(phi,b_Ku2,c='blue',linewidth=1.5)
#plt.plot(phi,b_Kr2,c='blue',linestyle='--',linewidth=2)

#plt.plot(phi,s_Ku2,c='darkgray',linewidth=1)
#plt.plot(phi,s_Kr2,c='darkgray',linestyle='--',linewidth=1)
#plt.plot(phi,s_Ku2,c='lightblue')
#plt.plot(phi,s_Kr2,c='lightblue')

plt.plot(phi,Ghs1,lw='2',c='r')
plt.plot(phi,Ghs2,lw='2',c='r')

plt.plot(phi,b_Gu1,c='black',linewidth=1.5)
#plt.plot(phi,b_Gr1,c='black',linestyle='--',linewidth=2)
plt.plot(phi,b_Gu2,c='blue',linewidth=1.5)
#plt.plot(phi,b_Gr2,c='blue',linestyle='--',linewidth=2)

#plt.plot(phi,s_Gu1,c='darkgray',linewidth=1)
#plt.plot(phi,s_Gr1,c='darkgray',linestyle='--',linewidth=1)
#plt.plot(phi,s_Gu2,c='lightblue')
#plt.plot(phi,s_Gr2,c='lightblue')

#-------------------------
ax = fig.add_subplot(222)
plt.xlabel('Phi')
plt.ylabel('V (km/s)')
plt.xlim(pmin,pmax)
plt.ylim(0.0,7.0)

plt.plot(phi,Vphs1,lw='1',c='r')
plt.plot(phi,Vphs2,lw='1',c='r')
plt.plot(phi,Vshs1,lw='1',c='r')
plt.plot(phi,Vshs2,lw='1',c='r')

plt.plot(phi,b_Vpu1,c='black',linewidth=1.5)
#plt.plot(phi,b_Vpr1,c='black',linestyle='--',linewidth=2)
plt.plot(phi,b_Vpu2,c='blue',linewidth=1.5)
#plt.plot(phi,b_Vpr2,c='blue',linestyle='--',linewidth=2)
#plt.plot(phi,b_Vpu3,c='brown',linewidth=1.5)

plt.plot(phi,b_Vsu1,c='black',linewidth=1.5)
#plt.plot(phi,b_Vsr1,c='black',linestyle='--',linewidth=2)
plt.plot(phi,b_Vsu2,c='blue',linewidth=1.5)
#plt.plot(phi,b_Vsr2,c='blue',linestyle='--',linewidth=2)
#plt.plot(phi,b_Vsu3,c='brown',linewidth=1.5)

#-------------------------
ax = fig.add_subplot(223)
plt.xlabel('Phi')
plt.ylabel('Q-1')
plt.xlim(pmin,pmax)
plt.ylim(0.0001,0.2)
plt.semilogy()
plt.plot(phi,b_qp1,lw=2,c='black')
plt.plot(phi,b_qp2,lw=2,c='blue')
plt.plot(phi,b_qs1,lw=1,c='black')
plt.plot(phi,b_qs2,lw=1,c='blue')

#-------------------------
ax = fig.add_subplot(224)
plt.xlabel('Phi')
plt.ylabel('Vp/Vs')
plt.xlim(pmin,pmax)
plt.ylim(1.5,3.5)
plt.plot(phi,b_Vpu1/b_Vsu1,c='black')
plt.plot(phi,b_Vpr1/b_Vsr1,c='black',linestyle='--',linewidth=2)
plt.plot(phi,b_Vpu2/b_Vsu2,c='blue')
plt.plot(phi,b_Vpr2/b_Vsr2,c='blue',linestyle='--',linewidth=2)

# Write plot to file
fileplot='sca_bvss_example.png'
plt.tight_layout()
plt.savefig(fileplot,dpi=300)

plt.show()



