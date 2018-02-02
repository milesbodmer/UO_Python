#!/opt/local/bin/python2.7
# This script calculates the effective elastic properties of 
# multi-phase media using functions from sca.py and schemeling85.py
# Plotting is done with Matplotlib
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
import critical

# Import common libraries
import numpy as np
import matplotlib.pyplot as plt
from math import pi,e,log,log10,sqrt

###############################################################

fluid='basalt'
# Elastic moduli are in GPa, density is in g/cm^2
if fluid == 'magma':
    vpvs = 1.73
    Ki = 1.27e1
    Gi = 1.e-12
    rof = 2.3
elif fluid == 'water':
    vpvs = 1.8
    Ki = 2.2
    Gi = 1.e-12
    rof = 1.05
if fluid == 'basalt':
    Ki = 1.491e1
    Gi = 1.e-12
    rof = 2.7
    
# Bulk and shear modulus of solid matrix
if fluid in ['magma','water']:
  vp0=6.5
  vs0=vp0/vpvs
  rom=2.7
  Gm = rom*vs0**2
  Km = rom*(vp0**2-4.0/3.0*vs0**2)
elif fluid == 'basalt':
## Use values from Mainprice
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

# Calculate density of the composite
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



# Use Berryman's method
#--------------------------------
# Material 1: distribution of aspect ratios
# aspect ratio of solid grains
am=0.1
# Size of array of aspect ratios of fluid inclusions
na=21
# Span of aspect ratios (orders of magnitude either side of central value)
# This is 2*sigma of the gaussian distribution
orders=1.0
# Central value of aspect ratio distribution (power of 10)
mu=-1.0
# setup arrays
alog=np.logspace(mu-orders,mu+orders,na)
alin=np.linspace(mu-orders,mu+orders,na)
aones=np.ones(na)
print alog

# Distribution of porosity over aspect ratios:
# Constant:
#ca = aa*0.+1./len(aa)
# Truncated normal
ca = np.exp( -(alin-mu)**2. / (2.*(orders/2.0)**2.) )
ca = ca/np.sum(ca,axis=0)
print ca

# Insert solid phase aspect ratio at start of array
aa=np.insert(alog,0,am)

# Generate K and G arrays of same size
KK=np.ones(na)*Ki
KK=np.insert(KK,0,Km)
GG=np.ones(na)*Gi
GG=np.insert(GG,0,Gm)

# Generate porosity matrix
cam=np.ones([nc,na])*ca
cam=cam.transpose()
pp=np.ones([na,nc])*phi*cam
phim=1.-np.sum(pp,axis=0)
pp=np.insert(pp,0,phim,axis=0)

# Unrelaxed, high frequency case
b_Ku1,b_Gu1 = sca.mod_multi(KK,GG,aa,pp)
b_Vpu1=np.sqrt((b_Ku1+4.0/3.0*b_Gu1)/roe)
b_Vsu1=np.sqrt(b_Gu1/roe)

# Relaxed, low frequency case - Gassman
KK=np.ones(na)*Ki
KK=np.insert(KK,0,1e-12)
K,G = sca.mod_multi(KK,GG,aa,pp)
b_Kr1,b_Gr1 = sca.gassman_2(K,G,Km,Gm,Ki,Gi,phi)
b_Vpr1=np.sqrt((b_Kr1+4.0/3.0*b_Gr1)/roe)
b_Vsr1=np.sqrt(b_Gr1/roe)


#--------------------------------
# Material 2: single aspect ratio
am=0.1
ai=0.1
# Unrelaxed, high frequency
b_Ku2,b_Gu2 = sca.mod_b(Km,Ki,Gm,Gi,am,ai,1.-phi,phi)
b_Vpu2=np.sqrt((b_Ku2+4.0/3.0*b_Gu2)/roe)
b_Vsu2=np.sqrt(b_Gu2/roe)

# Relaxed, low frequency - Gassman
K,G = sca.mod_b(Km,0.0,Gm,0.0,am,ai,1.-phi,phi)
b_Kr2,b_Gr2 = sca.gassman_2(K,G,Km,Gm,Ki,Gi,phi)
b_Vpr2=np.sqrt((b_Kr2+4.0/3.0*b_Gr2)/roe)
b_Vsr2=np.sqrt(b_Gr2/roe)



# Calculate attenuation
b_qp1 = sca.qimax(b_Ku1,b_Kr1)
b_qp2 = sca.qimax(b_Ku2,b_Kr2)
b_qs1 = sca.qimax(b_Gu1,b_Gr1)
b_qs2 = sca.qimax(b_Gu2,b_Gr2)


# ==== Plot =============
pmin=0
pmax=1

fig = plt.figure()
ax = fig.add_subplot(221)
plt.xlabel('Phi')
plt.ylabel('K')
plt.xlim(pmin,pmax)
plt.ylim(0,90)
plt.plot(phi,Khs1,lw=2,c='r')
plt.plot(phi,Khs2,lw=2,c='r')

plt.plot(phi,b_Ku1,c='black',linewidth=1.5)
plt.plot(phi,b_Kr1,c='black',linestyle='--',linewidth=1.5)
plt.plot(phi,b_Ku2,c='blue',linewidth=1.5)
plt.plot(phi,b_Kr2,c='blue',linestyle='--',linewidth=1.5)

plt.plot(phi,Ghs1,lw=2,c='r')
plt.plot(phi,Ghs2,lw=2,c='r')

plt.plot(phi,b_Gu1,c='black',linewidth=1.5)
plt.plot(phi,b_Gr1,c='black',linestyle='--',linewidth=1.5)
plt.plot(phi,b_Gu2,c='blue',linewidth=1.5)
plt.plot(phi,b_Gr2,c='blue',linestyle='--',linewidth=1.5)

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
plt.plot(phi,b_Vpr1,c='black',linestyle='--',linewidth=1.5)
plt.plot(phi,b_Vpu2,c='blue',linewidth=1.5)
plt.plot(phi,b_Vpr2,c='blue',linestyle='--',linewidth=1.5)

plt.plot(phi,b_Vsu1,c='black',linewidth=1.5)
plt.plot(phi,b_Vsr1,c='black',linestyle='--',linewidth=1.5)
plt.plot(phi,b_Vsu2,c='blue',linewidth=1.5)
plt.plot(phi,b_Vsr2,c='blue',linestyle='--',linewidth=1.5)

#-------------------------
ax = fig.add_subplot(223)
plt.xlabel('Phi')
plt.ylabel('Q-1')
plt.xlim(pmin,pmax)
plt.ylim(0.0001,1.0)
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
plt.ylim(1.5,2.5)

plt.plot(phi,b_Vpu1/b_Vsu1,c='black',linewidth=1.5)
plt.plot(phi,b_Vpu2/b_Vsu2,c='blue',linewidth=1.5)
plt.plot(phi,b_Vpr1/b_Vsr1,c='black',linestyle='--',linewidth=1.5)
plt.plot(phi,b_Vpr2/b_Vsr2,c='blue',linestyle='--',linewidth=1.5)


plt.tight_layout()
fileplot='sca_distribution.png'
plt.savefig(fileplot,dpi=200)

plt.show()


