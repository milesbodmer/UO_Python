# Self consistent approximation based on Schmeling 1985
# Copyright (c) 2012, Jan Philipp Kruse and Harro Schmeling - Universitat Frankfurt
# All rights reserved
# Translated from Matlab to Python by Michele Paulatto - Imperial College London
#
# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met:
# 
# Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
# OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


import numpy as np
from math import pi, acos, sqrt, isinf

def mod_s(K0,mu0,Kf,alpha1,kappa,alpha3,beta,nmax,c):
# K0 - bulk modulus matrix (default is 0.66 [M Pa])
# mu0 - shear modulus matrix (default is 0.4 [MPa])
# Kf - bulk modulus melt (default is 0.2 [MPa])
# alpha1 - aspect ratios of films (default is 0.01)
# kappa - tube shape, 0 for tapered, large for cylinders
# alpha3 - aspect ratios of spheroids - (default is 1)
# beta - range of melt fractions (array)
# nmax - iterations, try 10 - 100 (default is 1000)
# c_all - concentrations of different geometries, fraction of all melt, sum must be = 1
# Returns Ku, Kr, muu, mur

  betanum     = len(beta)
  missmax1    = 10**-11*K0
  missmax2    = 10**-11*mu0
									
# Set up concentrations
#c_all   = numpy.zeros(3,2)
#c_all[2,1] = 1.0
#    c_all[0,0]  c1i - fraction of films of all melt
#    c_all[0,1]  c1c - fraction of films of all melt 
#    c_all[1,0]  c2i - fraction of tubes of all melt 
#    c_all[1,1]  c2c - fraction of tubes of all melt 
#    c_all[2,0]  c3i - fraction of spheroids of all melt 
#    c_all[2,1]  c3c - fraction of spheroids of all melt 
                                
  Ku     = np.zeros_like(beta)
  muu    = np.zeros_like(beta)
  Au     = np.zeros_like(beta)
  nuu    = np.zeros_like(beta)
  Kbar   = np.zeros_like(beta)
  Abar   = np.zeros_like(beta)
  nubar  = np.zeros_like(beta)
  Kiso   = np.zeros_like(beta)
  muisoo = np.zeros_like(beta)
  Aiso   = np.zeros_like(beta)
  nuiso  = np.zeros_like(beta)
  Kr     = np.zeros_like(beta)
  mur    = np.zeros_like(beta)
  nur    = np.zeros_like(beta)
  phi_tmp   = phi(alpha3)
  g_tmp     = g(alpha3,phi_tmp)

# Loop over porosity array
  for i in range(len(beta)): 
# UNRELAXED MODULI #
# START ITERATION INDEX #
      n                   = 1
# SET INITIAL CONDITIONS #      
# K #
      Ku_tmp      = [0.0]*(nmax+1)
      Ku_tmp[0]   = K0
      Ku_tmp[1]   = 0.0
# MU #
      muu_tmp     = [0.0]*(nmax+1)
      muu_tmp[0]  = mu0      
# INITIAL ERROR #
      miss1       = K0
      miss2       = mu0
# ---------------------- #
# RECURSION FOR K_U & MU_U #
# ------------------------ # 
      while ( miss1 >= missmax1 or miss2 >= missmax2 ) and n <= nmax:         
          K_tmp           = 0.0
          mu_tmp          = 0.0
          for j in range(3):
              nu_tmp      = poisson(Ku_tmp[n-1],muu_tmp[n-1])
              R_tmp       = R(Ku_tmp[n-1],muu_tmp[n-1])
              B_tmp       = B(j,Kf,Ku_tmp[n-1],nu_tmp,kappa,muu_tmp[n-1])
              theta_tmp   = theta(j,Ku_tmp[n-1],muu_tmp[n-1],nu_tmp,alpha1,kappa,g_tmp,phi_tmp,R_tmp)
              D_tmp       = D(j,Kf,Ku_tmp[0],theta_tmp)
              A_tmp       = A(j,alpha1,alpha3,Ku_tmp[n-1],g_tmp,R_tmp,D_tmp,phi_tmp,B_tmp,muu_tmp[n-1],nu_tmp,kappa,Kf)
              mu_tmp      = mu_tmp+muunrelaxed(A_tmp,c[j,0]*beta[i],c[j,1]*beta[i])
              K_tmp       = K_tmp+KModulunrelaxed(Kf,K0,c[j,0]*beta[i],c[j,1]*beta[i],theta_tmp,Ku_tmp[n-1])
         
          Ku_tmp[n]       = 1./(1./K0+K_tmp)
          muu_tmp[n]      = 1./(1./mu0+mu_tmp)          
# ERROR FOR MU AND K #
          miss1           = abs(Ku_tmp[n]-Ku_tmp[n-1])
          miss2           = abs(muu_tmp[n]-muu_tmp[n-1])         
          n               = n+1            
# SET VALUES FOR UNRELAXED TERMS #
# ------------------------------ #
      Ku[i]               = Ku_tmp[n-1]
      muu[i]              = muu_tmp[n-1]
      Au[i]               = A_tmp
      nuu[i]              = nu_tmp
# ------------------------------ #
      
# KMODULI' and RELAXED SHEARMODULI #
# -------------------------------- #
      
# START INDEX #
      n                   = 1     
# SET INITIAL CONDITIONS #
# ---------------------- #
# K'                               #
      Kbar_tmp     = [0]*(nmax+1)
      Kbar_tmp[0]  = K0
      Kbar_tmp[1]  = 0
# MY_R                      #
      mubar_tmp    = [0]*(nmax+1)
      mubar_tmp[0] = mu0
# INITIAL ERROR             #
      miss1        = K0
      miss2        = mu0
# -------------------------------- #
# RECURSION FOR K' & MY_R          #
# -------------------------------- #
      while ( miss1 >= missmax1 or miss2>=missmax2 ) and n <= nmax:
          K_tmp           = 0.0
          mu_tmp          = 0.0
          for j in range(3):
              nu_tmp      = poisson(Kbar_tmp[n-1],mubar_tmp[n-1])
              R_tmp       = R(Kbar_tmp[n-1],mubar_tmp[n-1])
              B1_tmp      = B(j,Kf,Kbar_tmp[n-1],nu_tmp,kappa,mubar_tmp[n-1])
              B2_tmp      = B(j,0.0,Kbar_tmp[n-1],nu_tmp,kappa,mubar_tmp[n-1])
              theta_tmp   = theta(j,Kbar_tmp[n-1],mubar_tmp[n-1],nu_tmp,alpha1,kappa,g_tmp,phi_tmp,R_tmp)
              D1_tmp      = D(j,Kf,Kbar_tmp[0],theta_tmp)
              D2_tmp      = D(j,0.0,Kbar_tmp[0],theta_tmp)
              A1_tmp      = A(j,alpha1,alpha3,Kbar_tmp[n-1],g_tmp,R_tmp,D1_tmp,phi_tmp,B1_tmp,mubar_tmp[n-1],nu_tmp,kappa,Kf)
              A2_tmp      = A(j,alpha1,alpha3,Kbar_tmp[n-1],g_tmp,R_tmp,D2_tmp,phi_tmp,B2_tmp,mubar_tmp[n-1],nu_tmp,kappa,0.0)
              mu_tmp      = mu_tmp+murelaxed(A1_tmp,A2_tmp,c[j,0]*beta[i],c[j,1]*beta[i])
              K_tmp       = K_tmp+KModulbar(Kf,K0,theta_tmp,Kbar_tmp[n-1],c[j,0]*beta[i],c[j,1]*beta[i])
          
          Kbar_tmp[n]     = 1./(1./K0+K_tmp)
          mubar_tmp[n]    = 1./(1./mu0+mu_tmp)
          
# ERROR FOR MY_R AND K'            #
          miss1           = abs(Kbar_tmp[n]-Kbar_tmp[n-1])
          miss2           = abs(mubar_tmp[n]-mubar_tmp[n-1])
# -------------------------------- #        
          n               = n+1
            
# SET VALUES FOR BAR AND MY_RELAXED TERMS #
# --------------------------------------- #
      mur[i]              = mubar_tmp[n-1]
      Kbar[i]             = Kbar_tmp[n-1]
      Abar[i]             = A_tmp
      nubar[i]            = nu_tmp
# --------------------------------------- # 
      
# --------------------------------------- #
# KMODULI ISOLATED                        #
# --------------------------------------- #
      
# START INDEX                             #
      n                   = 1     
# SET INITIAL CONDITIONS                  #
# --------------------------------------- #
# KISO                                    #
      Kiso_tmp     = [0.0]*(nmax+1)
      Kiso_tmp[0]  = K0
      Kiso_tmp[1]  = 0.0
# MUISO                            #
      muiso_tmp     = [0.0]*(nmax+1)
      muiso_tmp[0] = mu0
# INITIAL ERROR                    #
      miss1        = K0
      miss2        = mu0
      
# --------------------------------------- #
# RECURSION FOR KISO & MUISO              #
# --------------------------------------- #
      
      while ( miss1 >= missmax1 or miss2 >= missmax2 ) and n <= nmax:
          K_tmp           = 0.0
          mu_tmp          = 0.0
          
          for j in range(3):
              betamod     = (c[j,0]*beta[i])/(1-c[j,1]*beta[i])
              nu_tmp      = poisson(Kiso_tmp[n-1],muiso_tmp[n-1])
              R_tmp       = R(Kiso_tmp[n-1],muiso_tmp[n-1])
              B_tmp       = B(j,Kf,Kiso_tmp[n-1],nu_tmp,kappa,muiso_tmp[n-1])
              theta_tmp   = theta(j,Kiso_tmp[n-1],muiso_tmp[n-1],nu_tmp,alpha1,kappa,g_tmp,phi_tmp,R_tmp)
              D_tmp       = D(j,Kf,Kiso_tmp[0],theta_tmp)
              A_tmp       = A(j,alpha1,alpha3,Kiso_tmp[n-1],g_tmp,R_tmp,D_tmp,phi_tmp,B_tmp,muiso_tmp[n-1],nu_tmp,kappa,Kf)
              mu_tmp      = mu_tmp+muiso(A_tmp,betamod)
              K_tmp       = K_tmp+KModuliso(Kf,K0,betamod,theta_tmp,Kiso_tmp[n-1])
          
          Kiso_tmp[n]     = 1./(1./K0+K_tmp)
          muiso_tmp[n]    = 1./(1./mu0+mu_tmp)
          
# ERROR FOR MY_ISO AND K_ISO              #
          miss1           = abs(Kiso_tmp[n]-Kiso_tmp[n-1])
          miss2           = abs(muiso_tmp[n]-muiso_tmp[n-1])
# --------------------------------------- #
          n               = n+1     
# SET VALUES FOR ISOLATED TERMS           #
# --------------------------------------- #
      Kiso[i]             = Kiso_tmp[n-1]
      muisoo[i]           = muiso_tmp[n-1]
      Aiso[i]             = A_tmp
      nuiso[i]            = nu_tmp
# --------------------------------------- #
    
# SET VALUES FOR RELAXED TERMS            #
# --------------------------------------- #
      Kr[i]               = KModulrelaxed(Kf,Kiso[i],Kbar[i],(c[0,1]+c[1,1]+c[2,1])*beta[i])
      nur[i]              = poisson(Kr[i],mur[i])
# --------------------------------------- #
    
  return Ku,muu,Kr,mur





#-------------------------------------------
def KModulbar(Kf,K0,theta,Kst,betai,betac):
  if Kf == 0 and theta != 0:
    Kbarr = (betai/theta**-1)+theta*betac
  elif theta == 0:
    Kbarr = 0.0
  else:
    Kbarr = (((1./Kf-1./K0)*betai)/(1+theta**-1*(1./Kf-1./Kst)))+theta*betac
  return Kbarr
  
#-------------------------------------------
def KModuliso(Kf,K0,betamod,theta,Ku):
  if Kf == 0 and theta != 0:
    Kiso = betamod/theta**-1
  else:
    Kiso = ((1./Kf-1./K0)*betamod)/(1+theta**-1*(1./Kf-1./Ku))
  return Kiso
  
#-------------------------------------------
def KModulrelaxed(Kf,Ku,Kst,betac):
  if betac == 0:
    F = 0.0
  else:
    F = ((Kf*(Ku-Kst))/(betac*(Ku-Kf)))
  Krel = Ku*((Kst+F)/(Ku+F))
  return Krel

#-------------------------------------------
def KModulunrelaxed(Kf,K0,betai,betac,theta,Ku):
  if Kf == 0:
    if betai == 0:
        t1 = 0.0
    else:
        t1 = betai*theta
    if betac == 0:
        t2 = 0.0
    else:
        t2 = betac*theta
    Kunr = t1+t2
  else:
    Kunr = (((1./Kf)-(1./K0))*(betai+betac))/(1.+theta**-1*((1./Kf)-(1./Ku)))
  return Kunr
  
#-------------------------------------------
def muiso(A,betamod):
  muis = A*betamod
  return muis

#-------------------------------------------
def murelaxed(A1,A2,betai,betac):
  muu = A1*betai+betac*A2
  return muu

#-------------------------------------------
def muunrelaxed(A,betai,betac):
  if betai == 0:
    A1 = 0.0
  else:
    A1 = A*betai
  if betac == 0:
    A2 = 0.0
  else:
    A2 = A*betac
  muu = A1+A2
  return muu

#====================================================================================
# Low level functions

#-------------------------------------------
def A(flag,alpha1,alpha3,K,g,R,D,phi,B,mu,nu,kappa,Kf):
  if flag == 0:
    A   = (8./(15*pi*mu))*((1-nu)/(2-nu))*((2-nu)*D+3)*(1./alpha1)
  elif flag == 1:
    if Kf > 0:
      s1  = ((2*(1-nu)*((2+kappa)**2+2))/((2+kappa)**2-2))-1+2*nu
      s2  = 2*(1-nu)*((2+kappa)**2+2)-(1-2*nu)*((2+kappa)**2-2)
      s3  = (-2*(1-nu)*((2+kappa)**2+2)+mu*((2+kappa)**2-2)*(1./K-1./Kf-((1-2*nu)**2./(2*mu*(1+nu)))))**-1
      A   = B + (1./(15*mu))*s1*s2*s3
    else:
      A   = 0.0
  elif flag == 2:
    if alpha3 == 1:
        A   = (5./mu)*(4*mu+3*K)/(8*mu+9*K)
    else:
        t1  = 2*(1-.5*(-(1+alpha3**2)*(g/alpha3**2)+R*(2-phi+(1+alpha3**2)*(g/alpha3**2))))**-1
        t2  = (1-(1./4)*(3*phi+g-R*(g-phi)))**-1
        t3  = (1-(1./4)*(3*phi+g-R*(g-phi)))
        t4  = (B*phi*(3-4*R)+g-R*(g+phi-4./3))
        t5  = 2*(R*(g+phi)-g+B*(1-phi)*(3-4*R))
        t6  = (1-(1./8)*(9*phi+3*g-R*(5*phi+3*g))+.5*B*phi*(3-4*R))
        t7  = 2*(B*(1-phi)*(3-4*R)-1+(3./2)*phi+.5*g-.5*R*(5*phi+g-4))
        t8  = (.5*B*phi*(3-4*R)+.5*g-.5*R*(g-phi))
        t9  = (1-(1./4)*(3*phi+g-R*(g-phi)))**-1
        t10 = (R*((3./2)*g+(5./2)*phi)-(3./2)*(g+phi)+B*(3-4*R)-.5*(3*B-1)*(3-4*R)*(g+phi-R*(g-phi+2*phi**2)))**-1
        A   = (1./(5*mu))*(t1+t2+(t3*t4+t5*t6-t7*t8)*t9*t10)
  return A

#-------------------------------------------
def B(flag,Kf,K,nu,kappa,mu):
  if flag == 0:
    B = 0.0
  elif flag == 1:
    B = ((2./15)*(1+nu)+(2./15)*(1-nu)* 
        (((2+kappa)**2+2)/((2+kappa)**2-2))+ 
        (8./5)*((3./2)-nu)* 
        ((2+kappa)**2./((2+kappa)**2-2)))*(1./mu)
  elif flag == 2:
    B = Kf/(3*K)
  return B
  
#-------------------------------------------
def D(flag,Kf,K0,theta1):
  if Kf == 0:
    D = 1
  else:
    if flag == 0:
#      print K0,Kf,theta1,theta1+1./Kf
      D = (1./Kf-1./K0)*(theta1+1./Kf)**-1
    else:
      D = 0.0
  return D
  
#-------------------------------------------
def g(alpha,phi):
  if alpha != 1:
    g = alpha**2*(1-alpha**2)**-1*(3*phi-2)
  else:
    g = 0.0
  return g

#-------------------------------------------
def phi(alpha3):
  if alpha3 != 1:
    phiu = alpha3*(1-alpha3**2)**(-3./2)*(acos(alpha3)-alpha3*sqrt(1-alpha3**2))
  else:
    phiu = 0.0
  return phiu

#-------------------------------------------
def poisson(Ku,mueu):
  if isinf(Ku):
    nuu = .5
  else:
    nuu = (3*Ku-2*mueu)/(6*Ku+2*mueu)
  return nuu

#-------------------------------------------
def R(Ku,mueu):
  R = (3*mueu)/(3*Ku+4*mueu)
  return R

#-------------------------------------------
def theta(flag,K,mu,nu,alpha1,kappa,g,phi,R):
  if flag == 0:
    if nu == .5:
      thetaout = (2*(1-nu))/(pi*alpha1*mu)
    else:
      thetaout = (4./(3*pi))*(1./K)*((1-nu**2)/(1-2*nu))*(1./alpha1)
  elif flag == 1:
    if nu == .5:
      thetaout = (4*(1-nu)*((2+kappa)**2+2))/(mu*((2+kappa)**2)-2)+(1-2*nu)/(3*K)
    else:
      thetaout = (2./(3*K))*(((2*(1-nu**2))/(1-2*nu))*(((2+kappa)**2+2)/((2+kappa)**2-2))+.5*(1-2*nu))
  elif flag == 2:
    if g == 0 and phi == 0:
        if nu == .5:
          thetaout = (9*(1-nu))/(4*mu*(1+nu))
        else:
          thetaout = (3./(2*K))*((1-nu)/(1-2*nu))
    else:
      d1 = (1-(3*(g+phi)/2-R*((3./2)*g+(5./2)*phi-4./3)))
      d2 = (1-(1+(3./2)*(g+phi)-R*((3./2)*g+(5./2)*phi))+.5*(3-4*R)*(g+phi-R*(g-phi+2*phi**2)))**-1
      thetaout = K**-1*d1*d2
  return thetaout

