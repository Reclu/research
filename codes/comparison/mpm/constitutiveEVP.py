# !/usr/bin/python

import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import root

## Isotropic hardening

def isoHardening(z,Selas,Pn,E,Sigy,H,eta,n,dt):
    S,P = z
    f = np.abs(S)-H*P-Sigy
    rhs=(0.5*(f+np.abs(f))/eta)**n
    return ((S-Selas+E*(P-Pn)*np.sign(S)),(P-Pn-dt*rhs))
  

def stress_update_iso(eps,Sn,Epn,Pn,E,Sigy,H,eta,n,dt):
    #initialization
    S = np.zeros(len(eps))
    EP = np.zeros(len(eps))
    P  = np.zeros(len(eps))
    #Loop on integration points
    for i,DEFO in enumerate(eps):
        #(i) Elastic prediction
        Selas = E*(DEFO-EPn[i])
        #(ii) Compute the criterion 
        f = np.abs(Selas)-H*Pn[i]-Sigy
        if (f<=0):
            #elastic step
            S[i] = Selas
            EP[i] = EPn[i]
            P[i] = Pn[i]
        else:
            #viscoplastic correction
            Sguess = (Selas+Sn[i])/2.0
            #pdb.set_trace()
            result = root(isoHardening,(Sguess,Pn[i]+1.0e-7),args=(Selas,Pn[i],E,Sigy,H,eta,n,dt),method='hybr')
            S[i] = result.x[0] ; P[i] = result.x[1]
            EP[i] = EPn[i] + (P[i]-Pn[i])*np.sign(S[i])
    return S,EP,P




######################

## Kinematic hardening

def kinHardening(z,Selas,EPn,E,Sigy,H,eta,n,dt):
    S,EP = z
    f = (np.abs(S-H*EP)-Sigy)
    rhs= ((0.5*(f+np.abs(f))/eta)**n)*np.sign(S-H*EP)
    return ((S-Selas+E*(EP-EPn)),(EP-EPn-dt*rhs))
    
def stress_update_kin(eps,Sn,EPn,E,Sigy,H,eta,n,dt):
    #initialization
    S = np.zeros(len(eps))
    EP = np.zeros(len(eps))
    P  = np.zeros(len(eps))
    #Loop on integration points
    for i,DEFO in enumerate(eps):
        #(i) Elastic prediction
        Selas = E*(DEFO-EPn[i])
        #(ii) Compute the criterion 
        f = np.abs(Selas-H*EPn[i])-Sigy
        if (f<=0):
            #elastic step
            S[i] = Selas
            EP[i] = EPn[i]
        else:
            #viscoplastic correction
            Sguess = (Selas+Sn[i])/2.0
            result = scipy.optimize.root(kinHardening,(Sguess,EPn[i]+1.0e-7),args=(Selas,EPn[i],E,Sigy,H,eta,n,dt),method='hybr')
            S[i] = result.x[0] ; EP[i] = result.x[1]
    return S,EP
