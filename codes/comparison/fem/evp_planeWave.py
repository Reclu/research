#!/usr/bin/python

import numpy as np
from scipy.linalg import solve
from matplotlib import pyplot as plt
import os
import sys
from scipy import optimize
sys.path.append(os.path.expandvars('fem'))
from matrix1D import *
from linear_form import *
from grad1D import *

def equations(z,Selas,EPn,lam,mu,Sigy,H,eta,n,dt):
    S,EP = z
    Seq = computeS(S,EP,lam,mu,H)
    f = criterion(Seq,Sigy)
    return ((S-Selas+2*mu*(EP-EPn)),(EP-EPn-dt*creep_law(f,Seq,eta,n)))
    
def positive_part(x):
    return 0.5*(x+np.abs(x))

def creep_law(f,Seq,eta,n):
    return ((positive_part(f)/eta)**n)*np.sign(Seq)

def criterion(Seq,Sigy):
    return (np.abs(Seq)-Sigy)

def computeS(S,EP,lam,mu,H):
    KK = 3.0*(H/2.0) +(mu*(3.0*lam+2.0*mu))/(lam+2.0*mu)
    #KK = H +0.5*(4.*mu*(lam+mu)-lam**2)/(lam+2.0*mu)
    return (((2.0*mu)/(lam+2.0*mu))*S-KK*EP)

def bilinear(x,u_n,u,Sn,EPn,Pn,lam,mu,Sigy,H,eta,n,dt):
    #initialization
    h = x[1:len(x)]-x[:(len(x)-1)]
    eps_n = (u_n[1:len(u_n)]-u_n[:(len(u_n)-1)])/h
    eps = (u[1:len(u)]-u[:(len(u)-1)])/h
    S = np.zeros(len(eps))
    EP = np.zeros(len(eps))
    P  = np.zeros(len(eps))
    #Loop on integration points
    for i,DEFO in enumerate(eps):
        #(i) Elastic prediction
        Sx_trial = (lam+2.0*mu)*DEFO-2.0*mu*EPn[i]
        Sr_trial = lam*DEFO+mu*EPn[i]
        SS = Sx_trial-Sr_trial-3.0*(H/2.0)*EPn[i]
        # Sr_trial = lam*DEFO+(mu-0.5*lam)*EPn[i]
        # SS = Sx_trial-Sr_trial-H*EPn[i]
        Seq_trial = np.abs(SS)
        #(ii) Compute the criterion 
        f = Seq_trial - Sigy
        if (f<=0):
            #elastic step
            S[i] = Sx_trial
            EP[i] = EPn[i]
            P[i] = Pn[i]
        elif (f>0):
            #viscoplastic flow
            Sguess = (Sx_trial+Sn[i])/2.0
            result = optimize.root(equations,(Sguess,EPn[i]+1.0e-7),args=(Sx_trial,EPn[i],lam,mu,Sigy,H,eta,n,dt),method='hybr')
            S[i] = result.x[0] ; EP[i] = result.x[1]
            #Post-processing
            Seq = computeS(S[i],EP[i],lam,mu,H)
            f= criterion(Seq,Sigy)
            P[i]=Pn[i]+dt*creep_law(f,Seq,eta,n)*np.sign(Seq)
    return S,P,EP
    
"""
Implementation of the central differences method to solve the plane wave elastic-viscoplastic set of equations in dynamics with a linear isotropic strain hardening
"""
E=Young

#Define 1D mesh
#material parameters
lam = (nu*E)/(((1.0+nu)*(1.0-2.0*nu)))
mu = E/(2.0*(1.0+nu))
c = np.sqrt((lam+2.0*mu)/rho)       #Sound speed

#mesh
DX=length/Nelem

x = np.linspace(0.,length,Nelem+1)
x = np.linspace(-DX/2.,length+DX/2.,Nelem+1)

dx = x[1]-x[0]
t = np.zeros(NTmaxi)
y=x+(dx/2.0)
#Integration points location
y=y[:(len(y)-1)]
centroids=y
#Mass matrix
M = rho*Matrix1D(x,0,0,[],[])
#Mass matrix lumped (mass matrix approximated)
M_lumped = np.diag(np.sum(M,axis=1))

#definition and allocation of arrays
r = np.zeros((len(x),NTmaxi))
fext = np.zeros(len(x))
fint = np.zeros(len(x))
u = np.zeros((len(x),NTmaxi))
v = np.zeros((len(x),NTmaxi))
v_half = np.zeros(len(x))
a = np.zeros((len(x),NTmaxi))
epsilon =  np.zeros((len(x)-1,NTmaxi))
sigma =  np.zeros((len(x)-1,NTmaxi))
p =  np.zeros((len(x)-1,NTmaxi))
epsp =  np.zeros((len(x)-1,NTmaxi))

#Initial conditions
u[:,0] = np.zeros(len(x))
v[:,0] = np.zeros(len(x))

#Initial conditions

v_g = np.zeros(len(y))
for i in range(len(y)):
    if (i<len(range(len(y)))/2):
        v_g[i] = v0
    else:
        v_g[i] = -v0
#Definition at nodes from gauss points
wg = 2.0
v[:,0] = np.concatenate((np.array([v_g[0]]),(wg*v_g[1:len(v_g)]+wg*v_g[:(len(v_g)-1)])/(2.0*wg),np.array([v_g[-1]])),axis=0)

epsilon[:,0] = Grad1D(x,u[:,0])

#External forces
HEL = ((lam+2.0*mu)/(2.0*mu))*Sigy
applied_stress = -sigd
Load = np.zeros(NTmaxi)
Load[:] = 1.0*applied_stress
fext[0] = Load[0]
#Internal forces
sigma[:,0],p[:,0],epsp[:,0] = bilinear(x,np.zeros(len(x)),u[:,0],sigma[:,0],epsp[:,0],p[:,0],lam,mu,Sigy,H,eta,n,dx/c)
fint = Linear_form(x,1,[],sigma[:,0])
#Initial acceleration
a[:,0] = solve(M_lumped,(fext-fint)) 

v_half = v[:,0]+0.5*(dx/c)*a[:,0]

for i in range(NTmaxi)[1:]:
    #Current time step
    dt = dx/c
    if ((t[i-1]+dt)>timeOut):
	dt = timeOut - t[i-1]
    t[i]=t[i-1]+dt
    t_half = 0.5*(t[i]+t[i-1])
    #Displacement update
    u[:,i] = u[:,i-1]+dt*v_half
    epsilon[:,i] = Grad1D(x,u[:,i])
    #Compute internal forces
    sigma[:,i],p[:,i],epsp[:,i] = bilinear(x,u[:,i-1],u[:,i],sigma[:,i-1],epsp[:,i-1],p[:,i-1],lam,mu,Sigy,H,eta,n,dt)  #compute stress
    fint = Linear_form(x,1,[],sigma[:,i])  #compute predicted internal forces
    fext[0] = Load[i]                        #Loading
    #Solve for acceleration
    a[:,i] = solve(M_lumped,(fext-fint))
    v_half+=dt*a[:,i]

