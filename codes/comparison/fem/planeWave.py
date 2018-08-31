#!/usr/bin/python

import numpy as np
import os
from scipy.linalg import solve
from matplotlib import pyplot as plt
from matplotlib import animation
import pdb
import sys

sys.path.append(os.path.expandvars('fem'))
from matrix1D import *
from linear_form import *
from grad1D import *

def bilinear(x,u_n,u,EPn,Pn,lam,mu,Sigy,H):
    #initialization
    h = x[1:len(x)]-x[:(len(x)-1)]
    eps_n = (u_n[1:len(u_n)]-u_n[:(len(u_n)-1)])/h
    eps = (u[1:len(u)]-u[:(len(u)-1)])/h
    S = np.zeros(len(eps))
    EP = np.zeros(len(eps))
    P  = np.zeros(len(eps))
    TM  = np.zeros(len(eps))
    #Loop on integration points
    for i,DEFO in enumerate(eps):
        #(i) Elastic prediction
        Sx_trial = (lam+2.0*mu)*DEFO-2.0*mu*EPn[i]
        Sr_trial = lam*DEFO+mu*EPn[i]
        SS = Sx_trial-Sr_trial-3.0*(H/2.0)*EPn[i]
        Seq_trial = np.abs(SS)
        #(ii) Compute the criterion 
        f = Seq_trial - Sigy
        if (f<=0):
            #elastic step
            S[i] = Sx_trial
            EP[i] = EPn[i]
            P[i] = Pn[i]
            TM[i] = lam+2.0*mu
        elif (f>0):
            #elastoplastic step
            dP = f/(3.0*(mu+(H/2.0)))
            P[i] = Pn[i]+dP
            EP[i] = EPn[i]+(P[i]-Pn[i])*np.sign(SS)
            S[i] = Sx_trial - 2.0*mu*dP*np.sign(SS)
            TM[i] = (lam+2.0*mu)-8*(mu**2/(3.0*H+6.0*mu))
    return S,P,EP,TM

CFL=1.
"""
Implementation of the central differences method to solve the 1D elastic-plastic set of equations
in dynamics with a linear isotropic strain hardening
"""
E=Young
rho=rho
lam = (nu*E)/(((1.0+nu)*(1.0-2.0*nu)))
mu= E/(2.0*(1.0+nu))
HT = (lam+2.0*mu)-8*(mu**2/(3.0*H+6.0*mu))

c = np.sqrt((lam+2.0*mu)/rho)
cp=np.sqrt(HT/rho)

#Define 1D mesh
# Nelem*=3
# NTmaxi*=3
# length*=3
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
sig =  np.zeros((len(x)-1,NTmaxi))
p =  np.zeros((len(x)-1,NTmaxi))
epsp =  np.zeros((len(x)-1,NTmaxi))
tangent_modulus = np.zeros(len(x)-1)

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
fext[0] = -sigd
#Internal forces
sig[:,0],p[:,0],epsp[:,0],tangent_modulus = bilinear(x,np.zeros(len(x)),u[:,0],epsp[:,0],p[:,0],lam,mu,Sigy,H)
fint = Linear_form(x,1,[],sig[:,0])
#Initial acceleration
a[:,0] = solve(M_lumped,(fext-fint)) 

v_half = v[:,0]+0.5*(dx/c)*a[:,0]

for i in range(NTmaxi)[1:]:
    #Current time step
    dt = CFL*dx/c
    if ((t[i-1]+dt)>timeOut):
	dt = timeOut - t[i-1]
    t[i]=t[i-1]+dt
    t_half = 0.5*(t[i]+t[i-1])
    #print 'time step '+str(i)+' t = '+str(t[i])
    #Displacement update
    u[:,i] = u[:,i-1]+dt*v_half
    epsilon[:,i] = Grad1D(x,u[:,i])
    #Compute internal forces
    sig[:,i],p[:,i],epsp[:,i],tangent_modulus = bilinear(x,u[:,i-1],u[:,i],epsp[:,i-1],p[:,i-1],lam,mu,Sigy,H)  #compute stress
    fint = Linear_form(x,1,[],sig[:,i])  #compute predicted internal forces
    if (t[i]<timeUnload):
        fext[0] = -sigd      #Loading
    else:
        fext[0] = 0.0;fext[-1] = 0.0
    #Solve for acceleration
    a[:,i] = solve(M_lumped,(fext-fint))
    v_half+=dt*a[:,i]
    if (t[i]==timeOut):
        break


