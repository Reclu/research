#!/usr/bin/python

import numpy as np
import os
from scipy.linalg import solve
import pdb
from matrix1D import *
from linear_form import *
from grad1D import *
from bilinear import *

def loadingFunction(t,tfinal,tf):
    return np.sin(period*np.pi*t/tfinal)*(t<=tf)

def smoothSolution(x,c,t,tfinal,tf):
    if c*t < x:
        val=0.
    else:
        val=loadingFunction((t-x/c),tfinal,tf)
    return val
"""
Implementation of the central differences method to solve the 1D elastic-plastic set of equations
in dynamics with a linear isotropic strain hardening
"""
length=1.
Nelem=Mp             # Number of elements
rho=7800.
E=2.e11
Sigy=400.0e6           
c=np.sqrt(E/rho)
H = 10e9
HT=E*H/(E+H)
cp=np.sqrt(HT/rho)

A_0 = 1.

#mesh
NTmaxi=2*Nelem
x = np.linspace(0.,length,Nelem+1)
dx = x[1]-x[0]
t = np.zeros(NTmaxi)
CFL=1.
timeOut=0.5*length/c
timeUnload=2.*timeOut
y=x+(dx/2.0)
#Integration points location
y=y[:(len(y)-1)]

#Mass matrix
M = rho*A_0*Matrix1D(x,0,0,[],[])
#Mass matrix lumped (mass matrix approximated)
M_lumped = np.diag(np.sum(M,axis=1))

#definition and allocation of arrays
r = np.zeros((len(x),NTmaxi))
fext = np.zeros(len(x))
fint = np.zeros(len(x))
u = np.zeros((len(x),NTmaxi))
v = np.zeros((len(x),NTmaxi))
vth = np.zeros((len(x),NTmaxi))
v_half = np.zeros(len(x))
a = np.zeros((len(x),NTmaxi))
epsilon =  np.zeros((len(x)-1,NTmaxi))
sigma =  np.zeros((len(x)-1,NTmaxi))
anal =  np.zeros((len(x)-1,NTmaxi))
p =  np.zeros((len(x)-1,NTmaxi))
epsp =  np.zeros((len(x)-1,NTmaxi))
tangent_modulus = np.zeros(len(x)-1)

#Initial conditions
u[:,0] = np.zeros(len(x))
v[:,0] = np.zeros(len(x))
epsilon[:,0] = Grad1D(x,u[:,0])

#External forces
applied_stress = -1.e-4*Sigy
#Load = np.zeros(NTmaxi)
#Load[:Nt] = -1.0*applied_stress*A_0
#Load[Nt] = -1.0*applied_stress*A_0/2.0
fext[0] = -1.0*applied_stress*loadingFunction(0.,timeOut,timeUnload)*A_0
#Internal forces
sigma[:,0],p[:,0],epsp[:,0],tangent_modulus = bilinear(x,np.zeros(len(x)),u[:,0],epsp[:,0],p[:,0],E,Sigy,H)
fint = Linear_form(x,1,[],sigma[:,0]*A_0)
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
    #print 'time step '+str(i)+' t = '+str(t[i])
    #Displacement update
    u[:,i] = u[:,i-1]+dt*v_half
    epsilon[:,i] = Grad1D(x,u[:,i])
    #Compute internal forces
    sigma[:,i],p[:,i],epsp[:,i],tangent_modulus = bilinear(x,u[:,i-1],u[:,i],epsp[:,i-1],p[:,i-1],E,Sigy,H)  #compute stress
    fint = Linear_form(x,1,[],sigma[:,i]*A_0)  #compute predicted internal forces
    if (t[i]<timeUnload):
        fext[0] = -1.0*applied_stress*loadingFunction(t[i],timeOut,timeUnload)*A_0      #Loading
    elif (t[i]==timeUnload):
        fext[0] = -.5*applied_stress*loadingFunction(t[i],timeOut,timeUnload)*A_0      #Loading
    
    fext[0] = -1.0*applied_stress*loadingFunction(t[i],timeOut,timeUnload)*A_0   
    #Solve for acceleration
    a[:,i] = solve(M_lumped,(fext-fint))
    v_half+=dt*a[:,i]
    v[:,i]=v_half

    for k in range(len(x)):
        vth[k,i]=-applied_stress*smoothSolution(x[k],c,t[i],timeOut,timeUnload)/(rho*c)
    

    
    if (t[i]==timeOut):
        break
Increments=i

"""
#Contour plot of the stress
fig = plt.figure()
X,Y = np.meshgrid(y,t)
isoS = np.linspace(np.min(sigma[:len(y),:]), np.max(sigma[:len(y),:]),40, endpoint=True)
CP = plt.contourf(X, Y,sigma[:len(y),:].T,isoS)
plt.xlabel('x')
plt.ylabel('t')
plt.colorbar(CP)
plt.title('Stress (Pa)')
plt.show()
"""


