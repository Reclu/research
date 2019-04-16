#!/usr/bin/python

import numpy as np
import os
from scipy.linalg import solve
from matplotlib import pyplot as plt
from matplotlib import animation
import pdb
import RHS

#Define 1D mesh
Nelem = Mp
Nnode=2*Nelem+2
#Intercells points
length = 1.0
x = np.linspace(0.,length,Nelem+1)
dx = x[1]-x[0]

#Define material parameters
rho=7800.
E=2.e11
Sy=400.0e6           
c=np.sqrt(E/rho)

t_order=1
#Define parameters
CFL = 0.1
timeOut = length/c
tf=0.5*timeOut
dt = CFL*dx/c
NTmaxi = int(timeOut/dt)
t = np.zeros(NTmaxi)

#Definition of arrays
#Ghost cells at both ends of the bar are used to prescribe boundary conditions
rho = np.zeros((Nnode,NTmaxi))
R = np.zeros(Nnode)

# Boundary conditions on left ghost node
r0=1.


R[:] = rho[:,0]

dim=1 # number of unknowns
computeFlux=RHS.Flux(dim,0.,0.,c,x)

#Definition of functions

def Jac(Xe):
    return ((Xe[1]-Xe[0])/2.0)

def Lin_shp(idx,xi,Xe):
    if (idx == 1):
        return np.array([(1.0-xi)/2.0,(1.0+xi)/2.0])
    elif (idx == 2):
        return (np.array([-np.ones(len(xi)),np.ones(len(xi))])/(Xe[1]-Xe[0]))

def Lin_eval(xi,Xe):
    return ((Xe[0]*(-xi+1.0)+Xe[1]*(xi+1.0))/2.0)

def computeMatrix(x,Nnodes):
    Nelem = (Nnodes-2)/2
    T10 = np.array([np.arange(1,Nnodes-1,2),np.arange(2,Nnodes,2)]).T
    xi_integ = np.array([0.774596669241483,0.0,-0.774596669241483])
    w_integ = np.array([0.555555555555556,0.888888888888889,0.555555555555556])
    M = np.zeros((Nnodes,Nnodes))
    for i in range(Nelem):
        mapp = T10[i,:]
        Xe = x[[i,i+1]]
        loc_jac = Jac(Xe)
        M1 = Lin_shp(1,xi_integ,Xe)
        M2 = Lin_shp(1,xi_integ,Xe)
        D = np.array(w_integ)
        M[mapp[0]:mapp[1]+1,mapp[0]:mapp[1]+1] = (loc_jac*np.dot(M1,np.dot(np.diag(D),M2.T))) 
    return M

def UpdateState(dt,x,R):
    Nnodes = np.shape(R)[0]
    #Compute mass matrix
    M=computeMatrix(x,Nnodes)
    Ml = np.sum(M,axis=1)
    f=computeFlux(R)
    for i in range(1,Nnodes-1):
        R[i]+=dt*(f[i]/Ml[i])
    

def UpdateStateRK2(dt,x,R):
    Nnodes = np.shape(R)[0]
    #Compute mass matrix
    M=computeMatrix(x,Nnodes)
    # Lumped mass matrix
    Ml = np.sum(M,axis=1)
    f=computeFlux(R)
    k1=np.zeros(Nnodes)
    k2=np.zeros(Nnodes)
    for i in range(1,Nnodes-1):
        k1[i]=dt*(f[i]/Ml[i])
    f=computeFlux(R+k1)
    for i in range(1,Nnodes-1):
        k2[i]=dt*(f[i]/Ml[i])
    R+=0.5*(k1+k2)

def UpdateStateRK4(dt,x,R):
    Nnodes = np.shape(R)[0]
    k1=np.zeros(Nnodes)
    k2=np.zeros(Nnodes)
    k3=np.zeros(Nnodes)
    k4=np.zeros(Nnodes)
    #Compute mass matrix
    M=computeMatrix(x,Nnodes)
    Ml = np.sum(M,axis=1)
    f=computeFlux(R)
    for i in range(1,Nnodes-1):
        k1[i]=dt*(f[i]/Ml[i])
    f=computeFlux(R+dt*k1/2.)
    for i in range(1,Nnodes-1):
        k2[i]=dt*(f[i]/Ml[i])
    f=computeFlux(R+dt*k2/2.)
    for i in range(1,Nnodes-1):
        k3[i]=dt*(f[i]/Ml[i])
    f=computeFlux(R+dt*k3)
    for i in range(1,Nnodes-1):
        k4[i]=dt*(f[i]/Ml[i])
    R+=(k1+2*k2+2*k3+k4)/6.
    

print "... computing ..."
for i in range(NTmaxi)[1:]:
    
    if ((t[i-1]+dt)>timeOut):
	dt = timeOut - t[i-1]
    t[i]=t[i-1]+dt

    ### Apply boundary conditions
    # From [Leveque]
    R[0]=r0*(t[i]<tf)
    
    #Update state
    if t_order==1:
        UpdateState(dt,x,R)
    elif t_order==2:
        UpdateStateRK2(dt,x,R)
    elif t_order==4:
        UpdateStateRK4(dt,x,R)
    
    
    #Store results
    rho[:,i] = R[:]
    
    if (t[i]==timeOut):
        break

"""
####Animated plot ###########################################################
# First set up the figure, the axis, and the plot element we want to animate
#Sigma
fig = plt.figure()
plt.grid()
X=np.zeros(Nnode)
for i in range(Nelem):
    X[2*i+1]=i*dx
    X[2*i+2]=(i+1)*dx
X[-1]=X[-2]

ax = plt.axes(xlim=(x[0],x[-1]), ylim=(-0.5*r0,1.5*r0))
line, = ax.plot([], [],'ro', lw=2)

time_text = ax.text(0.02, 0.95, 'middle', transform=ax.transAxes)
# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line,time_text

# animation function.  This is called sequentially
def animate(i):
    line.set_data(X,rho[:,i])
    time_text.set_text('System at time = '+str(t[i]))
    return line,time_text
 
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=rho.shape[1], interval=50, blit=True)

#Animation of the stress
plt.show()

"""
