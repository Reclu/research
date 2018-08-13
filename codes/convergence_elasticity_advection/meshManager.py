# !/usr/bin/python

import numpy as np
import math as m

def buildMesh(Mp,l,ppc):
# Mesh built by giving :
# 1-Number of elements in x-direction
# 2-Length of meshed domain
# 3-Number of particle per cell
    nex = Mp/ppc
    nnx=nex+1
    lmp=l/(Mp-1)
    dx = ppc*l/nex
    xn = np.linspace(-lmp/2,l+lmp/2,nex+1)
    connect = np.array([np.arange(0,nnx-1,1),np.arange(1,nnx,1)]).T
    return xn,connect

def bar(x1,x2,Mp):
    xp=np.zeros((Mp,2))
    xp[:,0]=np.linspace(x1,x2,Mp)
    return xp
    
def circle(c,r,nr,nt):
    pi=m.pi
    xp=np.zeros(((nr-1)*nt + 1,2))
    xp[0,:]=np.array([0,0])
    dr=r/(nr-1)
    dt=2.*pi/nt
    count=1
    for t in range(nt):
        for r in range(nr-1):
           xp[count,:]=np.array([(r+1)*dr*m.cos(t*dt),(r+1)*dr*m.sin(t*dt)])
           count+=1
    xp[:]+=c
    return xp

def rectangle(x0,Nx,Ny,lx,ly):
    xp=np.zeros((Nx*Ny,2))
    dx=lx/(Nx-1)
    if Ny!=1:
        dy=ly/(Ny-1)
    else :
        dy=0
    for iy in range(Ny):
        for ix in range(Nx):
            xp[iy*Nx + ix,0]=x0[0]+ix*dx
            xp[iy*Nx + ix,1]=x0[1]+iy*dy
    return xp
