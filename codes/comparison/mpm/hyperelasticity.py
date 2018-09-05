#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import *
from pylab import *
import scipy.optimize as optimize

L=length
Mp=Nelem*ppc
Nn=Mp/ppc +1       # Number of nodes
lx=L*ppc/(Mp-1)    # Length of cells

############ METHODES
def buildMesh(Mp,l,ppc):
    nex = Mp/ppc
    nnx=nex+1
    lmp=l/(Mp-1)
    dx = ppc*l/nex
    xn = np.linspace(-lmp/2,l+lmp/2,nex+1)
    connect = np.array([np.arange(0,nnx-1,1),np.arange(1,nnx,1)]).T
    return xn,connect

def buildApproximation1D(xp,xn,connect):
    Nn=len(xn)
    Np=np.shape(xp)[0]
    Dofs=np.zeros((Nn))
    Parent=np.zeros(Np)
    matpoints=np.zeros(np.shape(connect)[0])
    Map=np.zeros((Nn,Np))
    Grad=np.zeros((Nn,Np))
    Le=(xn[1]-xn[0])
    for Pt in range(Np):
        # detect parent element of current material point
        parent=np.ceil((xp[Pt,0]-xn[0])/Le)-1 
        #print Pt,parent
        d=np.array([connect[int(parent),:]]).astype(np.int64)
        # Active nodes' indices storage
        Dofs[d[0]]+=1
        Parent[Pt]=parent+1
        xi=(xp[Pt,0]-xn[d[0][0]])/Le;
        N=np.array([1-xi,xi])
        dN=np.array([-1/Le,1/Le])
        Map[d,Pt]+=N.T;
        Grad[d,Pt]+=dN.T;
    # Number of material point in each element to compute the particles' volumes
    for elem in range(np.shape(connect)[0]):
        matpoints[elem]=np.shape(Parent[Parent==elem+1])[0]
    d=[]
    for j in range(np.shape(Dofs)[0]):
        if Dofs[j]!=0:
            d.append(j)
    return Map,Grad,d

def celerity(E,rho0,J):
    return np.sqrt(E*(3.*J**2 -1.)/(2.*rho0))

################## END OF METHODES


xn,connect=buildMesh(Mp,L,ppc)
xp=np.zeros((Mp,2))
xp[:,0]=np.linspace(0.,L,Mp)

# Material properties
Sy  = Sigy
c=np.sqrt(C/rho)

m=rho*lx/ppc

# Define force
s0 =sigd

# Define imposed specific gradient
# Define imposed gradient corresponding to sd = PK1
Res = lambda x: x + 2.*sigd/C - x**3
f0=optimize.brentq(Res,np.sqrt(3.)/3.+0.00001,1.5)
#f0=np.sqrt(2.*s0/C +1)          
if sigd >0. :
    # The fastest wave is the shock
    c_dt =  celerity(C,rho,f0)
else :
    # The fastest wave is given by the external edge of the rarefaction fan
    c_dt =  celerity(C,rho,1.)


# Time discretization
CFL/=factor
print "mpm cfl=",CFL
tfinal=timeOut   #1.*L/c
tf=timeUnload    #1.2e-3
Dt=CFL*lx/c_dt
NTMaxi=2*int(tfinal/Dt)

T=0.
n=0
A=np.zeros(Mp)
V=np.zeros(Mp)


V[0:Mp/2]=v0
V[Mp/2:Mp]=-v0

Sig=np.zeros(Mp)
Map=np.zeros((Nn,Mp))
Grad=np.zeros((Nn,Mp))
Md=m*np.eye(Mp,Mp)

a=np.zeros(Nn)
v=np.zeros(Nn)

disp=np.zeros(Nn)
Fe=np.zeros(Nn)
Fi=np.zeros(Nn)

# Storage
Pi=np.zeros((Mp,NTMaxi))
Def=np.zeros((Mp,NTMaxi))
Sth=np.zeros((Mp,NTMaxi))
Velocity=np.zeros((Mp,NTMaxi))
Pos=np.zeros((Mp,NTMaxi))
time=np.zeros(NTMaxi)
Pi[:,0]=Sig
Velocity[:,0]=V
Pos[:,0]=np.copy(xp[:,0])
time[0]=T

# Build approximation matrices
Map,Grad,Dofs=buildApproximation1D(xp,xn,connect)

mg=np.dot(np.dot(Map,Md),Map.T)
md=np.diag(np.sum(mg,axis=1))
mv=np.sum(np.dot(np.dot(Map,Md),Map.T),axis=1)

alpha=1e0 # for effective mass matrix

def computeTimeStep(J,dx,c0,CFL):
    c=np.sqrt(J)*c0
    dt = CFL*dx/c
    return dt

J=np.ones(Mp)
Def[:,0]=np.copy(J)

algo = 'USF'
#print "DGMPM time step ",Dt*factor,C
for n in range(NTMaxi)[1:]:    
    
    # Convection
    v[Dofs]=np.dot(Map[Dofs,:],np.dot(Md,V))/mv[Dofs]
    v[-1]=0.

    
    Jmax=np.max(J)
    Dt=computeTimeStep(Jmax,lx,c,CFL)
    
    if ((time[n-1]+Dt)>tfinal):
        Dt = tfinal - time[n-1]
    time[n]=time[n-1]+Dt
    
    # Forces building
    Fe[0]=-s0*(time[n-1]<tf)
    Fi[Dofs]=-np.dot(Grad[Dofs,:],np.dot(Md,Sig))
    
    # Solve motion equation
    a[Dofs]=(Fe[Dofs]+Fi[Dofs])/mv[Dofs]
    a[-1]=0.
    v+=a*Dt
    v[-1]=0.
    
    # Gradient and constitutive model
    dF=Dt*np.dot(Grad[Dofs,:].T,v[Dofs])
    Def[:,n]=(1.+dF)*Def[:,n-1]
    
    for i in range (Mp):
        Sig[i]=C*(Def[i,n]**3 -Def[i,n])/(2.*rho) 
        
    
    if s0==0.:
        Sig[0]=0.
        Sig[-1]=0.
    
    # Lagrangian step
    A=np.dot(Map[Dofs,:].T,a[Dofs])
    V+=Dt*A
    #V=np.dot(Map[Dofs,:].T,v[Dofs])
    #V[-1]=0.
    
    #xp[:,0]+=V*Dt
    xp[:,0]+=np.dot(Map[Dofs,:].T,v[Dofs])*Dt
    """
    # Compute new mapping (convective phase)
    Map,Grad,Dofs=buildApproximation1D(np.asmatrix(xp),xn,connect)
    
    mg=np.dot(np.dot(Map[Dofs,:],Md),Map[Dofs,:].T)
    md=np.diag(np.sum(mg,axis=1))
    """
    a=np.zeros(Nn)
    v=np.zeros(Nn)

    Velocity[:,n]=V
    Pi[:,n]=rho*Sig
    Pos[:,n]=np.copy(xp[:,0])
    #Def[:,n]=np.copy(J)
    
    if (time[n]==tfinal):
        increments=n
        break
    
"""
plt.plot(time[:n],Velocity[-1,:n])
plt.grid()
plt.show()
"""
factor=1
pi=Pi[:,0:increments:factor]
defo=Def[:,0:increments:factor]
pos=Pos[:,0:increments:factor]
velo=Velocity[:,0:increments:factor]
temps=time[0:increments:factor]
time=temps[0:increments:factor]

