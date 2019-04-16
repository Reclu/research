#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import *
from pylab import *
#from scipy.linalg import solve
from DGmesh import *

##########################
def bar(x1,x2,Mp):
    xp=np.zeros((Mp,2))
    xp[:,0]=np.linspace(x1,x2,Mp)
    return xp

##########################


#print 'Initializing problem ...'
# Define geometry of the problem
L=1.               # Length of the bar
Nelem=Mp/ppc             # Number of elements

Nn=Nelem*2 + 2             # Number of elements

 
# Material properties

rho=7800.
E=2.e11
Sy=400.0e6           
c=math.sqrt(E/rho)


#print '       Mesh Definition'

mesh = DGmesh(Mp,L,ppc,c,rho)
dx=mesh.xn[1]-mesh.xn[0]
xp=bar(0.,L,Mp)

shift=0.*dx
xp[:,0]+=shift
#mesh.xn+=0.01
"""
plt.plot(xp[:,0],xp[:,1],'ro',label='Material points')
plt.plot(mesh.xn,np.zeros(len(mesh.xn)),'b+',label='Nodes')
plt.axis('equal')
plt.legend(loc='best',numpoints=1)
plt.show()
"""

mass=rho*dx/ppc


# Boundary condition
R0=1.e2               

# Define imposed specific quantity
r0=R0/rho               

#print '       Algorithmic parameters'
# Time discretization

Dt=CFL*dx/c 
tfinal=1.*L/c
tf=0.5*tfinal;
inc=round(tfinal/Dt)
t_order= 1

T=0.
n=0

# Material points' fields
Md=mass*np.eye(Mp,Mp)
U = np.zeros(Mp)


# Nodes' fields
u = np.zeros(Nn)

# Storage
Stress=np.zeros((Mp,int(inc)+2))
analytical=np.zeros((Mp,int(inc)+2))
FOU=np.zeros((Mp,int(inc)+2))
time=np.zeros(int(inc)+2)
Stress[:,0]=U[:]
time[0]=T

# Build approximation matrices
Map,Grad,Dofs,parent=mesh.buildApproximation(xp)


mg=np.dot(np.dot(Map[Dofs,:],Md),Map[Dofs,:].T)
md=np.diag(np.sum(mg,axis=1))
mass_vector = np.dot(np.dot(Map,Md),Map.T)
mass_vector = np.sum(mass_vector,axis=1)
K=np.dot(np.dot(Grad[Dofs,:],Md),Map[Dofs,:].T)
alpha=1.e0 # for effective mass matrix

dim=1 # Number of unknowns
mesh.setMapping(K)

def plotStress(sig,color):
    X=np.zeros(len(sig))
    for i in range(Nelem+1):
        X[2*i]=X[2*i+1]=i*dx
    plt.plot(X,sig,color,label='$\sigma$',lw=2)
    plt.plot(X,s0*np.ones(len(X)),'m--',label='applied stress',lw=2)
    plt.legend()
    plt.xlabel('x',fontsize=16)
    plt.ylabel('$\sigma$',fontsize=16)
    plt.grid()
    plt.title('Stress along the bar',fontsize=22)
    plt.show()

def UpdateState(dt,dofs,M,U,md,dim,limiter):
    Nnodes = np.shape(U)[0]
    f=mesh.computeFlux(U,dofs,md,dim,limiter)
    U[dofs]+=dt*solve(M,f[dofs])
    #return U
    return U,dt*solve(M,f[dofs])


def UpdateStateRK2(dt,dofs,Ml,U,md,dim):
    Nnodes = np.shape(U)[0]
    k1=np.zeros(Nnodes)
    k2=np.zeros(Nnodes)
    # first step : compute flux and intermediate state w
    f=mesh.computeFlux(U,dofs,md,dim)
    k1[dofs]+=dt*(solve(Ml,f[dofs]))
    w = U+k1*0.5
    # second step : compute flux and update U
    f=mesh.computeFlux(w,dofs,md,dim)
    k2[dofs]+=dt*(solve(Ml,f[dofs]))
    U+=k2
    return U

def loadingFunction(t,tfinal,tf):
    return np.sin(8.*np.pi*t/tfinal)*(t<tf)

def smoothSolution(x,c,t,tfinal,tf):
    if c*t < x:
        val=0.
    else:
        val=loadingFunction((t-x/c),tfinal,tf)
    return val

if algo=='test':
    alpha=CFL
    limiter=0
else:
    alpha=1.
    limiter=-1

print '... computing ...'
while T<tfinal:
    
    # Effective mass matrix
    mf=(1-alpha)*mg + alpha*md
    
    # Mapping from material points to nodes
    u[Dofs]=solve(mf,np.dot(Map[Dofs,:],np.dot(Md,U)))
    #u[Dofs]=np.dot(Map[Dofs,:],np.dot(Md,U))/mass_vector[Dofs]

    if sinusoidal:
        u[0]=r0*loadingFunction(T+Dt,tfinal,tf)
    else:
        u[0]=r0*(T<tf)
        
    if t_order==1 :
        u,du=UpdateState(Dt,Dofs,md,u,mass_vector,dim,limiter)
    elif t_order==2 :
        u=UpdateStateRK2(Dt,Dofs,md,u,mass_vector,dim)
    
    
    # Mapping back to the material points
    U=np.dot(Map.T,u)

    #print "Increment ",n," Time t=",T," s."
    
    n+=1
    T+=Dt
    Stress[:,n]=rho*U
    for i in range(Mp):
        if i>0:
            FOU[i,n]=FOU[i,n-1]-CFL*(FOU[i,n-1]-FOU[i-1,n-1])
        else :
            FOU[i,n]=FOU[i,n-1]-CFL*(FOU[i,n-1]-rho*r0*(T<=tf))
        if sinusoidal:
            analytical[i,n]=R0*smoothSolution(xp[i,0],c,T,tfinal,tf)
        else:
            analytical[i,n]=r0*(c*T>xp[i,0] and c*(T-tf)<=xp[i,0])*rho
        
    time[n]=T
    




