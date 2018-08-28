# !/usr/bin/python

from pylab import *
from matplotlib import animation
from matplotlib import pyplot as plt
from scipy.linalg import solve
import numpy as np
import math 
import pdb


L=length
Mp=Nelem*ppc
Nn=Nelem + 1             # Number of elements
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



def bilinear(eps_n,eps,EPn,Pn,E,Sigy,H):
    #initialization
    S = np.zeros(len(eps))
    EP = np.zeros(len(eps))
    P  = np.zeros(len(eps))
    TM  = np.zeros(len(eps))
    #Loop on integration points
    for i,DEFO in enumerate(eps):
        #(i) Elastic prediction
        Selas = E*(DEFO-EPn[i])
        #(ii) Compute the criterion 
        f = np.abs(Selas) - (Sigy+H*Pn[i])
        if (f<=0):
            #elastic step
            S[i] = Selas
            EP[i] = EPn[i]
            P[i] = Pn[i]
            TM[i] = E
        elif (f>0):
            #elastoplastic step: solve a nonlinear scalar equation
            dP = f/(E+H)
            P[i] = Pn[i]+dP
            EP[i] = EPn[i]+(P[i]-Pn[i])*np.sign(Selas)
            S[i] = E*(DEFO-EP[i])
            TM[i] = (E*H)/(E+H)
    return S,P,EP,TM

################## END OF METHODES

xn,connect=buildMesh(Mp,L,ppc)

xp=np.zeros((Mp,2))
xp[:,0]=np.linspace(0.,L,Mp)


# Material properties
E=Young
c = math.sqrt(E/rho)
m = rho*lx/ppc

#Bc_nodes=np.array([0]) # left extremity fixed
Bc_nodes=[Nn-1]

s0 =-sigd

# Time discretization
#factor=2.
CFL/=factor

tfinal=timeOut   #1.*L/c
tf=timeUnload    #1.2e-3
Dt=CFL*lx/c

inc=int(tfinal/Dt)

T=0.
n=0

# Material points' fields
A=np.zeros(Mp)
V=np.zeros(Mp)
#X=xp[:,0]
X=xp[:,0]
Sig=np.zeros(Mp)
Eps=np.zeros(Mp)
Fimp=np.zeros(Mp)
Fimp[0]=s0

# Grid nodes' fields
a=np.zeros(Nn)
v=np.zeros(Nn)
Fe=np.zeros(Nn)
Fi=np.zeros(Nn)

# Approximation builing
Map,Grad,Dofs=buildApproximation1D(np.asmatrix(xp),xn,connect)
#Map,Grad,Dofs=buildGIMPApproximation1D(np.asmatrix(xp),xn,lx/2)
MD=m*np.eye(Mp)

#v0=s0/(rho*c)
V[0:Mp/2]=v0
V[Mp/2:Mp]=-v0


# Nodal mass matrix
mn=np.dot(np.dot(Map[Dofs,:],MD),Map[Dofs,:].T)
# Lumped mass matrix
md=np.diag(np.sum(mn,axis=1))
mv=np.sum(np.dot(np.dot(Map,MD),Map.T),axis=1)

# Storage
Stress=np.zeros((Mp,int(inc)+2))
Sth=np.zeros((Mp,int(inc)+2))
Vth=np.zeros((Mp,int(inc)+2))
p=np.zeros((Mp,int(inc)+2))
Epsp=np.zeros((Mp,int(inc)+2))
dEpsp=np.zeros((Mp,int(inc)+2))
Pos=np.zeros((Mp,int(inc)+2))
Velocity=np.zeros((Mp,int(inc)+2))
time=np.zeros(int(inc)+2)
NRG=np.zeros(int(inc)+2)
kin=np.zeros(int(inc)+2)
strain=np.zeros(int(inc)+2)
Stress[:,n]=Sig[:]
Pos[:,n]=X[:]
Velocity[:,n]=V[:]
time[n]=T

Sig,p[:,0],Epsp[:,0],tangent_modulus = bilinear(np.zeros(len(Eps)),Eps,Epsp[:,0],p[:,0],E,Sigy,H)

alg=algo
alpha=1.
#if ppc==2: alpha=0.5

n=0
kin[n]=0.5*np.inner(np.dot(MD,V),V)
strain[n]=0.5*np.inner(np.dot(MD,Sig/rho),Eps)
NRG[n]=kin[n]+strain[n]
Vth[0:Mp/2,0]=v0
Vth[Mp/2:Mp,0]=-v0

while T<tfinal:
    #pdb.set_trace()
    mf=(1-alpha)*mn + alpha*md

    if alg=='USL':
        
        # Convection /!\
        v[Dofs]=np.dot(Map[Dofs,:],np.dot(MD,V))/mv[Dofs]
        #v[Dofs]=solve(mf[Dofs],np.dot(Map[Dofs,:],np.dot(MD,V)))
        
        # Internal force vector
        Fe[0]=s0*(T<tf)
        Fi[Dofs]=-np.dot(Grad[Dofs,:],np.dot(MD,Sig))/rho

        # Solve motion equation
        a[Dofs]=(Fe[Dofs]+Fi[Dofs])/mv[Dofs]
        #a[Dofs]=solve(mf[Dofs],Fe[Dofs]+Fi[Dofs])
        
        #Time integrator
        v+=Dt*a
           
        # Gradient and constitutive model
        Epsn=Eps
        Eps+=Dt*np.dot(Grad[Dofs,:].T,v[Dofs])
        Sig,p[:,n],Epsp[:,n],tangent_modulus = bilinear(Epsn,Eps,Epsp[:,n],p[:,n],E,Sigy,H)
    elif alg=='USF':
        # Convection /!\
        v[Dofs]=solve(mf,np.dot(Map[Dofs,:],np.dot(MD,V)))
        
        # Gradient and constitutive model
        Epsn=Eps
        Eps+=Dt*np.dot(Grad[Dofs,:].T,v[Dofs])
        Sig,p[:,n],Epsp[:,n],tangent_modulus = bilinear(Epsn,Eps,Epsp[:,n],p[:,n],E,Sigy,H)
        
        # Internal force vector
        Fi[Dofs]=-np.dot(Grad[Dofs,:],np.dot(MD,Sig))/rho
        Fe[0]=s0*(T<tf)
        
        # Solve motion equation
        a[Dofs]=solve(mf,Fi[Dofs]+Fe[Dofs])
        
        # Time integrator
        v+=Dt*a
        
    # Lagrangian step 
    if mpm_mapping :
        A=np.dot(Map[Dofs,:].T,a[Dofs])
        V+=Dt*A
    else :
        V=np.dot(Map[Dofs,:].T,v[Dofs])    


    if update_position:
        xp[:,0]+=Dt*V
        Map,Grad,Dofs=buildApproximation1D(np.asmatrix(xp),xn,connect)
        MD=m*np.eye(Mp)

        # Nodal mass matrix
        mn=np.dot(np.dot(Map[Dofs,:],MD),Map[Dofs,:].T)
        # Lumped mass matrix
        md=np.diag(np.sum(mn,axis=1))
        mv=np.sum(np.dot(np.dot(Map,MD),Map.T),axis=1)

    n+=1
    T+=Dt
    Stress[:,n]=Sig[:]
    for i in range(Mp):
        Sth[i,n]=-rho*c*v0*(abs(xp[i,0]-L/2.)<c*T)
        Vth[i,n]=-v0*(abs(xp[i,0]-L/2.)>c*T)*np.sign(xp[i,0]-L/2.)
    Pos[:,n]=xp[:,0]
    Velocity[:,n]=V[:]
    time[n]=T
    
    kin[n]=0.5*np.inner(np.dot(MD,V),V)
    strain[n]=0.5*np.inner(np.dot(MD,Sig/rho),Eps)
    NRG[n]=kin[n]+strain[n]
    increments=n
factor=1#int(factor)
#factor=1
sig=Stress[:,0:-1:factor]
velo=Velocity[:,0:-1:factor]
epsp=Epsp[:,0:-1:factor]
pos=Pos[:,0:-1:factor]
