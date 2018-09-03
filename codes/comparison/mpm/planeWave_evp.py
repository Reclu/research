# !/usr/bin/python

from pylab import *
from matplotlib import animation
from matplotlib import pyplot as plt
from scipy.linalg import solve
from scipy.optimize import fsolve
from scipy.optimize import root
import numpy as np
import pdb
import os
import sys
sys.path.append(os.path.expandvars('mpm'))


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

def bilinear(eps,Sn,EPn,Pn,lam,mu,Sigy,H,eta,n,dt):
    #initialization
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
            result = root(equations,(Sguess,EPn[i]+1.0e-7),args=(Sx_trial,EPn[i],lam,mu,Sigy,H,eta,n,dt),method='hybr')
            S[i] = result.x[0] ; EP[i] = result.x[1]
            #Post-processing
            Seq = computeS(S[i],EP[i],lam,mu,H)
            f= criterion(Seq,Sigy)
            P[i]=Pn[i]+dt*creep_law(f,Seq,eta,n)*np.sign(Seq)
    return S,P,EP
################## END OF METHODES

xn,connect=buildMesh(Mp,L,ppc)

xp=np.zeros((Mp,2))
xp[:,0]=np.linspace(0.,L,Mp)


# Material properties
E=Young
c = np.sqrt(E/rho)
m = rho*lx/ppc
power=n
lam = (nu*E)/(((1.0+nu)*(1.0-2.0*nu)))
mu = E/(2.0*(1.0+nu))
c = np.sqrt((lam+2.0*mu)/rho)       #Sound speed

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

alg=algo
alpha=1.
#if ppc==2: alpha=0.5

n=0
kin[n]=0.5*np.inner(np.dot(MD,V),V)
strain[n]=0.5*np.inner(np.dot(MD,Sig/rho),Eps)
NRG[n]=kin[n]+strain[n]
#compute stress
Sig,p[:,n],Epsp[:,n] = bilinear(np.zeros(Mp),Sig,Epsp[:,n-1],p[:,n-1],lam,mu,Sigy,H,eta,power,Dt)


while T<tfinal:
    #pdb.set_trace()
    mf=(1-alpha)*mn + alpha*md

    if alg=='USL':
        
        # Convection /!\
        v[Dofs]=np.dot(Map[Dofs,:],np.dot(MD,V))/mv[Dofs]
        
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
        Sig,p[:,n+1],Epsp[:,n+1] = bilinear(Eps,Sig,Epsp[:,n-1],p[:,n-1],lam,mu,Sigy,H,eta,power,Dt)
            
    elif alg=='USF':
        # Convection /!\
        v[Dofs]=solve(mf,np.dot(Map[Dofs,:],np.dot(MD,V)))
        
        # Gradient and constitutive model
        Epsn=Eps
        Eps+=Dt*np.dot(Grad[Dofs,:].T,v[Dofs])
        if hardening=='kinematic':
            Sig,Epsp[:,n+1]=stress_update_kin(Eps,Stress[:,n],Epsp[:,n],E,Sigy,H,eta,power,Dt)  
        elif hardening=='isotropic':
            Sig,Epsp[:,n+1],p[:,n+1]=stress_update_iso(Eps,Stress[:,n],Epsp[:,n],p[:,n],E,Sigy,H,eta,power,Dt)
            
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
    Stress[:,n]=Sig
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
