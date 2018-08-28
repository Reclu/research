#!/usr/bin/python

import numpy as np
from pylab import *
import sys
import os
sys.path.append(os.path.expandvars('dgmpm'))
from class_elastoplasticity import *

##########################
def bar(x1,x2,Mp):
    xp=np.zeros((Mp,2))
    xp[:,0]=np.linspace(x1,x2,Mp)
    return xp

## Kinematic hardening
def bilinear(eps,EPn,Pn,E,Sigy,H):
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
        if (f<=0.):
            #elastic step
            S[i] = Selas
            EP[i] = EPn[i]
            P[i] = Pn[i]
            TM[i] = E
        elif (f>0.):
            #elastoplastic step: solve a nonlinear scalar equation
            dP = f/(E+H)
            P[i] = Pn[i]+dP
            EP[i] = EPn[i]+(P[i]-Pn[i])*np.sign(Selas)
            S[i] = E*(DEFO-EP[i])
            TM[i] = (E*H)/(E+H)
    return S,P,EP,TM

def computeCourantNumber(order,parent,Map):
    ## CFL number for Euler algorithm
    sol=[]
    # sum over mesh cells to compute the lowest Courant number by solving the minimum amplification factor equation
    
    cell2=np.where(parent==1)[0]
    for alpha in cell2:
        if parent[alpha]==0:
            continue
        # Number of material points in the same cell
        matpointN = np.where(parent==parent[alpha])[0]
        Nmp = len(matpointN)
        # Number of material points in cell i-1
        matpointM = np.where(parent==parent[alpha]-1)[0]
        Nmpp = len(matpointM)
        # Number of material points in cell i+1
        matpointL = np.where(parent==parent[alpha]+1)[0]
        L = len(matpointL)
        # nodes indices
        n1 = 2*parent[alpha]+1 ; n2=2*parent[alpha]+2
        S1 = Map[n1,matpointN] ; S2 = Map[n2,matpointN]
        Sum1=np.sum(S1) ; Sum2=np.sum(S2)
        
        Sp1 = Map[n1-2,matpointM] ; Sp2= Map[n1-1,matpointM] 
        Sump1=np.sum(Sp1) ; Sump2=np.sum(Sp2)

        Courant = symbols('Courant')
        Res=0.
        Dmu=0.
        for mu in range(Nmp): # sum over particles contained in the same cell
            Dmu = Map[n1,alpha]*S1[mu]/Sum1 + Map[n2,alpha]*S2[mu]/Sum2 + Courant*( Map[n2,alpha]/Sum2 - Map[n1,alpha]/Sum1 -Nmp*Map[n2,alpha]*S2[mu]/(Sum2**2) )
            ## Second order contributions
            if order==2:
                Dmu += 0.5*Nmp*(Courant**2)*((S2[mu]/Sum2)*(Map[n1,alpha]/Sum1-Map[n2,alpha]/Sum2) + (Map[n2,alpha]/(Sum2**2))*(Nmp*S2[mu]/Sum2-1.) )
            Res+=np.abs(Dmu)
        for mu in range(Nmpp): # sum over particles contained in previous cell
            Dmu=Courant*Nmp*Map[n1,alpha]*Sp2[mu]/(Sum1*Sump2)
            ## Second order contributions
            if order==2:
                Dmu +=0.5*Nmp*(Courant**2)*( Map[n1,alpha]/(Sum1*Sump2)*(1.-Nmpp*Sp2[mu]/Sump2) -(Sp2[mu]/Sump2)*(Map[n1,alpha]/Sum1-Map[n2,alpha]/Sum2) )
            Res+=np.abs(Dmu)
        Residual = lambdify((Courant),Res-1.)
        #sol.append(optimize.root(Residual,1.,method='hybr',options={'xtol':1.e-12}).x)
        sol.append(optimize.newton(Residual,1.))
    solution_CFL=np.min(sol)
    print "Courant number set to :",solution_CFL
    return solution_CFL

##########################


# Define geometry of the problem
L=length 
Mp=Nelem*ppc          
Nn=Nelem*2 + 2



# Material properties
E=Young
HT=E*H/(E+H)
cp=np.sqrt(HT/rho)
c=np.sqrt(E/rho)
Sy=Sigy

mesh = DGmesh(Mp,L,ppc,c,cp,rho,Sy,H)
dx=mesh.xn[1]-mesh.xn[0]
xp=bar(0.,L,Mp)


mass=rho*dx/ppc

# Define applied stress
sd=sigd     

# Define imposed specific stress
s0=sd/rho


# Time discretization
if compute_CFL:
    CFL=computeCourantNumber(t_order,parent,Map)
else:
    if t_order==2:
        CFL=1.
    else:
        CFL=float(t_order)/float(ppc)
        
Dt=CFL*dx/c 
tfinal=timeOut
tf=timeUnload
inc=round(tfinal/Dt)

limit=False

T=0.
n=0

# Material points' fields
Md=mass*np.eye(Mp,Mp)
U = np.zeros((Mp,2))


U[0:Mp/2,1]=v0
U[Mp/2:Mp,1]=-v0

# Auxilary variables vector (for boundary conditions)
W = np.copy(U)
W[0:Mp/2,1]=v0
W[Mp/2:Mp,1]=-v0

# Nodes' fields
u = np.zeros((Nn,2))
ep = np.zeros(Nn)
# Auxilary variables vector (for boundary conditions)
w = np.zeros((Nn,2))

# Storage
sig=np.zeros((Mp,int(inc)+2))
epsp = np.zeros((Mp,int(inc)+2))
p=np.zeros((Mp,int(inc)+2))
Epsilon=np.zeros((Mp,int(inc)+2))
dEpsp=np.zeros((Mp,int(inc)+2))

Velocity=np.zeros((Mp,int(inc)+2))
pos=np.zeros((Mp,int(inc)+2))
time=np.zeros(int(inc)+2)
sig[:,0]=U[:,0]
Velocity[:,0]= U[:,1]
pos[:,0]=xp[:,0]
time[0]=T

# Build approximation matrices
Map,Grad,Dofs,parent=mesh.buildApproximation(xp)

#mesh.set_sigy(np.dot(Map,np.ones(Mp)*Sigy))
mg=np.dot(np.dot(Map[Dofs,:],Md),Map[Dofs,:].T)
md=np.diag(np.sum(mg,axis=1))
mass_vector = np.dot(np.dot(Map,Md),Map.T)
mass_vector = np.sum(mass_vector,axis=1)
K=np.dot(np.dot(Grad[Dofs,:],Md),Map[Dofs,:].T)
alpha=1e0 # for effective mass matrix

dim=2 # Number of unknowns
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

def UpdateState(dt,dofs,Ml,U,W,md,ep,limit):
    Nnodes = np.shape(U)[0]
    f=mesh.computeFlux(U,W,dofs,md,ep)
    for i in range(len(dofs)):
        if md[dofs[i]]!=0.:
            U[dofs[i],:]+=dt*f[dofs[i],:]/md[dofs[i]]
    return U

def UpdateStateRK2(dt,dofs,Ml,U,W,md,ep):
    Nnodes = np.shape(U)[0]
    k1=np.zeros((Nnodes,2))
    k2=np.zeros((Nnodes,2))
    # first step : compute flux and intermediate state w
    f=mesh.computeFlux(U,W,dofs,md,ep)
    for i in range(len(dofs)):
        if md[dofs[i]]!=0.:
            k1[dofs[i],:]+=dt*f[dofs[i],:]/md[dofs[i]]
    u12 = U+k1*0.5
    w12 = np.copy(u12)
    w12[:,0] = rho*u12[:,0]*E
    w12[:,1] = u12[:,1]
    
    w12[0,:] = prescribeStress(w12[1,:],[ep[0],ep[1]],rho,(c,c),(cp,cp),(Sy,Sy),(H,H),sd*(T<=timeUnload))
    # second step : compute flux and update U
    f=mesh.computeFlux(u12,w12,dofs,md,ep)
    for i in range(len(dofs)):
        if md[dofs[i]]!=0.:
            k2[dofs[i],:]+=dt*f[dofs[i],:]/md[dofs[i]]
    U+=k2
    return U

while T<tfinal:
    
    # Mapping from material points to nodes
    
    Um=np.dot(Md,U)
    Wm=np.dot(Md,W)
    EPm=np.dot(Md,p[:,n])
    for i in range(len(Dofs)):
        if mass_vector[Dofs[i]]!=0.:
            u[Dofs[i],:]=np.dot(Map[Dofs[i],:],Um)/mass_vector[Dofs[i]]
            w[Dofs[i],:]=np.dot(Map[Dofs[i],:],Wm)/mass_vector[Dofs[i]]
            ep[Dofs[i]]=np.dot(Map[Dofs[i],:],EPm)/mass_vector[Dofs[i]]
    
    ep[0]=ep[1] ; ep[-1]=ep[-2]

    s0=0.
    # Apply load on first node
    w[2*parent[0],0]= 2.*s0*(T<tf) - w[2*parent[0]+1,0] 
    w[2*parent[0],1]= w[2*parent[0]+1,1]
    # Transmissive boundary conditions
    w[2*parent[-1]+3,0] = 2.*s0*(T<tf) - w[2*parent[-1]+2,0]
    w[2*parent[-1]+3,1] = w[2*parent[-1]+2,1]
    
    if t_order==1 :
        u=UpdateState(Dt,Dofs,md,u,w,mass_vector,ep,limit)
    elif t_order==2 :
        u=UpdateStateRK2(Dt,Dofs,md,u,w,mass_vector,ep)
    
    # Mapping back to the material points
    U=np.dot(Map.T,u)
    
    u = np.zeros((Nn,2))
    w = np.zeros((Nn,2))

    Eps=U[:,0]*rho
    Sig,p[:,n+1],epsp[:,n+1],tangent_modulus = bilinear(Eps,epsp[:,n],p[:,n],E,Sy,H)


    #print 'Increment =', n, 't = ', T,' s.'
    n+=1
    T+=Dt
    Velocity[:,n]=U[:,1]
    sig[:,n]=Sig
    Epsilon[:,n]=U[:,0]*rho
    W[:,0]=Sig
    W[:,1]=np.copy(U[:,1])

    pos[:,n]=xp[:,0]
    time[n]=T
    
increments=n
