#!/usr/bin/python

import numpy as np
from pylab import *
import sys
import os
import scipy.optimize
import pdb
from matplotlib import animation
from matplotlib import pyplot as plt
from class_AcousticSolver import *
from constitutiveModels import *
from material import *
from state import *
from hardeningModels import *
from kinematics import *

##########################
def bar(x1,x2,Mp):
    xp=np.zeros((Mp,2))
    xp[:,0]=np.linspace(x1,x2,Mp)
    return xp


def constitutiveUpdate(eps,EPn,Pn,integrator):
    #initialization
    S = np.zeros(len(eps))
    EP = np.zeros(len(eps))
    P  = np.zeros(len(eps))
    #Loop on integration points
    for i,DEFO in enumerate(eps):
        zeros=np.zeros((3,3))
        state0=state(zeros,zeros,[np.array([[EPn[i],0.,0.],[0.,-EPn[i]/2.,0.],[0.,0.,-EPn[i]/2.]]),Pn[i]],integrator.hardeningModel)
        state1=state(np.array([[DEFO,0.,0.],[0.,0.,0.],[0.,0.,0.]]),zeros,[zeros,0.],integrator.hardeningModel)

        dp=integrator.constitutiveUpdate(state0,state1)
        S[i]=state1.flux[0,0]
        P[i]=state1.get("CUMULATED_PLASTIC_STRAIN")
        EP[i]=state1.get("PLASTIC_STRAIN")[0,0]
    return S,P,EP

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


### PARAMETERS
NTmaxi = 300
length = 2.0
ppc=1
Nelem = 100
E = 2.0e11
nu=0.3
Sigy = 400.0e6
#Sigy = 200.0e7
lamb = (nu*E)/(((1.0+nu)*(1.0-2.0*nu)))
mu= E/(2.0*(1.0+nu))
H = 0.18*E
#H = 10.e9
IH_power=1.
rho = 7800.0
c = np.sqrt((lamb+2.0*mu)/rho)
sigd =0.
HEL = ((lamb+2.0*mu)/(2.0*mu))*Sigy
v0=3.5*HEL/(rho*c)
#v0=9.*Sigy/(rho*c)

kin=planeWaveKinematic()
mat=plasticMaterial(E,nu,rho,Sigy)
hardening = isotropicHardening(mat,H,1./IH_power)
hardening = kinematicHardening(mat,H)
integrator=J2Plasticity(mat,kin,hardening)

timeOut = 1.25*length/(np.sqrt(E/rho))
t_order=1
timeUnload = 2*timeOut

# Define geometry of the problem
L=length 
Mp=Nelem*ppc          
Nn=Nelem*2 + 2



mesh = DGmesh(Mp,L,ppc,c,rho)
dx=mesh.xn[1]-mesh.xn[0]
xp=bar(0.,L,Mp)


mass=rho*dx/ppc

# Define applied stress
sd=sigd         

# Define imposed specific stress
s0=sd/rho

# Time discretization
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
Sth=np.zeros((Mp,int(inc)+2))
EPth = np.zeros((Mp,int(inc)+2))
p=np.zeros((Mp,int(inc)+2))
Epsilon=np.zeros((Mp,int(inc)+2))

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


def UpdateState(dt,dofs,Ml,U,W,md,limit):
    Nnodes = np.shape(U)[0]
    f=mesh.computeFlux(U,W,dofs,md)
    for i in range(len(dofs)):
        if md[dofs[i]]!=0.:
            U[dofs[i],:]+=dt*f[dofs[i],:]/md[dofs[i]]
    if limit :
        Umean=np.zeros((Nelem,2))
        # Valeur moyenne de U sur chaque element
        for i in range(Nelem):
            Umean[i,:]=(U[2*i,:]+U[2*i+1,:])*0.5
            # peut etre fait dans la boucle d'update
        vl=np.zeros(Nelem)
        vr=np.zeros(Nelem)
        for i in range(Nelem):
            # add a ghost cell if limiters used !
            vl[i]=limited_flux(Umean[i,:]-U[2*i,:],Umean[i,:]-Umean[i-1,:],Umean[i+1,:]-Umean[i,:])
            vr[i]=limited_flux(Umean[i,:]-U[2*i,:],Umean[i,:]-Umean[i-1,:],Umean[i+1,:]-Umean[i,:])
            
        # Pour chaque element, comparer avec gauche et droite
    return U

def UpdateStateRK2(dt,dofs,Ml,U,W,md):
    print "RK2 requires the computation of constitutive equations !!!!!"
    Nnodes = np.shape(U)[0]
    k1=np.zeros((Nnodes,2))
    k2=np.zeros((Nnodes,2))
    # first step : compute flux and intermediate state w
    f=mesh.computeFlux(U,dofs,md)
    for i in range(len(dofs)):
        k1[dofs[i],:]+=dt*f[dofs[i],:]/(2.*md[dofs[i]])
    w = U+k1
    # second step : compute flux and update U
    f=mesh.computeFlux(w,dofs,md)
    for i in range(len(dofs)):
        k2[dofs[i],:]+=dt*f[dofs[i],:]/md[dofs[i]]
    U+=k2
    return U



while T<tfinal:
    
    # Mapping from material points to nodes
    
    Um=np.dot(Md,U)
    Wm=np.dot(Md,W)
    for i in range(len(Dofs)):
        if mass_vector[Dofs[i]]!=0.:
            u[Dofs[i],:]=np.dot(Map[Dofs[i],:],Um)/mass_vector[Dofs[i]]
            w[Dofs[i],:]=np.dot(Map[Dofs[i],:],Wm)/mass_vector[Dofs[i]]
            
    
    # Apply load on first node
    w[2*parent[0],0]=2.*s0*(T<tf) - w[2*parent[0]+1,0] 
    w[2*parent[0],1]=w[2*parent[0]+1,1]
    # Transmissive boundary conditions
    w[2*parent[-1]+3,0] =2.*s0*(T<tf) - w[2*parent[-1]+2,0]
    w[2*parent[-1]+3,1] =w[2*parent[-1]+2,1]
    
    if t_order==1 :
        u=UpdateState(Dt,Dofs,md,u,w,mass_vector,limit)
    elif t_order==2 :
        u=UpdateStateRK2(Dt,Dofs,md,u,w,mass_vector,ep)
    
    # Mapping back to the material points
    U=np.dot(Map.T,u)
    
    u = np.zeros((Nn,2))
    w = np.zeros((Nn,2))

    
    Eps=U[:,0]*rho
    Sig,p[:,n+1],epsp[:,n+1]= constitutiveUpdate(Eps,epsp[:,n],p[:,n],integrator)
    
    #pdb.set_trace()
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

plt.figure(figsize=(3,3*2.2178))

X,Y=np.meshgrid(xp[:Mp/2,0],time[:n])
plt.contourf(X,Y,-sig[Mp/2:,:increments].T,cmap=plt.cm.binary)
point=length/2.+0.4
coor=abs(xp[:,0]-point)
index=np.where(coor==min(coor))[0][0]
plt.plot([xp[index,0]-length/2.,xp[index,0]-length/2.],[0.,T],'r')
plt.xlim((0.,length/2.))
# plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
# plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
# plt.xlabel('x (m)')
# plt.ylabel('t (s)')

#plt.legend()
# plt.tick_params(
#     axis='x',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     bottom=False,      # ticks along the bottom edge are off
#     top=False,         # ticks along the top edge are off
#     labelbottom=False) # labels along the bottom edge are off
# plt.tick_params(
#     axis='y',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     left=False,      # ticks along the bottom edge are off
#     right=False,         # ticks along the top edge are off
#     labelleft=False) # labels along the bottom edge are off

# plt.figure()
# plt.plot(time[:increments],p[index,:increments])
# plt.grid()

plt.figure()
plt.plot(time[:increments],sig[index,:increments])
plt.grid()

plt.show()





fig, (ax1,ax2) = plt.subplots(2,1)

# intialize two line objects (one in each axes)
line1, = ax1.plot([], [],'g', ms=2.5,label='dgmpm')
line2, = ax2.plot([], [],'g', ms=2.5,label='dgmpm')
line = [line1,line2]

ax1.grid()
ax2.grid()
ax1.set_xlabel('x (m)', fontsize=18)
ax1.set_ylabel(r'$\sigma$', fontsize=18)
ax2.set_xlabel('x (m)', fontsize=18)
ax2.set_ylabel(r'$p$', fontsize=18)

ax1.set_xlim(0.,length)
ax1.set_ylim(1.1*np.min(sig),1.1*np.max(sig))
ax2.set_xlim(0.,length)
ax2.set_ylim(1.1*np.min(p),1.1*np.max(p))
ax1.legend(numpoints=1)
def init():
    line[0].set_data([], [])
    line[1].set_data([], [])
    return line

def animate(i):
    line[0].set_data(pos[:,i],sig[:,i])
    line[1].set_data(pos[:,i],p[:,i])
    return line


anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=increments, interval=100, blit=True)

plt.show()

