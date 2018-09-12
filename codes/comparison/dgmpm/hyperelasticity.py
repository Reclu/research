#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import scipy.optimize as optimize
import sys
import os
sys.path.append(os.path.expandvars('dgmpm'))
from class_hyperelasticity import *
import hyperelasticExactSolution as analytic
##########################
def bar(x1,x2,Mp):
    xp=np.zeros((Mp,2))
    xp[:,0]=np.linspace(x1,x2,Mp)
    return xp

##########################


# Define geometry of the problem
L=length               # Length of the bar
Mp=Nelem*ppc             # Number of Material points
Nn=Nelem*2 + 2             # Number of elements



# Material properties

rho0=rho
Sy=Sigy           
c0=np.sqrt(C/rho0)


exact_solver=False

mesh = DGmesh(Mp,length,ppc,rho0,C,exact_solver)
dx=mesh.xn[1]-mesh.xn[0]
xp=bar(0.,length,Mp)
xp0=np.copy(xp[:,0])


mass=rho0*dx/ppc


# Define applied stress
sd= sigd              

# Define imposed gradient corresponding to sd = PK1
Res = lambda x: x + 2.*sd/C - x**3
f0=optimize.brentq(Res,np.sqrt(3.)/3.+0.00001,1.5)

qR=np.array([1.,0.])
qL=analytic.applyBoundaryCondition(C,rho0,f0,qR)
q_star = analytic.compute_stationnary_solution(qL,qR,C,rho0)

if t_order==1:
    CFL=1./ppc
elif t_order==2:
    CFL=1.
# Time discretization
dt=CFL*dx/c0 
tfinal=timeOut
tf=timeUnload;
NTMaxi=2*int(tfinal/dt)
#t_order= 1
updated_lagrangian=True
# limit = 0 : minmod // limit = 1 : superbee // limit = 2 : MUSCL
limit=-1 


T=0.
n=0

# Material points' fields
Md=mass*np.eye(Mp,Mp)
# Conserved quantities vector
U = np.zeros((Mp,2))
# Auxilary variables vector (for boundary conditions)
W = np.zeros((Mp,2))

U[:,0]=1./rho0
## Initial velocity
U[0:Mp/2,1]=v0
U[Mp/2:Mp,1]=-v0

# Conserved quantities vector
u = np.zeros((Nn,2))
# Auxilary variables vector (for boundary conditions)
w = np.zeros((Nn,2))

# Storage
Pi=np.zeros((Mp,int(NTMaxi)))
Pi_th=np.zeros((Mp,int(NTMaxi)))
V_th=np.zeros((Mp,int(NTMaxi)))
velo=np.zeros((Mp,int(NTMaxi)))
pos=np.zeros((Mp,int(NTMaxi)))
time=np.zeros(int(NTMaxi))
velo[:,0]= np.copy(U[:,1])
pos[:,0]=np.copy(xp[:,0])
time[0]=T


# Build approximation matrices
Map,Grad,Dofs,parent=mesh.buildApproximation(xp)


mg=np.dot(np.dot(Map[Dofs,:],Md),Map[Dofs,:].T)
md=np.diag(np.sum(mg,axis=1))
mass_vector0 = np.dot(np.dot(Map,Md),Map.T)
mass_vector0 = np.sum(mass_vector0,axis=1)
K0=np.dot(np.dot(Grad[Dofs,:],Md),Map[Dofs,:].T)
mass_vector=np.copy(mass_vector0)
mesh.setMapping(K0)


def plotStress(sig,color):
    X=np.zeros(len(sig))
    for i in range(Nelem+1):
        X[2*i]=X[2*i+1]=i*dx
    plt.plot(X,sig,color,label='$\sigma$',lw=2)
    plt.legend()
    plt.xlabel('x',fontsize=16)
    plt.ylabel('$\sigma$',fontsize=16)
    plt.grid()
    plt.title('Stress along the bar',fontsize=22)
    plt.show()

def UpdateState(dt,dofs,Ml,U,W,md,limit):
    Nnodes = np.shape(U)[0]
    if limit!=-1 : boolean=True
    else : boolean=False
    f=mesh.computeFlux(U,W,dofs,boolean,limit)
    for i in range(len(dofs)):
        if md[dofs[i]]!=0.:
            U[dofs[i],:]+=dt*f[dofs[i],:]/md[dofs[i]]
    return U

def UpdateStateRK2(dt,dofs,Ml,u,w,md,sd):
    Nnodes = np.shape(u)[0]
    k1=np.zeros((Nnodes,2))
    k2=np.zeros((Nnodes,2))
    # first step : compute flux and intermediate state w
    f=mesh.computeFlux(u,w,dofs,False,-1)
    for i in range(len(dofs)):
        k1[dofs[i],:]+=dt*f[dofs[i],:]/(2.*md[dofs[i]])
    u12 = u+k1
    w12 = np.zeros((Nnodes,2))

    u12[0,0]=np.copy(u12[1,0])
    u12[-2,0] = np.copy(u12[-1,0])
    
    w12[1:-1,0] =computePK1(C,rho0,u12[1:-1,0]*rho0)
    w12[:,1] = np.copy(u12[:,1])
    w12[0,0]=2.*sd - w12[1,0]
    w12[0,1]=np.copy(w12[1,1])
    #w12[-2,0] = np.copy(w12[-1,0])
    #w12[-2,1] = -np.copy(w12[-1,1])
    # second step : compute flux and update U
    f=mesh.computeFlux(u12,w12,dofs,False,-1)

    for i in range(len(dofs)):
        k2[dofs[i],:]+=dt*f[dofs[i],:]/md[dofs[i]]
    u+=k2
    return u

def computeTimeStep(rho,C,J,dx,CFL):
    c=np.sqrt(C*(3.*J**2 -1.)/(2.*rho))
    dt = CFL*dx/c
    return dt

def computePK1(C,rho0,F):
    pi = C*F*(F**2-1.)/2.
    return pi

def computeCelerity(C,rho0,F):
    cel=np.zeros(len(F))
    for i in range(len(F)):
        cel[i]=np.sqrt(C*(3.*F[i]**2-1.)/(2.*rho0))
    return cel

pi0 = computePK1(C,rho0,U[:,0]*rho0)

u0=computePK1(C,rho0,np.ones(Nn))


for n in range(NTMaxi)[1:]:    
    
    # Mapping from material points to nodes
    Um=np.dot(Md,U)
    Wm=np.dot(Md,W)
    for i in range(len(Dofs)):
        if mass_vector[Dofs[i]]!=0.:
            u[Dofs[i],:]=np.dot(Map[Dofs[i],:],Um)/mass_vector[Dofs[i]]
            w[Dofs[i],:]=np.dot(Map[Dofs[i],:],Wm)/mass_vector[Dofs[i]]

    
    # Apply load on first node
    u[2*parent[0],0]=np.copy(u[2*parent[0]+1,0])
    w[2*parent[0],0]=2.*sd*(time[n-1]<tf) - w[2*parent[0]+1,0]
    w[2*parent[0],1]=np.copy(w[2*parent[0]+1,1])
    
    # Reflective boundary conditions
    u[2*parent[-1]+3,0] = np.copy(u[2*parent[-1]+2,0])
    w[2*parent[-1]+3,0] = -np.copy(w[2*parent[-1]+2,0])
    w[2*parent[-1]+3,1] = np.copy(w[2*parent[-1]+2,1])
    

    
    Jmax=np.max(u[:,0]*rho0)
    dt=computeTimeStep(rho0,C,Jmax,dx,CFL)

    if ((time[n-1]+dt)>tfinal):
        dt = tfinal - time[n-1]
    time[n]=time[n-1]+dt
    
    if t_order==1 :
        u=UpdateState(dt,Dofs,md,u,w,mass_vector0,limit)
    elif t_order==2 :
        u=UpdateStateRK2(dt,Dofs,md,u,w,mass_vector0,sd)
    
    u0=computePK1(C,rho0,u[:,0]*rho0)
    
    # Mapping back to the material points
    U=np.dot(Map[Dofs,:].T,u[Dofs,:])
    
    xp[:,0]+=dt*U[:,1]
    
    if updated_lagrangian :
        # Compute new mapping (convective phase)
        Map,Grad,Dofs,parent=mesh.buildApproximation(np.asmatrix(xp))
        K=np.dot(np.dot(Grad[Dofs,:],Md),Map[Dofs,:].T)
        mesh.setMapping(K)
        mg=np.dot(np.dot(Map[Dofs,:],Md),Map[Dofs,:].T)
        md=np.diag(np.sum(mg,axis=1))
        mass_vector = np.dot(np.dot(Map,Md),Map.T)
        mass_vector = np.sum(mass_vector,axis=1)
        
    
    u= np.zeros((Nn,2))
    W[:,0] = computePK1(C,rho0,U[:,0]*rho0)
    W[:,1] = np.copy(U[:,1])
    
    #print 'Increment =', n, 't = ', time[n],' s.'
    velo[:,n]=np.copy(W[:,1])
    Pi[:,n]=np.copy(W[:,0])
    pos[:,n]=xp[:,0]
    
    for i in range(Mp):
        Pi_th[i,n],V_th[i,n]= analytic.computeSolution(time[n],pos[i,0],qL,qR,q_star,C,rho0)
    
    if (time[n]==tfinal):
        increments=n
        break

x=mesh.xn
time=time[0:increments]

#Sigma
fig = plt.figure()
plt.grid()
ax = plt.axes(xlim=(0.,L), ylim=(1.1*np.min(Pi),1.1*np.max(Pi)))
ax.plot(mesh.xn,np.zeros(len(mesh.xn)),'b+', lw='2.', ms='8.')
lineList = []
line1, = ax.plot([], [],'r+', lw='2.', ms='8.',label='DG_MPM')
lineList.append(line1)
line2, = ax.plot([], [],'ro', lw='2.',ms='8.',label='MPM')
lineList.append(line2)

fig.legend((lineList),('DG_MPM','MPM'),'upper right',numpoints=1)

time_text = ax.text(0.05, 0.95,'', transform=ax.transAxes)

ax.set_xlabel('x (m)', fontsize=18)
ax.set_ylabel(r'$\sigma$ (Pa)', fontsize=18)
ax.set_title('Elastic wave in 1D bar')

# initialization function: plot the background of each frame
def init():
    time_text.set_text('')
    for line in lineList:
        line.set_data([], [])
    return tuple(lineList)+(time_text,)

# animation function.  This is called sequentially
def animate(i):
    #lineList[0].set_data(DGMPM["Pos"][:,i],DGMPM["Pi"][:,i])
    lineList[1].set_data(pos[:,i],Pi[:,i])
    time_text.set_text('Stress ')
    return tuple(lineList)+(time_text,)

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=increments, interval=20, blit=True)
#Animation of the stress
plt.grid()
#anim.save('StressBar.mp4', extra_args=['-vcodec', 'libx264'])
plt.show()
