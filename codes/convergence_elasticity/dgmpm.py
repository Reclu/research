#!/usr/bin/python

import numpy as np
from pylab import *
from DGmesh import *
import matplotlib.pyplot as plt
from scipy import optimize
from sympy import *
from matplotlib import animation
import pdb

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
def bar(x1,x2,Mp):
    xp=np.zeros((Mp,2))
    xp[:,0]=np.linspace(x1,x2,Mp)
    return xp

##########################


#print 'Initializing problem ...'
# Define geometry of the problem
L=25.               # Length of the bar
Nelem=Mp            # Number of elements
Mp=Nelem*ppc
Nn=Nelem*2 + 2             # Number of elements


# Material properties

rho=1.
E=100.
Sy=400.0e6           
c=np.sqrt(E/rho)
H = 10e9
HT=E*H/(E+H)
cp=np.sqrt(HT/rho)


#print '       Mesh Definition'
#mesh = DGmesh(2*Mp,2*L,ppc,c,rho)
mesh = DGmesh(Mp,L,ppc,c,rho)
dx=mesh.xn[1]-mesh.xn[0]
xp=bar(0.,L,Mp)


"""
plt.plot(mesh.xn,np.zeros(len(mesh.xn)),'b+')
plt.plot(xp[:,0],xp[:,1],'ro')
plt.show()
"""

# Build approximation matrices
Map,Grad,Dofs,parent=mesh.buildApproximation(xp)

mass=rho*dx/ppc


# Define applied stress
sd=0.
v0=0.1

# Define imposed specific stress
s0=sd/rho               

#print '       Algorithmic parameters'
# Time discretization
CFL=computeCourantNumber(t_order,parent,Map)
if CFL==1.:
    #print "CFL=1 changed to CFL=0.99999999999"
    CFL=0.99999999999
Dt=CFL*dx/c 
tfinal=0.75*L/c
tf=2.*tfinal#0.75*L/c;
inc=round(tfinal/Dt)

update_position=False


T=0.
n=0

# Material points' fields
Md=mass*np.eye(Mp,Mp)
U = np.zeros((Mp,2))

U[:,1] = v0*np.sin(np.pi*xp[:,0]/L)/rho

# Nodes' fields
u = np.zeros((Nn,2))

# Storage
Stress=np.zeros((Mp,int(inc)+2))
Sth=np.zeros((Mp,int(inc)+2))
Vth=np.zeros((Mp,int(inc)+2))
Velocity=np.zeros((Mp,int(inc)+2))
Pos=np.zeros((Mp,int(inc)+2))
time=np.zeros(int(inc)+2)
energy=np.zeros(int(inc)+2)
Stress[:,0]=U[:,0]
Velocity[:,0]= U[:,1]
Vth[:,0]= U[:,1]
Pos[:,0]=xp[:,0]
time[0]=T

# Build approximation matrices
#Map,Grad,Dofs,parent=mesh.buildApproximation(xp)


mg=np.dot(np.dot(Map[Dofs,:],Md),Map[Dofs,:].T)
md=np.diag(np.sum(mg,axis=1))
mass_vector = np.dot(np.dot(Map,Md),Map.T)
mass_vector = np.sum(mass_vector,axis=1)
K=np.dot(np.dot(Grad[Dofs,:],Md),Map[Dofs,:].T)

dim=2 # Number of unknowns
mesh.setMapping(K)


def UpdateState(dt,dofs,M,U,md,limiter):
    if limiter!=-1:
        limit=True
    else:
        limit=False
    Nnodes = np.shape(U)[0]
    f=mesh.computeFlux(U,dofs,md,limit,limiter)
    """
    for i in range(len(dofs)):
        U[dofs[i],:]+=dt*f[dofs[i],:]/md[dofs[i]]
    """
    U[dofs,:]+=dt*np.linalg.solve(M,f[dofs])
    return U

def UpdateStateRK2(dt,dofs,M,U,md,limiter):
    Nnodes = np.shape(U)[0]
    k1=np.zeros((Nnodes,2))
    k2=np.zeros((Nnodes,2))
    if limiter!=-1:
        limit=True
    else:
        limit=False
    # first step : compute flux and intermediate state w
    f=mesh.computeFlux(U,dofs,md,limit,limiter)
    k1[dofs,:]+=dt*np.linalg.solve(M,f[dofs])
    w = U+0.5*k1
    # second step : compute flux and update U
    f=mesh.computeFlux(w,dofs,md,limit,limiter)
    k2[dofs,:]+=dt*np.linalg.solve(M,f[dofs])
    U+=k2
    return U



# Effective mass matrix coefficient
if ppc !=1:
    alpha=1.
else:
    alpha=1.
    """
limiter =-1 : none
          0 : minmod
          1 : superbee
          2 : muscl
"""

limiter=-1

#print '... computing ...'
while T<tfinal:
    
    # Effective mass matrix
    mf=(1-alpha)*mg + alpha*md
    
    # Mapping from material points to nodes
    
    Um=np.dot(Md,U)
    """
    for i in range(len(Dofs)):
        if mass_vector[Dofs[i]]!=0.:
            u[Dofs[i],:]=np.dot(Map[Dofs[i],:],Um)/mass_vector[Dofs[i]]
    """
   
    u[Dofs,:]=np.linalg.solve(mf,np.dot(Map[Dofs,:],Um))
    # Apply load on first node
    u[2*parent[0],0]= u[1,0] 
    u[2*parent[0],1]=-u[1,1]
    # Transmissive boundary conditions
    u[2*parent[-1]+3,0] = u[2*parent[-1]+2,0]
    u[-1,1] = -u[-2,1]

    
    if t_order==1 :
        u=UpdateState(Dt,Dofs,mf,u,mass_vector,limiter)
    elif t_order==2 :
        u=UpdateStateRK2(Dt,Dofs,mf,u,mass_vector,limiter)
    
    # Mapping back to the material points
    U=np.dot(Map[Dofs,:].T,u[Dofs,:])
    
    if update_position :
        xp[:,0]+=Dt*U[:,1]
        # Compute new mapping (convective phase)
        Map,Grad,Dofs,parent=mesh.buildApproximation(np.asmatrix(xp))
        mg=np.dot(np.dot(Map[Dofs,:],Md),Map[Dofs,:].T)
        md=np.diag(np.sum(mg,axis=1))
        mass_vector = np.dot(np.dot(Map,Md),Map.T)
        mass_vector = np.sum(mass_vector,axis=1)
        K=np.dot(np.dot(Grad[Dofs,:],Md),Map[Dofs,:].T)
        u = np.zeros((Nn,2))
        mesh.setMapping(K)


    #print 'Increment =', n, 't = ', T,' s.'
    n+=1
    T+=Dt
    Velocity[:,n]=U[:,1]
    Stress[:,n]=rho*U[:,0]
    
    for i in range(Mp):
        Sth[i,n]=E*v0*np.sin(np.pi*c*T/L)*np.cos(np.pi*xp[i,0]/L)/c
        Vth[i,n]=v0*np.cos(np.pi*c*T/L)*np.sin(np.pi*xp[i,0]/L)
    Pos[:,n]=xp[:,0]
    time[n]=T
    """
    plt.plot(xp[:,0],Stress[:,n],'r-o',lw=2.5)
    plt.plot(xp[:,0],Sth[:,n],'k-',lw=1.5)
    plt.grid()
    plt.show()
    """
Increments=n
"""
fig, (ax2) = plt.subplots(1,1)
plt.title('Elastic wave in a bar')

line1, = ax2.plot([], [],'b-s' , lw=2.,ms=1.5)
line2, = ax2.plot([], [],'r-s' , lw=2.,ms=1.5)
line3, = ax2.plot([], [],'k--' , lw=1.,ms=1.5)
line4, = ax2.plot([], [],'k--' , lw=1.,ms=1.5)
line = [line1,line2, line3,line4]

ax2.set_xlim(0.,L);
ax2.set_ylim(np.min(Stress),np.max(Stress));

ax2.grid()
ax2.set_xlabel('x (m)', fontsize=18)
ax2.set_ylabel('Stress', fontsize=18)
fig.legend((line),("V","S"),'upper right',numpoints=1)

def init():
    line[0].set_data([], [])
    line[1].set_data([], [])
    return line

def animate(i):
    line[0].set_data(Pos[:,i],Velocity[:,i]*rho*c)
    line[1].set_data(Pos[:,i],Stress[:,i])
    line[2].set_data(Pos[:,i],Vth[:,i]*rho*c)
    line[3].set_data(Pos[:,i],Sth[:,i])
    return line

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=Increments, interval=50, blit=True)
plt.show()
"""
