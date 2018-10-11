#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import *
from pylab import *
import scipy.optimize as optimize
import pdb
L=length
Mp=Nelem*ppc
meshFactor=3
Nelem*=meshFactor
Nn=Nelem +1       # Number of nodes
lx=L*ppc/(Mp-1)    # Length of cells
############ METHODES
def buildMesh(Nelem,Mp,solid_length,ppc,meshFactor):
    dx_mp=solid_length/(Mp-1)
    mesh_length=meshFactor*(solid_length + dx_mp)
    xn = np.linspace(0.,mesh_length,Nelem+1)-dx_mp/2.
    connect = np.array([np.arange(0,Nelem,1),np.arange(1,Nelem+1,1)]).T
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
    
    return Map,Grad,d,Parent.astype(int)-1

def celerity(E,rho0,J):
    return np.sqrt(E*(3.*J**2 -1.)/(2.*rho0))

################## END OF METHODES


xn,connect=buildMesh(Nelem,Mp,L,ppc,meshFactor)
xp=np.zeros((Mp,2))
xp[:,0]=np.linspace(0.,L,Mp)

"""
plt.plot(xp[:,0],xp[:,1],'ro')
plt.plot(xn,np.zeros(len(xn)),'b+')
plt.grid()
plt.show()
"""
# Material properties
Sy  = Sigy
c=np.sqrt(C/rho)

m=rho*lx/ppc

# Define force
s0 = sigd

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
Map,Grad,Dofs,parent=buildApproximation1D(xp,xn,connect)

mg=np.dot(np.dot(Map,Md),Map.T)
md=np.diag(np.sum(mg,axis=1))
mv=np.sum(np.dot(np.dot(Map,Md),Map.T),axis=1)

alpha=1e0 # for effective mass matrix

# def computeTimeStep(J,dx,c0,CFL):
#     c=np.sqrt(J)*c0
#     dt = CFL*dx/c
#     return dt
def computeTimeStep(rho,Tangent,J,dx,CFL):
    c=np.sqrt(Tangent*(3.*J**2 -1.)/(2.*rho))
    dt = CFL*dx/c
    return dt

J=np.ones(Mp)
Def[:,0]=np.copy(J)

algo = 'USF'
density=rho*np.ones(Mp)
#print "DGMPM time step ",Dt*factor,C
for n in range(NTMaxi)[1:]:    
    
    # Convection
    v[Dofs]=np.dot(Map[Dofs,:],np.dot(Md,V))/mv[Dofs]
    #v[connect[parent[-1],1]]=0.

    
    # Jmax=np.max(J)
    # Dt=computeTimeStep(Jmax,lx,c,CFL)
    Jmax=np.max(Def[:,n-1])
    point=np.where(Def[:,n-1]==Jmax)[0]
    if len(point)>1:
        point=point[0]
    #if updated_lagrangian: c =  celerity(C,density[point],f0)
    Dt=computeTimeStep(density[point],C,Jmax,lx,CFL)
    
    if ((time[n-1]+Dt)>tfinal):
        Dt = tfinal - time[n-1]
    time[n]=time[n-1]+Dt
    
    # Forces building
    Fe[connect[parent[0],0]]=-s0*(time[n-1]<tf)
    Fi[Dofs]=-np.dot(Grad[Dofs,:],np.dot(Md,Sig))
    
    # Solve motion equation
    a[Dofs]=(Fe[Dofs]+Fi[Dofs])/mv[Dofs]
    #a[connect[parent[-1],1]]=0.
    v+=a*Dt
    #v[connect[parent[-1],1]]=0.
    
    # Gradient and constitutive model
    ## Velocity gradient
    L=np.dot(Grad[Dofs,:].T,v[Dofs])
    dF=Dt*L
    if updated_lagrangian:
        Def[:,n]=(1.+dF)*Def[:,n-1]
    else:
        Def[:,n]=Def[:,n-1]+dF
        
    for i in range (Mp):
        Sig[i]=C*(Def[i,n]**3 -Def[i,n])/(2.*density[i]) 
        
    
    
    # Lagrangian step
    A=np.dot(Map[Dofs,:].T,a[Dofs])
    V+=Dt*np.dot(Map[Dofs,:].T,a[Dofs])
    
    xp[:,0]+=np.dot(Map[Dofs,:].T,v[Dofs])*Dt
    if updated_lagrangian:
        # Compute new mapping (convective phase)
        prevParent=parent
        Map,Grad,Dofs,parent=buildApproximation1D(np.asmatrix(xp),xn,connect)
        if (parent!=prevParent).any():
            print "********************************"
            print "*    GRID CROSSING OCCURING    *"
            print "********************************"
            
        mg=np.dot(np.dot(Map[Dofs,:],Md),Map[Dofs,:].T)
        md=np.diag(np.sum(mg,axis=1))
    
    a=np.zeros(Nn)
    v=np.zeros(Nn)

    Velocity[:,n]=V
    Pi[:,n]=density*Sig
    Pos[:,n]=np.copy(xp[:,0])
    #Def[:,n]=np.copy(J)
    if updated_lagrangian :
        density= rho/Def[:,n]
        
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
"""
#Sigma
fig = plt.figure()
plt.grid()
ax = plt.axes(xlim=(0.,L), ylim=(np.min(Pi),1.5*np.max(Pi)))
ax.plot(xn,np.zeros(len(xn)),'b+', lw='2.', ms='8.')
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
    lineList[1].set_data(Pos[:,i],Pi[:,i])
    time_text.set_text('Stress ')
    return tuple(lineList)+(time_text,)

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=increments, interval=20, blit=True)
#Animation of the stress
plt.grid()
#anim.save('StressBar.mp4', extra_args=['-vcodec', 'libx264'])
plt.show()
"""
