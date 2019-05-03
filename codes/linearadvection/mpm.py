#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import *
from pylab import *
import math

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
        parent=math.ceil((xp[Pt,0]-xn[0])/Le)-1 
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

def bar(x1,x2,Mp):
    xp=np.zeros((Mp,2))
    xp[:,0]=np.linspace(x1,x2,Mp)
    return xp


print 'Initializing problem ...'
# Define geometry of the problem
L=1.               # Length of the bar
Mp=60             # Number of Material points
ppc=1
Ne=Mp/ppc            # Number of elements
Nn=Ne+1              # NUmber of nodes 
Le=L*ppc/(Mp-1)    # Length of cells
lmp = L/(Mp-1)

print '       Mesh Definition'
#xn,connect=buildMesh(Nn,L+Le)
xn =np.linspace(-lmp/2,L+lmp/2,Nn)
connect = np.array([np.arange(0,Nn-1,1),np.arange(1,Nn,1)]).T

xp=bar(0.,L,Mp)
#xp[:,0]+=0.002

plt.plot(xp[:,0],xp[:,1],'ro',label='Material points')
plt.plot(xn,xn*0.,'b+',label='Nodes')
plt.axis('equal')
plt.legend(loc='best',numpoints=1)
plt.show()


# Material properties
rho=7800.
E=2.e11
Sy  = 400.0e6           #Traction yield stress
c=math.sqrt(E/rho)

mass=rho*Le/ppc

# Define force
Q0=0.5*Sy/rho

print '       Algorithmic parameters'
# Time discretization
CFL = 0.1#/ppc
Dt=CFL*Le/c 
tfinal=0.75*L/c
inc=round(tfinal/Dt)

T=0.
n=0
Q=np.zeros(Mp)


X=np.zeros(Mp)
X[:]=xp[:,0]
Map=np.zeros((Nn,Mp))
Grad=np.zeros((Nn,Mp))
Md=mass*np.eye(Mp,Mp)

q=np.zeros(Nn)
Fint=np.zeros(Nn)

# Storage
Stress=np.zeros((Mp,int(inc)+2))
time=np.zeros(int(inc)+2)
Stress[:,0]=Q*rho
time[0]=T
# Build approximation matrices
Map,Grad,Dofs=buildApproximation1D(xp,xn,connect)

mg=np.dot(np.dot(Map,Md),Map.T)
md=np.diag(np.sum(mg,axis=1))
K=np.dot(Grad,np.dot(Md,Map.T))
print K
def updateStateRK2(Dt,xp,xn,connect,rho,Md,V,Sig,f):
    a=np.zeros(Nn)
    v=np.zeros(Nn)
    Fe=np.zeros(Nn)
    Fi=np.zeros(Nn)
    
    Map,Grad,Dofs=buildApproximation1D(np.asmatrix(xp),xn,connect)
    mg=np.dot(np.dot(Map[Dofs,:],Md),Map[Dofs,:].T)
    md=np.diag(np.sum(mg,axis=1))
    
    # Convection
    v[Dofs]=solve(md,np.dot(Map[Dofs,:],np.dot(Md,V)))
    k1=np.copy(v)

    # Forces building
    Fe[0]=f
    Fi[Dofs]=-np.dot(Grad[Dofs,:],np.dot(Md,Sig))/rho
    
    # Solve motion equation
    a[Dofs]=solve(md,Fe[Dofs]+Fi[Dofs])
    k1+=a*Dt/2.
        
    k1[-1]=0.
    # Gradient and constitutive model
    dEps=Dt*np.dot(Grad[Dofs,:].T,k1[Dofs])
    Sig+=E*dEps
        
    k2=np.copy(k1)
    Fi=np.zeros(Nn)
    Fi[Dofs]=-np.dot(Grad[Dofs,:],np.dot(Md,Sig))/rho
    a[Dofs]=solve(md,Fe[Dofs]+Fi[Dofs])
    k2+=a*Dt
    
    v+=k2
    # Lagrangian step
    A=np.dot(Map[Dofs,:].T,a[Dofs])
    V+=Dt*A/2.
    # Gradient and constitutive model
    dEps=Dt*np.dot(Grad[Dofs,:].T,k1[Dofs])
    Sig+=E*dEps
    
    return V,Sig

alpha=1.e0 # for effective mass matrix
alg='USL' # USL or USF (Update Stress Last or First)
Q[0]=Q0

print '========================================================'
print '=                  COMPUTATION STARTED                 ='
print '========================================================'
while T<tfinal:
    
    # Effective mass matrix
    mf=(1-alpha)*mg + alpha*md
    
    # Convection
    q=solve(mf,np.dot(Map,np.dot(Md,Q)))
    
    # Forces building
    Fint=np.dot(K,c*q)
    
    # Solve advection equation on the mesh
    dq=solve(mf,Fint)
    q+=Dt*dq
        
    # Map back to material points
    Q+=Dt*np.dot(Map.T,dq)
    #Q=np.dot(Map.T,q)
    Q[0]=Q0

    q=np.zeros(Nn)
    
    print 'Increment =', n, 't = ', T,' s.'
    n+=1
    T+=Dt
    Stress[:,n]=Q[:]*rho
    time[n]=T
    
print '========================================================'
print '=                 COMPUTATION COMPLETED                ='
print '========================================================'
print '... waiting for animation building'


frames = [0,1,2,3,4,5,6,7,8,9,10]
"""
for i in frames:
    plt.plot(xp[:,0],Stress[:,i],'ro',label='t =  '+str(time[i]))
    plt.grid()
    plt.show()
"""


print 'Animation finished !'
####Animated plot ###########################################################
# First set up the figure, the axis, and the plot element we want to animate
pas=int(1)

fig = plt.figure()
plt.grid()
#ax = plt.axes(xlim=(xp[0],xp[-1]), ylim=(-1.5))
ax = plt.axes(xlim=(Le/2,L+Le/2), ylim=(np.min(Stress),np.max(Stress)))
ax.grid()
ax.set_xlabel('x (m)', fontsize=18)
ax.set_ylabel(r'$\sigma$ (Pa)', fontsize=18)
ax.set_title('Stress wave propagation in a bar (Courant Number = '+str(0.5)+')', fontsize=16)
line, = ax.plot([], [],'r', lw=2.)
line2,= ax.plot([], [],'b--', lw=1.5)
fig.legend((line,line2),('MPM','Analytical'),'upper right',numpoints=1)
time_text = ax.text(0.02, 0.95, 'middle', transform=ax.transAxes)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line,time_text

# animation function.  This is called sequentially
def animate(i):
    line.set_data(xp[:,0],Stress[:,i])
    #time_text.set_text('Stress (Pa) at time = '+str(time[i]))
    return line,time_text
 
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=Stress.shape[1], interval=50, blit=True)

#Animation of the stress
plt.show()
#anim.save('StressBar.mp4', extra_args=['-vcodec', 'libx264'])


