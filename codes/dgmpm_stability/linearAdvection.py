#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import *
from pylab import *
from scipy import optimize
from DGmesh import *
import pdb
from sympy import *

##########################
def bar(x1,x2,Mp):
    xp=np.zeros((Mp,2))
    xp[:,0]=np.linspace(x1,x2,Mp)
    return xp

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

def UpdateState(dt,dofs,Ml,U,md,dim,limiter):
    Nnodes = np.shape(U)[0]
    f=mesh.computeFlux(U,dofs,md,dim,limiter)
    U[dofs]+=dt*(np.linalg.solve(Ml,f[dofs]))
    return U

def UpdateStateRK2(dt,dofs,Ml,U,md,dim,limiter):
    Nnodes = np.shape(U)[0]
    k1=np.zeros(Nnodes)
    k2=np.zeros(Nnodes)
    # first step : compute flux and intermediate state w
    f=mesh.computeFlux(U,dofs,md,dim,limiter)
    k1[dofs]+=dt*(np.linalg.solve(Ml,f[dofs]))
    w = U+k1*0.5
    # second step : compute flux and update U
    f=mesh.computeFlux(w,dofs,md,dim,limiter)
    k2[dofs]+=dt*(np.linalg.solve(Ml,f[dofs]))
    U+=k2
    return U

def buildDiscreteOperator(Map,parent,invParent,t_order):
    CFL = symbols('CFL')
    H_matrix=zeros(Map.shape[1])
    for p in range(np.shape(Map)[1]):
        # loop over the material points
        # 1) find the connected material points
        ## a- in the same cell
        elem = parent[p]
        if elem==0: continue
        sharing=invParent[elem]
        ## b- in the previous cell
        previous=invParent[elem-1]
        # 2) build the Hij terms
        # index of the left and right nodes
        n1 = 2*elem+1 ; n2 = n1+1
        for k in sharing:
            H_matrix[p,k] = Map[n1,p]*Map[n1,k]/np.sum(Map[n1,:]) + Map[n2,p]*Map[n2,k]/np.sum(Map[n2,:]) + CFL*(Map[n2,p]/np.sum(Map[n2,:]) - Map[n1,p]/np.sum(Map[n1,:]) - len(sharing)*Map[n2,p]*Map[n2,k]/(np.sum(Map[n2,:])**2) )
            if t_order==2:
                H_matrix[p,k] += 0.5*len(sharing)*(CFL**2)*((Map[n2,k]/np.sum(Map[n2,:]))*(Map[n1,p]/np.sum(Map[n1,:])-Map[n2,p]/np.sum(Map[n2,:])) + (len(sharing)/np.sum(Map[n2,:]))*(Map[n2,k]/np.sum(Map[n2,:])-1.))
                
        # index of the left and right nodes in the previous element
        np1 = 2*(elem-1)+1 ; np2 = np1+1
        for k in previous:
            H_matrix[p,k] = CFL*len(sharing)*Map[n1,p]*Map[np2,k]/(np.sum(Map[n1,:])*np.sum(Map[np2,:]))
            if t_order==2:
                H_matrix[p,k] +=0.5*len(sharing)*(CFL**2)*( Map[n1,p]/(np.sum(Map[n1,:])*np.sum(Map[np2,:]))*(1.-len(sharing)*Map[np2,k]/np.sum(Map[np2,:])) -(Map[np2,k]/np.sum(Map[np2,:]))*(Map[n1,p]/np.sum(Map[n1,:])-Map[n2,p]/np.sum(Map[n2,:])) )
    # Deal with the first element
    for p in invParent[0]:
        sharing=invParent[0]
        # index of the left and right nodes
        n1 = 1 ; n2 = 2
        for k in sharing:
            H_matrix[p,k] = Map[n1,p]*Map[n1,k]/np.sum(Map[n1,:]) + Map[n2,p]*Map[n2,k]/np.sum(Map[n2,:]) + CFL*(Map[n2,p]/np.sum(Map[n2,:]) - Map[n1,p]/np.sum(Map[n1,:]) - len(sharing)*Map[n2,p]*Map[n2,k]/(np.sum(Map[n2,:])**2) )
            if t_order==2:
                H_matrix[p,k] += 0.5*len(sharing)*(CFL**2)*((Map[n2,k]/np.sum(Map[n2,:]))*(Map[n1,p]/np.sum(Map[n1,:])-Map[n2,p]/np.sum(Map[n2,:])) + (len(sharing)/np.sum(Map[n2,:]))*(Map[n2,k]/np.sum(Map[n2,:])-1.))#+ (Map[n2,p]/np.sum(Map[n2,:]))*(len(sharing)*Map[n2,k]/np.sum(Map[n2,:])-1.)/np.sum(Map[n2,:]))
                
    Hoperator = lambdify((CFL),H_matrix)
    return H_matrix,Hoperator

def UpdateStateDiscreteOperator(U,H,BC,CFL,invParent,t_order):
    U_updated=np.zeros(np.shape(U))
    for p in range(np.shape(U)[0]):
        for k in range(np.shape(U)[0]):
            U_updated[p]+=H[p,k]*U[k]
    ## next, enforce the BC at the left points
    for p in invParent[0]:
        sharing=invParent[0]
        # index of the left and right nodes
        n1 = 1 ; n2 = 2
        for k in sharing:
            Hpk = CFL*len(sharing)*Map[n1,p]*Map[n2,k]/(np.sum(Map[n1,:])*np.sum(Map[n2,:]))
            if t_order==2:
                Hpk += 0.5*len(sharing)*(CFL**2)*( Map[n1,p]/(np.sum(Map[n1,:])*np.sum(Map[n2,:]))*(1.-len(sharing)*Map[n2,k]/np.sum(Map[n2,:])) -(Map[n2,k]/np.sum(Map[n2,:]))*(Map[n1,p]/np.sum(Map[n1,:])-Map[n2,p]/np.sum(Map[n2,:])) )
                # Hpk += 0.5*len(sharing)*(CFL**2)*((Map[n2,k]/np.sum(Map[n2,:]))*(Map[n1,p]/np.sum(Map[n1,:])-Map[n2,p]/np.sum(Map[n2,:])) + (Map[n2,p]/np.sum(Map[n2,:]))*(len(sharing)*Map[n2,k]/np.sum(Map[n2,:])-1.)/np.sum(Map[n2,:]))
            U_updated[p]+=Hpk*BC
    return U_updated

def gridSearch(function,tol=1.e-7):
    samples=100000
    # Find the bigest root of the residual by grid search algorithm
    CFL=np.linspace(0.,1.,samples)
    for i in CFL:
        if i==CFL[samples-1]: return i
        a0=function(i)
        if a0<tol:
            continue
        else:
            return i

def computeCriticalCFL(Mp,H,invParent):
    CFL=symbols('CFL')
    sol=[]
    for p in invParent[1]:
        res = 0
        for k in range(Mp):
            res+=np.abs(H[p,k])
        # solve the residual
        residual=lambdify((CFL),res-1.)
        #solution=gridSearch(residual)
        solution=optimize.root(residual,1.,method='hybr',options={'xtol':1.e-4}).x[0]
        if abs(residual(solution))>1.e-3: print "CAUTION: residual norm after solution is", abs(residual(solution))
        print "CFL solution for point ",p," ",solution
        sol.append(solution)
    Courant=min(sol)
    print "Critical Courant number set to ",Courant
    return Courant

def computeLpNorm(Unum,Uexact,dx,p):
    return (((dx*((np.abs(Unum-Uexact))**p)).sum())**(1.0/p))

def computeRelativeError(Unum,Uexact,dx,p):
    return (computeLpNorm(Unum,Uexact,dx,p)/computeLpNorm(np.zeros(len(Unum)),Uexact,dx,p))

##########################


print 'Initializing problem ...'
# Define geometry of the problem
L=1.               # Length of the bar
Mp=100            # Number of Material points
ppc=2
Nelem=Mp/ppc             # Number of elements

Nn=Nelem*2 + 2             # Number of elements



# Material properties

rho=7800.
E=2.e11
Sy=400.0e6           
c=np.sqrt(E/rho)


print '       Mesh Definition'

mesh = DGmesh(Mp,L,ppc,c,rho)
dx=mesh.xn[1]-mesh.xn[0]
xp=bar(0.,L,Mp)

coor=np.zeros(Nn)
for i in range(Nn):
    coor[i]=(i/2)*dx -dx/(2.*ppc) 

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


# Build approximation matrices
Map,Grad,Dofs,parent=mesh.buildApproximation(xp)
invParent=[]
for i in range(Nelem):
    invParent.append(np.where(parent==i)[0])


print '       Algorithmic parameters'
# Time discretization
# Build the discrete operator
t_order= 1
Hsym,HOperator=buildDiscreteOperator(Map,parent,invParent,t_order)
CFL=computeCriticalCFL(Mp,Hsym,invParent)
Dt=CFL*dx/c 
tfinal=1.*L/c
tf=2.0*tfinal;
inc=round(tfinal/Dt)
tunload=5000*Dt
T=0.
n=0

# Material points' fields
Md=mass*np.eye(Mp,Mp)
U = np.zeros(Mp)
Uh = np.zeros(Mp)


# Nodes' fields
u = np.zeros(Nn)

# Storage
Stress=np.zeros((Mp,int(inc)+2))
Stressh=np.zeros((Mp,int(inc)+2))
analytical=np.zeros((Mp,int(inc)+2))
time=np.zeros(int(inc)+2)
Stress[:,0]=U[:]
Stressh[:,0]=U[:]
time[0]=T


mg=np.dot(np.dot(Map[Dofs,:],Md),Map[Dofs,:].T)
md=np.diag(np.sum(mg,axis=1))
mass_vector = np.dot(np.dot(Map,Md),Map.T)
mass_vector = np.sum(mass_vector,axis=1)
K=np.dot(np.dot(Grad[Dofs,:],Md),Map[Dofs,:].T)
alpha=1.e0 # for effective mass matrix

dim=1 # Number of unknowns
mesh.setMapping(K)


alpha=1.
"""
limiter = 0 : minmod
          1 : superbee
          2 : muscl
"""
limiter=-1

print '... computing ...'


while T<tfinal:
    
    # Effective mass matrix
    mf=(1-alpha)*mg + alpha*md
    
    
    # Mapping from material points to nodes
    u[Dofs]=np.linalg.solve(mf,np.dot(Map[Dofs,:],np.dot(Md,U)))
    
    u[0]=r0*(T<tunload)
    
    # plt.plot(coor,u*rho,'b-+',lw=2.)
    # plt.plot(xp[:,0],U*rho,'r-o',lw=2.)
    # plt.grid()
    # plt.show()
    
    if t_order==1 :
        u=UpdateState(Dt,Dofs,md,u,mass_vector,dim,limiter)
    elif t_order==2 :
        u=UpdateStateRK2(Dt,Dofs,md,u,mass_vector,dim,limiter)
    
    Uh=UpdateStateDiscreteOperator(Uh,HOperator(CFL),r0,CFL,invParent,t_order)
        
    # Mapping back to the material points
     
    U=np.dot(Map.T,u)
    
    #print U
    #xp[:,0]+=Dt*U[:,1]
    # Compute new mapping (convective phase)
    """
    Map,Grad,Dofs,parent=mesh.buildApproximation(np.asmatrix(xp))
    
    mg=np.dot(np.dot(Map[Dofs,:],Md),Map[Dofs,:].T)
    md=np.diag(np.sum(mg,axis=1))
    K=np.dot(np.dot(Grad[Dofs,:],Md),Map[Dofs,:].T)
    u = np.zeros((Nn,2))
    mesh.setMapping(K)
    """

    
    n+=1
    T+=Dt
    Stress[:,n]=rho*U
    Stressh[:,n]=rho*Uh
    for i in range(Mp):
        analytical[i,n]=r0*(c*T>xp[i,0] and c*(T-tunload)<xp[i,0])*rho
    
    time[n]=T
    """
    plt.plot(Pos[:,n],Stress[:,n],'r-o',label='DGMPM',lw =2.)
    plt.legend()
    plt.grid()
    plt.show()
    """
## Compute the error between the two numerical solutions
error = computeRelativeError(Stress[:,n],Stressh[:,n],dx,2)
print "Error between the two numerical procedures: ",error
print '... building animation ...'
print 'Animation finished !'
####Animated plot ###########################################################
# First set up the figure, the axis, and the plot element we want to animate

pas=int(1/CFL)

fig = plt.figure()
plt.grid()
#ax = plt.axes(xlim=(xp[0],xp[-1]), ylim=(-1.5))
ax = plt.axes(xlim=(0,L), ylim=(np.min(Stress),1.1*np.max(Stress))) #np.max(Stress)
ax.grid()
ax.set_xlabel('x (m)', fontsize=18)
ax.set_ylabel(r'$\sigma$ (Pa)', fontsize=18)
ax.set_title('Stress wave propagation in a bar', fontsize=16)
line, = ax.plot([], [],'ro', lw=2.)
line2,= ax.plot([], [],'k--', lw=1.5)
line3,= ax.plot([], [],'y-*', lw=1.5)
fig.legend((line,line2),('DGMPM','Analytical'),'upper right',numpoints=1)
time_text = ax.text(0.02, 0.95, 'middle', transform=ax.transAxes)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    time_text.set_text('')
    return line,line2,time_text

# animation function.  This is called sequentially
def animate(i):
    line.set_data(xp[:,0],Stress[:,i])
    line2.set_data(xp[:,0],analytical[:,i])
    line3.set_data(xp[:,0],Stressh[:,i])
    #time_text.set_text('Stress (Pa) at time = '+str(time[i]))
    return line,line2,line3,time_text
 
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=Stress.shape[1], interval=50, blit=True)

#Animation of the stress
plt.show()
#anim.save('StressBar.mp4', extra_args=['-vcodec', 'libx264'])



