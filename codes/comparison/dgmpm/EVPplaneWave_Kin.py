#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import *
from pylab import *
import scipy.integrate
from scipy import optimize as opt
#from scipy.integrate import ode
from scipy import optimize
from sympy import *
sys.path.append(os.path.expandvars('dgmpm'))
from class_elasticity import *


##########################

def positive_part(x):
    return 0.5*(x+np.abs(x))

def creep_law(f,Seq,eta,n):
    return ((positive_part(f)/eta)**n)*np.sign(Seq)

def criterion(Seq,Sigy):
    return (np.abs(Seq)-Sigy)

def computeSeq(S,EP,lam,mu,H):
    KK = 3.0*(H/2.0) +(mu*(3.0*lam+2.0*mu))/(lam+2.0*mu)
    #KK = H +0.5*(4.*mu*(lam+mu)-lam**2)/(lam+2.0*mu)
    return (((2.0*mu)/(lam+2.0*mu))*S-KK*EP)


def computeS(t,U,EP,lam,mu,H,Sigy,eta,n):
    S = np.zeros(len(U))
    Seq = computeSeq(U[0],EP,lam,mu,H)
    f = criterion(Seq,Sigy)
    S[0] = -(2.0*mu)*creep_law(f,Seq,eta,n)
    S[1] = 0
    return S

def integrateODE(dt,U,EP,H,sigy,eta,n,mu,lam):
    # subdivision of time step
    sub_steps = 1
    sig_i = np.copy(U)
    s1 = np.copy(U)
    DT = dt/sub_steps
    bound = sigy + 0.5*3.*H*EP
    theta = 1.0
    for i in range(sub_steps):
        ## Nonlinear solution procedure
        R=lambda x: x -s1 + 2.0*mu*DT*(theta*creep_law(criterion(computeSeq(x,EP,lam,mu,H),sigy),computeSeq(x,EP,lam,mu,H),eta,n)+(1.0-theta)*creep_law(criterion(computeSeq(s1,EP,lam,mu,H),sigy),computeSeq(s1,EP,lam,mu,H),eta,n))        
        #sig = opt.brentq(R,bound,s1)
        sig = opt.fsolve(R,s1)
        s1 = sig
    return sig

def bar(x1,x2,Mp):
    xp=np.zeros((Mp,2))
    xp[:,0]=np.linspace(x1,x2,Mp)
    return xp

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
        sol.append(optimize.newton(Residual,1.))
    solution_CFL=np.min(sol)
    print "Courant number set to :",solution_CFL
    return solution_CFL

##########################


# Define geometry of the problem
L=length               # Length of the bar
Mp=Nelem*ppc             # Number of Material points
Nn=Nelem*2 + 2             # Number of elements


# Material properties
E=Young
Sy=Sigy           
c=np.sqrt(E/rho)

lam= (E*nu)/((1.0+nu)*(1.0-2.0*nu)) 
mu = E/(2.0*(1.0+nu))
c = np.sqrt((lam+2.0*mu)/rho)

power=n
# relaxation time
tau=(eta/Sigy)**n

mesh = DGmesh(Mp,L,ppc,c,rho)
dx=mesh.xn[1]-mesh.xn[0]
xp=bar(0.,L,Mp)

# Build approximation matrices
Map,Grad,Dofs,parent=mesh.buildApproximation(xp)

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


T=0.
n=0

# Material points' fields
Md=mass*np.eye(Mp,Mp)
U = np.zeros((Mp,2))

#v0=Sy/(2*rho*c)
U[0:Mp/2,1]=v0
U[Mp/2:Mp,1]=-v0


# Nodes' fields
u = np.zeros((Nn,2))

# Storage
sig=np.zeros((Mp,int(inc)+2))
epsp = np.zeros((Mp,int(inc)+2))
p = np.zeros((Mp,int(inc)+2))
velo=np.zeros((Mp,int(inc)+2))
pos=np.zeros((Mp,int(inc)+2))
time=np.zeros(int(inc)+2)
NRG=np.zeros(int(inc)+2)
kin=np.zeros(int(inc)+2)
strain=np.zeros(int(inc)+2)
sig[:,0]=U[:,0]
velo[:,0]= U[:,1]
pos[:,0]=np.copy(xp[:,0])
time[0]=T


mg=np.dot(np.dot(Map[Dofs,:],Md),Map[Dofs,:].T)
md=np.diag(np.sum(mg,axis=1))
mass_vector = np.dot(np.dot(Map,Md),Map.T)
mass_vector = np.sum(mass_vector,axis=1)
K=np.dot(np.dot(Grad[Dofs,:],Md),Map[Dofs,:].T)
alpha=1.e0 # for effective mass matrix

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

def UpdateState(dt,dofs,Ml,U,md,limit):
    Nnodes = np.shape(U)[0]
    if limit!=-1 : boolean=True
    else : boolean=False
    f=mesh.computeFlux(U,dofs,md)
    dU=np.zeros(np.shape(U))
    for i in range(len(dofs)):
        if md[dofs[i]]!=0.:
            dU[dofs[i],:]=dt*f[dofs[i],:]/md[dofs[i]]
            U[dofs[i],:]+=dt*f[dofs[i],:]/md[dofs[i]]
    return U,dU

def UpdateStateRK2(dt,dofs,Ml,U,md):
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


n=0
kin[n]=0.5*np.inner(np.dot(Md,velo[:,n]),velo[:,n])
strain[n]=0.5*np.inner(np.dot(Md,sig[:,n]/rho),sig[:,n]/E)
NRG[n]=kin[n]+strain[n]



while T<tfinal:
    
    # Effective mass matrix
    mf=(1-alpha)*mg + alpha*md


    # ## integration of ODE on material points
    # for j in range(Mp):
    #     if (Dt/tau)<10.:
    #         ## Non stiff problem
    #         r = scipy.integrate.ode(computeS).set_integrator('vode', method='adams',order=5)
    #     else:
    #         ## Stiff problem
    #         r = scipy.integrate.ode(computeS).set_integrator('vode', method='bdf',order=12)
    #     r.set_initial_value(rho*U[j,:],time[n]+Dt).set_f_params(epsp[j,n],H,E,Sy,eta,power,tau)
    #     r.integrate(r.t+Dt/2.)
    #     if r.successful():
    #         U[j,:] = r.y/rho
    for j in range(Mp):
        U[j,0] = integrateODE(Dt/2.,rho*U[j,0],epsp[j,n],H,Sigy,eta,power,mu,lam)/rho
        
    # Mapping from material points to nodes
    Um=np.dot(Md,U)
    for i in range(len(Dofs)):
        if mass_vector[Dofs[i]]!=0.:
            u[Dofs[i],:]=np.dot(Map[Dofs[i],:],Um)/mass_vector[Dofs[i]]
            
    # Apply load on first node
    u[2*parent[0],0]=2.*s0*(T<tf) - u[2*parent[0]+1,0] 
    u[2*parent[0],1]=u[2*parent[0]+1,1]
    # Transmissive boundary conditions
    u[2*parent[-1]+3,0] =- u[2*parent[-1]+2,0]
    u[2*parent[-1]+3,1] =u[2*parent[-1]+2,1]


    if t_order==1:
        u,du=UpdateState(Dt,Dofs,md,u,mass_vector,limit)
    elif t_order==2:
        u=UpdateStateRK2(Dt,Dofs,md,u,mass_vector)
    
    # Mapping back to the material points
    U=np.dot(Map[Dofs,:].T,u[Dofs,:])

    
    #  integration of ODE on material points
    # for j in range(Mp):
    #     if (Dt/tau)<1.e2:
    #         ## Non stiff problem
    #         r = scipy.integrate.ode(computeS).set_integrator('vode', method='adams',order=5)
    #     else:
    #         ## Stiff problem
    #         r = scipy.integrate.ode(computeS).set_integrator('vode', method='bdf',order=12)
    #     #r = scipy.integrate.ode(computeS).set_integrator('vode', method='bdf')
    #     r.set_initial_value(rho*U[j,:],time[n]+Dt).set_f_params(epsp[j,n],H,E,Sy,eta,power,tau)
    #     r.integrate(r.t+Dt)
    #     if r.successful():
    #         U[j,:] = r.y/rho
    for j in range(Mp):
        U[j,0] = integrateODE(Dt,rho*U[j,0],epsp[j,n],H,Sigy,eta,power,mu,lam)/rho

            
    if update_position:
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
    
    n+=1
    T+=Dt
    velo[:,n]=U[:,1]
    sig[:,n]=rho*U[:,0]

    Seq=computeSeq(sig[:,n],epsp[:,n-1],lam,mu,H)
    f = criterion(Seq,Sy)
    epsp[:,n] = epsp[:,n-1] +  creep_law(f,Seq,eta,power)*Dt
    
        
    pos[:,n]=xp[:,0]
    time[n]=T

    kin[n]=0.5*np.inner(np.dot(Md,velo[:,n]),velo[:,n])
    strain[n]=0.5*np.inner(np.dot(Md,sig[:,n]/rho),sig[:,n]/E)
    NRG[n]=kin[n]+strain[n]

    increments=n

x=mesh.xn


