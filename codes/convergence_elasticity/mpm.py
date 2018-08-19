# !/usr/bin/python

from pylab import *
from scipy.linalg import solve
import numpy as np
from matplotlib import animation
import pdb

L=1.
Nelem=Mp
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

rho=1.e3
E=1.e7
Sy=400.0e6           
c=np.sqrt(E/rho)
H = 10e9
HT=E*H/(E+H)
cp=np.sqrt(HT/rho)
m = rho*lx/ppc


s0 =-0.*Sy
v0=0.1

# Time discretization
CFL=0.4
Dt=CFL*lx/c
tfinal=0.5*L/c

tfinal=0.02
Dt=1.e-3
tf=2.*tfinal#0.75*L/c

inc=int(tfinal/Dt)+1

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

G=1.e-4
V = c*np.pi*G*np.sin(np.pi*xp[:,0])

# Grid nodes' fields
a=np.zeros(Nn)
v=np.zeros(Nn)
Fe=np.zeros(Nn)
Fi=np.zeros(Nn)

# Approximation builing
Map,Grad,Dofs=buildApproximation1D(np.asmatrix(xp),xn,connect)
MD=m*np.eye(Mp)

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
Stress[:,n]=Sig[:]
Pos[:,n]=X[:]
Velocity[:,n]=V[:]
time[n]=T

Sig,p[:,0],Epsp[:,0],tangent_modulus = bilinear(np.zeros(len(Eps)),Eps,Epsp[:,0],p[:,0],E,Sy,H)

alpha=1.

def loadingFunction(t,tfinal,tf):
    return np.sin(period*np.pi*t/tfinal)*(t<tf)

def smoothSolution(x,c,t,tfinal,tf):
    if c*t < x:
        val=0.
    else:
        val=loadingFunction((t-x/c),tfinal,tf)
    return val

#alg='USF'
for n in range(int(inc)+2)[1:]:
    mf=(1-alpha)*mn + alpha*md

    T+=Dt
    time[n]=T
    
    if time[n]>tfinal:
        Increments=n-1
        break
    if alg=='USL':
        
        # Convection /!\
        v[Dofs]=np.dot(Map[Dofs,:],np.dot(MD,V))/mv[Dofs]
        #v[Dofs]=solve(mf[Dofs],np.dot(Map[Dofs,:],np.dot(MD,V)))
        v[0]=0.
        v[-1]=0.
        
        # Internal force vector
        Fe[0]=0.#-s0*loadingFunction(T+Dt,tfinal,tf)
        Fi[Dofs]=-np.dot(Grad[Dofs,:],np.dot(MD,Sig))/rho

        # Solve motion equation
        a[Dofs]=(Fe[Dofs]+Fi[Dofs])/mv[Dofs]
        #a[Dofs]=solve(mf[Dofs],Fe[Dofs]+Fi[Dofs])
        
        
        # Time integrator
        a[0]=0.
        a[-1]=0.
        v+=Dt*a
           
        
        A=np.dot(Map[Dofs,:].T,a[Dofs])
        V+=Dt*A

        # Gradient and constitutive model
        Epsn=Eps
        
        Eps+=Dt*np.dot(Grad[Dofs,:].T,v[Dofs])
        Sig,p[:,n],Epsp[:,n],tangent_modulus = bilinear(Epsn,Eps,Epsp[:,n],p[:,n],E,Sy,H)
        
    elif alg=='USF':
        # Convection /!\
        v[Dofs]=np.dot(Map[Dofs,:],np.dot(MD,V))/mv[Dofs]
        #v[Dofs]=solve(mf,np.dot(Map[Dofs,:],np.dot(MD,V)))
        v[0]=0.
        v[-1]=0.
        
        
        # Gradient and constitutive model
        Epsn=Eps
        Eps+=Dt*np.dot(Grad[Dofs,:].T,v[Dofs])
        Sig,p[:,n],Epsp[:,n],tangent_modulus = bilinear(Epsn,Eps,Epsp[:,n],p[:,n],E,Sy,H)
        
        # Internal force vector
        Fi[Dofs]=-np.dot(Grad[Dofs,:],np.dot(MD,Sig))/rho
        Fe[0]=-s0*loadingFunction(T+Dt,tfinal,tf)
        
        # Solve motion equation
        a[Dofs]=(Fe[Dofs]+Fi[Dofs])/mv[Dofs]
        #a[Dofs]=solve(mf,Fi[Dofs]+Fe[Dofs])
        a[0]=0.
        a[-1]=0.
        
        # Time integrator
        v+=Dt*a
        A=np.dot(Map[Dofs,:].T,a[Dofs])
        V+=Dt*A
        
    a=np.zeros(Nn)
    v=np.zeros(Nn)
    Fi=np.zeros(Nn)

    
    
    Stress[:,n]=Sig[:]
    Pos[:,n]=xp[:,0]
    Velocity[:,n]=V[:]
    
    for i in range(Mp):
        Sth[i,n]=E*np.pi*G*np.sin(np.pi*c*T)*np.cos(np.pi*xp[i,0])
        #Vth[i,n]=v0*np.cos(np.pi*c*T/L)*np.sin(np.pi*xp[i,0]/L)
        Vth[i,n]=c*np.pi*G*np.cos(np.pi*c*T)*np.sin(np.pi*xp[i,0])
        
    # if T+Dt>tfinal:
    #     break
"""
fig, (ax2) = plt.subplots(1,1)
plt.title('Elastic wave in a bar')
line1, = ax2.plot([], [],'b-s' , lw=2.,ms=1.5)
line2, = ax2.plot([], [],'r-s' , lw=2.,ms=1.5)
line = [line1,line2]
ax2.set_xlim(0.,L);
ax2.set_ylim(np.min(Velocity),np.max(Velocity));
 
ax2.grid()
ax2.set_xlabel('x (m)', fontsize=18)
ax2.set_ylabel('Stress', fontsize=18)
fig.legend((line),("V","S"),'upper right',numpoints=1)
def init():
    line[0].set_data([], [])
    line[1].set_data([], [])
    return line
def animate(i):
    line[0].set_data(Pos[:,i],Velocity[:,i])
    line[1].set_data(Pos[:,i],Vth[:,i])
    return line

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=Increments, interval=50, blit=True)
plt.show()
"""
