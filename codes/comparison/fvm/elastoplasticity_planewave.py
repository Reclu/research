#!/usr/bin/python

import numpy as np


"""
Implementation of high order shock capturing methods to solve the 1D elastoplasticity set of equations in dynamics with a kinematic hardening
"""
E_A=Young
rho_A=rho
Sigy_A=Sigy
H_A=H
nu_A=nu
#Define 1D mesh
DX=length/Nelem
x = np.linspace(0.,length,Nelem+1)
x = np.linspace(-DX/2.,length+DX/2.,Nelem+1)
dx = x[1]-x[0]
y=x+(dx/2.0)
#Centres of cells
y=y[:(len(y)-1)]
centroids=y

#Define parameters
courantNumber = 1.0
t = np.zeros(NTmaxi)
#Define material parameters
lam_A = (nu_A*E_A)/(((1.0+nu_A)*(1.0-2.0*nu_A)))
mu_A = E_A/(2.0*(1.0+nu_A))
c_A = np.sqrt((lam_A+2.0*mu_A)/rho_A)

HT_A = (lam_A+2.0*mu_A)-8*(mu_A**2/(3.0*H_A+6.0*mu_A))
HT_A2 = (mu_A*(3.0*lam_A+2.0*mu_A)+(3.0/2.0)*H_A*(2.0*mu_A+lam_A))/(3.0*(mu_A+H_A/2.0))
HT_A3 = (lam_A+2.0*mu_A+((2.0*mu_A*lam_A)/(mu_A+(3.0*H_A/2.0)))) \
     / (1.0+ ((2.0*mu_A)/(mu_A+(3.0*H_A/2.0))))

cp_A = np.sqrt(HT_A/rho_A)

HEL = ((lam_A+2.0*mu_A)/(2.0*mu_A))*Sigy_A

"""Limiters of waves:
mthlim(i)<0: no limiter
         =0: minmod limiter
         =1: superbee limiter
         =2: van Leer limiter
         =3: Monotonized centered
         =4: Beam-Warming
         =5: Lax-Wendroff
"""
limit = True
mthlim = np.zeros(4,dtype=int)
for idx,i in enumerate(mthlim):
    mthlim[idx] = fvmlimiter   #Choose of the Superbee limiter

#Definition of arrays
#2 Ghost cells at both ends of the bar are used to prescribe boundary conditions
v = np.zeros((len(y)+4,NTmaxi))
disp = np.zeros((len(y)+4,NTmaxi))
eps = np.zeros((len(y)+4,NTmaxi))
sig = np.zeros((len(y)+4,NTmaxi))
EP = np.zeros((len(y)+4,NTmaxi))
EPeq = np.zeros((len(y)+4,NTmaxi))
U = np.zeros((len(y)+4,2))
matC = np.zeros(len(y)+4)
matCp = np.zeros(len(y)+4)
rho = np.zeros(len(y)+4)
Sigy = np.zeros(len(y)+4)
H = np.zeros(len(y)+4)
lam = np.zeros(len(y)+4)
mu = np.zeros(len(y)+4)

#Initial conditions and definition properties
for i in range(len(y)+4):
    if (i<len(range(len(y)+4))/2):
        v[i,0] = v0
        matC[i] = c_A 
        matCp[i] = cp_A 
        rho[i] = rho_A
        Sigy[i] = Sigy_A
        H[i] = H_A
        lam[:] = lam_A
        mu[:] = mu_A
    else:
        v[i,0] = -v0
        matC[i] = c_A 
        matCp[i] = cp_A
        rho[i] = rho_A
        Sigy[i] = Sigy_A
        H[i] = H_A
        lam[:] = lam_A
        mu[:] = mu_A

U[:,0] = sig[:,0]
U[:,1] = v[:,0]

#Definition of functions
def computeAlpha(dU,rho1,rho2,c1,c2):
    alpha1 = (dU[0]+rho2*c2*dU[1])/(rho1*c1+rho2*c2)
    alpha2 = -(dU[0]-rho1*c1*dU[1])/(rho1*c1+rho2*c2)
    return [alpha1,alpha2]

def computeLimitedWaves(waves,s,mthlim,dtdx):
    cqxx = np.zeros((Nelem+3,U.shape[1]))
    #Loop on interfaces
    for i in range(Nelem+2)[1:]:
        #Apply limiter to waves 
        modWaves = limiter(s[i,:],waves[i,:,:],waves[i-1,:,:],waves[i+1,:,:],mthlim,i)
        #Second order corrections of normal fluxes
        for j in range(len(s[i,:])):
            cqxx[i,:] += np.abs(s[i,j])*(1.0-(dtdx)*np.abs(s[i,j]))*modWaves[:,j]
    return cqxx

def limiter(s,wave,wave_L,wave_R,mthlim,j):
    modWaves = np.zeros((2,4))
    for i in range(len(s)):
        if (mthlim[i]>=0) and (not (wave[:,i]==np.zeros(2)).all()):
            #Compute the L2-norm of wave and projections
            wnorm2 = np.dot(wave[:,i],wave[:,i])
            if (wnorm2 == 0.0):
                print 'wnorm2 = 0 pour l onde',i,'de l edge',j
                print 'wave[:,:,',j,'] = ',wave
                break
            wwR = np.dot(wave[:,i],wave_R[:,i])
            wwL = np.dot(wave[:,i],wave_L[:,i])
            #Upwind and compute wave strength
            if (s[i]>0.0):
                wlimiter = philim(wnorm2,wwL,mthlim[i])
            elif (s[i]<0.0):
                wlimiter = philim(wnorm2,wwR,mthlim[i])
        else:
            wlimiter =0.0
        #Apply the limiter
        modWaves[:,i] = wlimiter*wave[:,i]
    return modWaves

def philim(a,b,meth):
    theta = b/a
    if (meth==0):
        philim = np.max([0.0,np.min([1.0,theta])])  #minmod
    elif (meth==1):
        philim = np.max([0.0,np.min([1.0,2.0*theta]),np.min([2.0,theta])])  #superbee
    elif (meth==2):
        philim = (theta+np.abs(theta))/(1.0+np.abs(theta))   #van Leer
    elif (meth==3):
        c=(1.0+theta)/2.0
        philim = np.max([0.0,np.min([c,2.0,2.0*theta])])   #Monotonized Centered
    elif (meth==4): 
        philim = theta                                   #Beam-Warming
    elif (meth==5):       
        philim = 1.0                                   #Lax-Wendroff
    return philim


def computeSeq(S,EP,lam,mu,H):
    KK = 3.0*(H/2.0) +(mu*(3.0*lam+2.0*mu))/(lam+2.0*mu)
    return ((2.0*mu)/(lam+2.0*mu))*S-KK*EP

def EP_RiemannSolver(U,EP,rho,lam,mu,Sigy,H):
    #Allocation of fluctuations arrays
    amdq = np.zeros((Nelem+3,U.shape[1]))
    apdq = np.zeros((Nelem+3,U.shape[1]))
    waves = np.zeros((Nelem+3,U.shape[1],2*U.shape[1]))
    s = np.zeros((Nelem+3,2*U.shape[1]))
    matC = np.sqrt((lam+2.0*mu)/rho)
    HT = (lam+2.0*mu)-8*(mu**2/(3.0*H+6.0*mu)) ; matCp = np.sqrt(HT/rho)
    #Loop on interfaces
    for i in range(Nelem+3):
        UL = U[i,:]; SL = U[i,0]
        UR = U[i+1,:] ; SR = U[i+1,0]
        dU = UR-UL
        #Elastic prediction: elastic Riemann solver
        alpha = computeAlpha(dU,rho[i],rho[i+1],matC[i],matC[i+1])
        Strial = SL+alpha[0]*rho[i]*matC[i]
        #Tests on the criterion
        SS_star_L = computeSeq(Strial,EP[i],lam[i],mu[i],H[i])
        SS_star_R = computeSeq(Strial,EP[i+1],lam[i+1],mu[i+1],H[i+1])
        fL = np.abs(SS_star_L)-Sigy[i]
        fR = np.abs(SS_star_R)-Sigy[i+1]
        if (fL>0.0) and (fR<0.0):
            #Leftward plastic wave
            KK_L = 3.0*(H[i]/2.0) +(mu[i]*(3.0*lam[i]+2.0*mu[i]))/(lam[i]+2.0*mu[i])
            S_starL = ((lam[i]/(2.0*mu[i]))+1.0)*(Sigy[i]*np.sign(SS_star_L)+ KK_L*EP[i])
            alpha1 = (S_starL-SL)/(rho[i]*matC[i])
            R = U[i+1,:]-U[i,:]-alpha1*np.array([rho[i]*matC[i],1.0])
            alphaA = computeAlpha(R,rho[i],rho[i+1],matCp[i],matC[i+1])
            waves[i,:,0] = alpha1*np.array([rho[i]*matC[i],1.0])
            waves[i,:,1] = alphaA[0]*np.array([rho[i]*matCp[i],1.0])
            waves[i,:,3] = alphaA[1]*np.array([-rho[i+1]*matC[i+1],1.0])
            s[i,0] = -matC[i] ; s[i,1] = -matCp[i] ; s[i,3] = matC[i+1]
            amdq[i,:] = s[i,0]*waves[i,:,0]+s[i,1]*waves[i,:,1]
            apdq[i,:] = s[i,3]*waves[i,:,3]
        elif (fL<0.0) and (fR>0.0):
            #Rightward plastic wave
            KK_R = 3.0*(H[i+1]/2.0) +(mu[i+1]*(3.0*lam[i+1]+2.0*mu[i+1])) \
                   /(lam[i+1]+2.0*mu[i+1])
            S_starR = ((lam[i+1]/(2.0*mu[i+1]))+1.0)*(Sigy[i+1]*np.sign(SS_star_R) \
                                                      + KK_R*EP[i+1])
            alpha2 = (S_starR-SR)/(rho[i+1]*matC[i+1])
            R = U[i+1,:]-U[i,:]-alpha2*np.array([-rho[i+1]*matC[i+1],1.0])
            alphaA = computeAlpha(R,rho[i],rho[i+1],matC[i],matCp[i+1])
            waves[i,:,0] = alphaA[0]*np.array([rho[i]*matC[i],1.0])
            waves[i,:,2] = alphaA[1]*np.array([-rho[i+1]*matCp[i+1],1.0])
            waves[i,:,3] = alpha2*np.array([-rho[i+1]*matC[i+1],1.0])
            s[i,0] = -matC[i] ; s[i,2] = matCp[i+1] ; s[i,3] = matC[i+1]
            amdq[i,:] = s[i,0]*waves[i,:,0]
            apdq[i,:] = s[i,2]*waves[i,:,2]+s[i,3]*waves[i,:,3]
        elif (fL>0.0) and (fR>0.0):
            #Four waves
            KK_L = 3.0*(H[i]/2.0) +(mu[i]*(3.0*lam[i]+2.0*mu[i]))/(lam[i]+2.0*mu[i])
            S_starL = ((lam[i]/(2.0*mu[i]))+1.0)*(Sigy[i]*np.sign(SS_star_L) \
                                                  + KK_L*EP[i])
            KK_R = 3.0*(H[i+1]/2.0) +(mu[i+1]*(3.0*lam[i+1]+2.0*mu[i+1])) \
                   /(lam[i+1]+2.0*mu[i+1])
            S_starR = ((lam[i+1]/(2.0*mu[i+1]))+1.0)*(Sigy[i+1]*np.sign(SS_star_R) \
                                                      + KK_R*EP[i+1])
            alpha1 = (S_starL-SL)/(rho[i]*matC[i])
            alpha2 = (S_starR-SR)/(rho[i+1]*matC[i+1])
            R = U[i+1,:]-U[i,:]-alpha1*np.array([rho[i]*matC[i],1.0])-alpha2*np.array([-rho[i+1]*matC[i+1],1.0])
            alphaP = computeAlpha(R,rho[i],rho[i+1],matCp[i],matCp[i+1])
            waves[i,:,0] = alpha1*np.array([rho[i]*matC[i],1.0])
            waves[i,:,1] = alphaP[0]*np.array([rho[i]*matCp[i],1.0])
            waves[i,:,2] = alphaP[1]*np.array([-rho[i+1]*matCp[i+1],1.0])
            waves[i,:,3] = alpha2*np.array([-rho[i+1]*matC[i+1],1.0])
            s[i,0] = -matC[i] ; s[i,1] = -matCp[i] ; s[i,2] = matCp[i+1] ; s[i,3] = matC[i+1]
            amdq[i,:] = s[i,0]*waves[i,:,0]+s[i,1]*waves[i,:,1]
            apdq[i,:] = s[i,2]*waves[i,:,2]+s[i,3]*waves[i,:,3]
        else:
            #Elastic fluctuations
            waves[i,:,0] = alpha[0]*np.array([rho[i]*matC[i],1])
            waves[i,:,3] = alpha[1]*np.array([-rho[i+1]*matC[i+1],1])
            s[i,0] = -matC[i] ; s[i,3] = matC[i+1]
            amdq[i,:] = s[i,0]*waves[i,:,0]
            apdq[i,:] = s[i,3]*waves[i,:,3]
    return amdq,apdq,waves,s

def conservativeUpdate(U,dt0dx,amdq,apdq,cqxx,limit):
    for i in range(Nelem+2)[2:]:
        U[i,:] -= dt0dx*(apdq[i-1,:]+amdq[i,:])
        if (limit):
            U[i,:] -= (dt0dx/2.0)*(cqxx[i,:]-cqxx[i-1,:])


def computeEP(U,EPn,lam,mu,H,Sigy):
    EP = np.copy(EPn)
    for i in range(Nelem+3):
        SS = computeSeq(U[i,0],EPn[i],lam[i],mu[i],H[i])
        f = np.abs(SS)-Sigy[i]
        if (f>0):
            KK = 3.0*(H[i]/2.0) +(mu[i]*(3.0*lam[i]+2.0*mu[i]))/(lam[i]+2.0*mu[i])
            EP[i] = (-Sigy[i]*np.sign(SS)+((2.0*mu[i])/(lam[i]+2.0*mu[i]))*U[i,0])/KK
    return EP

CFL = 1.0

for i in range(NTmaxi)[1:]:
    #Apply Courant condition
    dt = CFL* dx/c_A
    if ((t[i-1]+dt)>timeOut):
        dt = timeOut - t[i-1]
    t[i]=t[i-1]+dt
    ##Apply boundary conditions: 
    #Left
    U[0,0] = -U[3,0] ; U[0,1] = U[3,1]
    U[1,0] = -U[2,0] ; U[1,1] = U[2,1]
    #Right
    U[Nelem+1+1,0] = -U[Nelem+1,0] ; U[Nelem+1+1,1] = U[Nelem+1,1]
    U[Nelem+2+1,0] = -U[Nelem-1+1,0] ; U[Nelem+2+1,1] = U[Nelem-1+1,1]
    ###Elastic-plastic Riemann solver and computation of fluctuations
    amdq,apdq,waves,s = EP_RiemannSolver(U,EP[:,i-1],rho,lam,mu,Sigy,H)
    #Compute limited waves
    cqxx = computeLimitedWaves(waves,s,mthlim,dt/dx)
    #Conservative update and update of the plastic strain
    conservativeUpdate(U,dt/dx,amdq,apdq,cqxx,limit)
    #Compute Updated Plastic Strain
    EP[:,i] = computeEP(U,EP[:,i-1],lam,mu,H,Sigy)
    #Store results
    sig[:,i] = U[:,0]
    v[:,i] = U[:,1]
    #Update of other fields
    dEPeq = (EP[:,i] - EP[:,i-1])*np.sign(sig[:,i])
    EPeq[:,i] = EPeq[:,i-1] + dEPeq
    eps[:,i] = sig[:,i]/E_A + EP[:,i]
    disp[:,i] = disp[:,i-1] + v[:,i]*dt
    if (t[i]==timeOut):
        increments=i
        break

increments=i
sig=sig[2:-2,:increments]
EP=EP[2:-2,:increments]
v=v[2:-2,:increments]
