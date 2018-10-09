#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy.integrate import ode
import scipy.optimize
import pdb
from buildTeXFiles import *
import os
import sys

directory=os.path.basename(__file__)[:22]
if not os.path.exists('pgf_'+str(directory)+'/'):
    os.system('mkdir pgf_'+str(directory)+'/')
path='pgf_'+str(directory)+'/'

def export2pgfPlot3D(fileName,field1,field2,field3,dico={"xlabel":'x',"ylabel":'y',"zlabel":'z'}):
    #pdb.set_trace()
    dataFile=open(fileName,"w")
    # 2D export
    xlabel=dico["xlabel"]
    ylabel=dico["ylabel"]
    zlabel=dico["zlabel"]
    dataFile.write('# Curve ('+str(xlabel)+';'+str(ylabel)+';'+str(zlabel)+') '+str(field1.shape[0])+' points.\n')
    for i,x in enumerate(field1):
        dataFile.write(str(x)+' '+str(field2[i])+' '+str(field3[i])+' i\n')
    dataFile.close()

def export2pgfPlotFile(fileName,fields,*kwargs):
    #pdb.set_trace()
    dataFile=open(fileName,"w")
    # 2D export
    n_fields = np.shape(fields)[0]
    n_labels = np.shape(kwargs)[0]
    if n_fields != n_labels : print n_fields-n_labels," missing in export to pgf File !"
    labels=[]
    line1 = ' # Curve ('
    line2 = ''
    for i in range(n_labels):
        labels.append(str(kwargs[i]))
        line1+=str(kwargs[i])+' ; '
        line2+=str(kwargs[i])+' '
        if i==n_labels-1 :
            line1+=str(np.shape(fields)[1])+' points.)\n'
            line2+=' \n'
        
    dataFile.write(line1)
    dataFile.write(line2)
    for i,x in enumerate(fields[0,:]):
        dataFile.write(str(x))
        for j in range(n_fields-1):
            dataFile.write(' '+str(fields[j+1,i]))
        dataFile.write(' i\n')
    dataFile.close()



def vonMisesYieldSurface(sigy):
    radius=np.sqrt((2./3.)*sigy**2)
    theta=np.linspace(0,2*np.pi,50)
    s2 = radius*np.cos(theta)
    s3 = radius*np.sin(theta)
    s1=0.
    c=np.sqrt(2.)/2.;
    s=np.sqrt(2.)/2.;
    P2=np.array([[c,-c,0.],[c,c,0.],[0.,0.,1.]])
    P1=np.array([[c,0.,-c],[0.,1.,0.],[c,0.,c]])
    c=np.cos(np.arctan(1./np.sqrt(2.0)))
    s=np.sin(np.arctan(1./np.sqrt(2.0)))
    P1=np.array([[c,0.,-s],[0.,1.,0.],[s,0.,c]])
    cylindre=np.zeros((3,len(s2)))

    for i in range(len(s2)):
        cylindre[:,i] = np.dot(P2,np.dot(P1,np.array([s1,s2[i],s3[i]])))
    return cylindre

def tangentModulus(sigma,lamb,mu,beta,tangent):
    H=np.zeros((3,3))
    #    |H1111 H1112 H1122|
    # H =|H1211 H1212 H1222|
    #    |H2211 H2212 H2222|
    # sigma = [sig11 , sig12 , sig22 , sig33 ]
    sigDev = computeDeviatoricPart(sigma)
    sigdnorm2=np.dot(sigDev.T,sigDev)
    BETA=beta/sigdnorm2
    s11=sigDev[0];s12=sigDev[1]/np.sqrt(2.);s22=sigDev[2];s33=sigDev[3]
    
    ## Plane stress tangent modulus Hijkl = Hijkl - Hij33*H33kl/H3333
    H1133=(lamb -BETA*s11*s33)
    H1233=(-BETA*s12*s33)
    H1122=(lamb -BETA*s11*s22)
    H2222=(lamb+2.*mu -BETA*s22**2)
    H1222=(-BETA*s12*s22)
    H2233=(lamb-BETA*s22*s33)
    H3333=(lamb+2.*mu-BETA*s33*s33)
    if tangent=='planeStress':
        H[0,0]=lamb+2.*mu - BETA*s11**2 -H1133*H1133/H3333
        H[0,1]=-BETA*s11*s12 -H1133*H1233/H3333
        H[0,2]=lamb-BETA*s11*s22 -H1133*H2233/H3333
        H[1,0]=-BETA*s12*s11-H1233*H1133/H3333
        H[1,1]=mu-BETA*s12**2 -H1233*H1233/H3333
        H[1,2]=-BETA*s12*s22-H1233*H2233/H3333
        H[2,0]=lamb - BETA*s11*s22 -H2233*H1133/H3333
        H[2,1]=-BETA*s22*s12 -H2233*H1233/H3333
        H[2,2]=lamb+2.*mu-BETA*s22**2 -H2233*H2233/H3333
    elif tangent=='thinWalled':
        H[0,0]=lamb+2.*mu - BETA*s11**2 -H1122*(H1122+H1133)/(H2233+H2222)
        H[0,1]=-BETA*s11*s12  -H1222*(H1122+H1133)/(H2233+H2222)
        H[0,2]=lamb-BETA*s11*s22 
        H[1,0]=-BETA*s12*s11-H1122*(H1222+H1233)/(H2233+H2222)
        H[1,1]=mu-BETA*s12**2-H1222*(H1222+H1233)/(H2233+H2222)
        H[1,2]=-BETA*s12*s22
        H[2,0]=lamb - BETA*s11*s22 
        H[2,1]=-BETA*s22*s12
        H[2,2]=lamb+2.*mu-BETA*s22**2
    else :
        H[0,0]=lamb+2.*mu - BETA*s11**2 
        H[0,1]=-BETA*s11*s12
        H[0,2]=lamb-BETA*s11*s22
        H[1,0]=-BETA*s12*s11
        H[1,1]=mu-BETA*s12**2
        H[1,2]=-BETA*s12*s22
        H[2,0]=lamb-BETA*s11*s22
        H[2,1]=-BETA*s12*s22
        H[2,2]=lamb+2.*mu-BETA*s22**2
    return H

def acousticTensor(H,n):
    n1=n[0] ; n2=n[1]
    C11 = H[0,0]*n1**2 + H[1,1]*n2**2 + 2.*H[0,1]*n1*n2
    C12 = H[0,1]*n1**2 + H[1,2]*n2**2 + (H[0,2]+H[1,1])*n1*n2
    C22 = H[1,1]*n1**2 + H[2,2]*n2**2 + 2.*H[2,1]*n1*n2
    return np.array([C11,C12,C22])


def acousticEigenStructure(C):
    C11=C[0];C12=C[1];C22=C[2]
    ## omega1,w1 associated to cf
    ## omega2,w2 associated to cs
    omega1=0.5*(C11+C22 + np.sqrt((C11-C22)**2+4.*C12**2))
    omega2=0.5*(C11+C22 - np.sqrt((C11-C22)**2+4.*C12**2))
    #w1=np.array([-C12,C11-omega1])
    w1=np.array([C22-omega1,-C12])
    w2=np.array([-C12,C11-omega2])
    return [omega1,w1],[omega2,w2]

def computeDeviatoricPart(T):
    # T = [T11 T21 T22 T33]
    Pdev=np.array([[1.-1/3.,0.,-1./3.,-1./3.],[0.,1.,0.,0.],[-1./3.,0.,1.-1./3.,-1./3.],[-1./3.,0.,-1./3.,1.-1./3.]])
    Tdev=np.dot(Pdev,T)
    return np.array([Tdev[0],np.sqrt(2.)*Tdev[1],Tdev[2],Tdev[3]])

def criterion(s12,s22,p,sigy,H):
    f=np.sqrt(3./2.)*np.sqrt(2.*s12**2 +((2./3.)*s22)**2 + 2.*((-1./3.)*s22)**2 ) - (sigy+H*p)
    return f


def computeCliftonTangentSlowTau(sig12,sig22,lamb,mu,h):
    #deviatoric stresses
    s12=sig12 ; s22 = (2./3.)*sig22 ; s11 = -(1./3.)*sig22
    sigdnorm2=2.*s12**2 + s22**2 + 2.*s11**2
    thet = np.sqrt(3.)
    E = mu*(3.*lamb+2.*mu)/(lamb+mu)
    H=2./(h*sigdnorm2)
    M = 1./mu + H*(thet*sig12)**2 ; N = 1./E + H*(sig22/thet)**2
    L = 1./(mu*E) + (1./mu)*H*(sig22/thet)**2 + (1./E)*H*(thet*sig12)**2
    rcf2 = 0.5*(M+N+np.sqrt((M-N)**2+4.*(H*sig12*sig22)**2))/L
    rcs2 = 0.5*(M+N-np.sqrt((M-N)**2+4.*(H*sig12*sig22)**2))/L
    B = (1./mu)+H*(sig12*thet)**2-1./rcf2
    A = H*sig22*sig12
    return A/B

def computeEigenStresses(sig):
    #    | sig11 sig12   0   |
    #sig=| sig12 sig22   0   |
    #    |   0     0   sig33 |
    s3=sig[2,2]
    delta=(sig[0,0]-sig[1,1])**2+4.*sig[0,1]**2
    s1=0.5*(sig[0,0]+sig[1,1]-np.sqrt(delta))
    s2=0.5*(sig[0,0]+sig[1,1]+np.sqrt(delta))
    return np.array([s1,s2,s3])

def computeTangentSlowTau(sig12,sig11,lamb,mu,beta):
    n1=1.;n2=0.
    H=tangentModulus(np.array([sig11,sig12,0.,0.]),lamb,mu,beta,'thinWalled')
    C=acousticTensor(H,np.array([n1,n2]))
    eigenf,eigens=acousticEigenStructure(C)
    w1=eigenf[1][0];w2=eigenf[1][1]
    psi11=-w2/(w1)
    return psi11

def integrateODE(dtau,sig0,tau0,lamb,mu,beta):
    sigma=sig0
    # computeTangentSlowTau(sig12,sig11,lamb,mu,beta)
    # subdivision of time step
    #pdb.set_trace()
    sub_steps = 1
    dTAU = dtau/sub_steps
    theta = 1.0
    for i in range(sub_steps):
        ## Nonlinear solution procedure
        ## R = s^{n+1} - s^{n} - RHS
        R=lambda x: x - sigma - dTAU*(theta*computeTangentSlowTau(tau0+dTAU,x,lamb,mu,beta)+(1.0-theta)*computeTangentSlowTau(tau0,sigma,lamb,mu,beta))
        #pdb.set_trace()
        solution = scipy.optimize.fsolve(R,sigma)
        sigma = solution
    return solution

rho = 7800.
E = 2.e11
nu = 0.3
mu = 0.5*E/(1.+nu)
kappa = E/(3.*(1.-2.*nu))
lamb = kappa-2.*mu/3.
sigy = 100.0e6        
H = 100.0e6
beta=(6.*mu**2)/(3.*mu+H)
Niter=1000


# purple to red
col=['#781C81','#3F60AE','#539EB6','#6DB388','#CAB843','#E78532','#D92120']
col=['#332288','#88CCEE','#44AA99','#117733','#999933','#DDCC77','#CC6677','#882255','#AA4499']


sig0=np.array([0.,0.25*sigy,0.5*sigy,0.75*sigy,0.85*sigy,0.9*sigy,0.99999*sigy])

tau0=np.sqrt((sigy**2-sig0**2)/3.)

parameter="tau"

tauEnd=0.75*sigy
sigEnd=1.*sigy


TAU=np.zeros((len(sig0),Niter))    
SIG=np.zeros((len(sig0),Niter))    
TAUC=np.zeros((len(sig0),Niter))    
SIGC=np.zeros((len(sig0),Niter))    

## analysis in deviator plane
sigdev1=np.zeros((len(sig0),Niter))
sigdev2=np.zeros((len(sig0),Niter))
sigdev3=np.zeros((len(sig0),Niter))
sigdev1C=np.zeros((len(sig0),Niter))
sigdev2C=np.zeros((len(sig0),Niter))
sigdev3C=np.zeros((len(sig0),Niter))

              
## LOADING PATHS PLOTS
pgfFilesList=[]
yields11_s12=[]
deviatorPlots=[]

## Initial yield surface
Samples=100
s12=np.linspace(0.,np.sqrt(1./3.)*sigy ,Samples)
s22=np.zeros(len(s12))
s22=np.zeros(Samples)
for i in range(Samples-1):
    fvm = lambda x : criterion(s12[i],x,0.,sigy,0.)
    s22[i] = scipy.optimize.brentq(fvm,0.*sigy,2.5*sigy)
fileName=path+'TWslow_yield0.pgf'
export2pgfPlotFile(fileName,np.array([s12,s22]),'sigma_12','sigma_11')
#pgfFilesList.append(fileName)
cylindre=vonMisesYieldSurface(sigy)
fileName=path+'TWCylindreDevPlane.pgf'
export2pgfPlot3D(str(fileName),cylindre[0,:],cylindre[1,:],cylindre[2,:])
#deviatorPlots.append(fileName)


for k,s in enumerate(sig0):
    # sig = np.matrix([[s,tau0[k],0.],[tau0[k],0.,0.],[0.,0.,0.]])
    # eigSigDev=computeEigenStresses(sig)
    # sigDev=computeDeviatoricPart(np.array([s,tau0[k],0.,0.]))
    # sigma = np.matrix([[sigDev[0],sigDev[1]/np.sqrt(2.),0.],[sigDev[1]/np.sqrt(2.),sigDev[2],0.],[0.,0.,sigDev[3]]])
    sigma = np.matrix([[s,tau0[k],0.],[tau0[k],0.,0.],[0.,0.,0.]])
    eigSigDev=computeEigenStresses(sigma)
    sigdev1[k,0]=eigSigDev[0]
    sigdev2[k,0]=eigSigDev[1]
    sigdev3[k,0]=eigSigDev[2]
    sigdev1C[k,0]=eigSigDev[0]
    sigdev2C[k,0]=eigSigDev[1]
    sigdev3C[k,0]=eigSigDev[2]
    
    dtau=(tauEnd-tau0[k])/Niter
    TAU[k,:]=np.linspace(tau0[k],tauEnd,Niter)
    SIG[k,0]=s
    
    
    for i in range(Niter-1):
        # r.set_f_params(SIG[k,i],rho,nu,E,lamb,mu,beta,H)
        # r.integrate(r.t+dtau)
        
        # SIG[k,i+1]=r.y   
        SIG[k,i+1]=integrateODE(dtau,SIG[k,i],TAU[k,i],lamb,mu,beta)

        # sigDev=computeDeviatoricPart(np.array([SIG[k,i+1],TAU[k,i+1],0.,0.]))
        # sigma = np.matrix([[sigDev[0],sigDev[1]/np.sqrt(2.),0.],[sigDev[1]/np.sqrt(2.),sigDev[2],0.],[0.,0.,sigDev[3]]])
        sigma = np.matrix([[SIG[k,i+1],TAU[k,i+1],0.],[TAU[k,i+1],0.,0.],[0.,0.,0.]])
        eigSigDev=computeEigenStresses(sigma)

        # sig = np.matrix([[SIG[k,i+1],TAU[k,i+1],0.],[TAU[k,i+1],0.,0.],[0.,0.,0.]])
        # eigSigDev=computeEigenStresses(sig)
        sigdev1[k,i+1]=eigSigDev[0]
        sigdev2[k,i+1]=eigSigDev[1]
        sigdev3[k,i+1]=eigSigDev[2]
    
    # Clifton solution
    TAUC[k,:]=np.linspace(tau0[k],tauEnd,Niter)
    SIGC[k,0]=s
    
    r = ode(computeCliftonTangentSlowTau).set_integrator('vode',method='bdf',order=5)
    r.set_initial_value(SIGC[k,0],TAUC[k,0]).set_f_params(lamb,mu,H)
    
    for i in range(Niter-1):
        r.set_f_params(SIGC[k,i],lamb,mu,H)
        r.integrate(r.t+dtau)
        
        SIGC[k,i+1]=r.y

        sigDev=computeDeviatoricPart(np.array([SIG[k,i+1],TAU[k,i+1],0.,0.]))
        sigma = np.matrix([[sigDev[0],sigDev[1]/np.sqrt(2.),0.],[sigDev[1]/np.sqrt(2.),sigDev[2],0.],[0.,0.,sigDev[3]]])
        sigma = np.matrix([[SIG[k,i+1],TAU[k,i+1],0.],[TAU[k,i+1],0.,0.],[0.,0.,0.]])
        eigSigDev=computeEigenStresses(sigma)

        # sig = np.matrix([[SIGC[k,i+1],TAUC[k,i+1],0.],[TAUC[k,i+1],0.,0.],[0.,0.,0.]])
        # eigSigDev=computeEigenStresses(sig)
        sigdev1C[k,i+1]=eigSigDev[0]
        sigdev2C[k,i+1]=eigSigDev[1]
        sigdev3C[k,i+1]=eigSigDev[2]
    
    print path
    fileName=path+'slowStressPlane_Stress'+str(k)+'.pgf'
    export2pgfPlotFile(fileName,np.array([TAU[k,0:-1:Niter/100],SIG[k,0:-1:Niter/100]]),'sigma_12','sigma_11')
    pgfFilesList.append(fileName)
    fileName=path+'TWslowStressPlane_Stress'+str(k)+'.pgf'
    export2pgfPlotFile(fileName,np.array([TAUC[k,0:-1:Niter/100],SIGC[k,0:-1:Niter/100]]),'sigma_12','sigma_11')
    pgfFilesList.append(fileName)
    fileName=path+'slowDevPlane_Stress'+str(k)+'.pgf'
    dico={"xlabel":r'$s_1$',"ylabel":r'$s_2$',"zlabel":r'$s_3$'}
    export2pgfPlot3D(fileName,sigdev1[k,0:-1:Niter/100],sigdev2[k,0:-1:Niter/100],sigdev3[k,0:-1:Niter/100],dico)
    fileName=path+'TWslowDevPlane_Stress'+str(k)+'.pgf'
    dico={"xlabel":r'$s_1$',"ylabel":r'$s_2$',"zlabel":r'$s_3$'}
    export2pgfPlot3D(fileName,sigdev1C[k,0:-1:Niter/100],sigdev2C[k,0:-1:Niter/100],sigdev3C[k,0:-1:Niter/100],dico)
    deviatorPlots.append(fileName)

###################################### POST PROCESSING
from mpl_toolkits.mplot3d import proj3d
def orthogonal_proj(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    return np.array([[1,0,0,0],
                        [0,1,0,0],
                        [0,0,a,b],
                        [0,0,0,zback]])
proj3d.persp_transformation = orthogonal_proj

Samples=100
s11=0.*sigy
s12=np.linspace(0.,np.sqrt(1./3.)*sigy ,Samples)
s22=np.zeros(len(s12))
s22=np.zeros(Samples)
for i in range(Samples-1):
    fvm = lambda x : criterion(s12[i],x,0.,sigy,0.)
    s22[i] = scipy.optimize.brentq(fvm,0.*sigy,2.5*sigy)

ax2d = plt.subplot2grid((1,2),(0,0))
ax2d.set_xlabel(r'$\sigma$',size=24.)
ax2d.set_ylabel(r'$\tau$',size=24.)
ax2d.set_ylim([0.,sigy])
ax2d.set_xlim([0.,2.*sigy])

##3d subplot settings
ax = plt.subplot2grid((1,2),(0,1),projection='3d')
ax.set_aspect("equal")
ax.set_xlabel(r'$\sigma_1$',size=24.)
ax.set_ylabel(r'$\sigma_2$',size=24.)
ax.set_zlabel(r'$\sigma_3$',size=24.)

radius=np.sqrt((2./3.)*sigy**2)

ax.set_xlim(-1.*radius,1.*radius)
ax.set_ylim(-1.*radius,1.*radius)
ax.set_zlim(-1.*radius,1.*radius)

## 3d and 2d subplots
ax2d.plot(s22,s12,'k',lw=1.,label='Initial yield surface (vM)')
# draw von Mises Criterion
r=np.sqrt((2./3.)*sigy**2)
    
theta=np.linspace(0,2*np.pi,50)
s2 = r*np.cos(theta)
s3 = r*np.sin(theta)
s1=0.
c=np.sqrt(2.)/2.;
s=np.sqrt(2.)/2.;
P2=np.array([[c,-c,0.],[c,c,0.],[0.,0.,1.]])
P1=np.array([[c,0.,-c],[0.,1.,0.],[c,0.,c]])
c=np.cos(np.arctan(1./np.sqrt(2.0)))
s=np.sin(np.arctan(1./np.sqrt(2.0)))
P1=np.array([[c,0.,-s],[0.,1.,0.],[s,0.,c]])
cylindre=np.zeros((3,len(s2)))

for i in range(len(s2)):
    cylindre[:,i] = np.dot(P2,np.dot(P1,np.array([s1,s2[i],s3[i]])))
ax.plot(cylindre[0,:],cylindre[1,:],cylindre[2,:], color="k")
ax.plot([0.,sigy],[0.,sigy],[0.,sigy], color="k")

ax2d.grid()
ax2d.legend(numpoints=1,loc='best')
elevation_Angle_radian=np.arctan(1./np.sqrt(2.0))
angle_degree= 180.*elevation_Angle_radian/np.pi
ax.view_init(angle_degree,45.)
for k in range(len(sig0)):
    ax2d.plot(SIG[k,:],TAU[k,:],color=col[k],lw=2.5,linestyle="--",label="Slow wave ")
    ax2d.plot(SIGC[k,:],TAUC[k,:],color=col[k],lw=1.,label="Slow wave (Clifton)")
    ax.plot(sigdev1[k,:],sigdev2[k,:],sigdev3[k,:],color=col[k],lw=2.5,linestyle="--")
    ax.plot(sigdev1C[k,:],sigdev2C[k,:],sigdev3C[k,:],color=col[k],lw=1.)

#ax2d.legend(loc='best',numpoints=1)
plt.show()

