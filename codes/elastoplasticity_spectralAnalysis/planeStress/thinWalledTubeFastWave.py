#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy.integrate import ode
import scipy.optimize
import pdb

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

def computeCliftonTangentFastSig(sig22,sig12,lamb,mu,h):
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
    
    return -A/B

def computeCliftonTangentFastTau(sig12,sig22,lamb,mu,h):
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
    B = (1./E)+H*(sig22/thet)**2-1./rcf2
    A = H*sig22*sig12
    return -A/B

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

# def computeTangentFastSig(sig11,sig12,rho,nu,E,lamb,mu,beta,h):
#     #deviatoric stresses
#     s12=sig12 ; s11 = (2./3.)*sig11 ; s22 = -(1./3.)*sig11
#     sigdnorm2=2.*s12**2 + s11**2 + 2.*s22**2
#     BETA=beta/sigdnorm2

#     C1111=lamb+2.*mu-BETA*s11**2 ; C1122=lamb-BETA*s11*s22 ; C1212 = mu-BETA*s12**2
#     C1222=-BETA*s12*s22 ; C1112=-BETA*s11*s12
#     H2211=-(nu/E+sig11**2/(3.*h*sigdnorm2)) 
#     H2212=-0.5*sig11*s12/(h*sigdnorm2)
#     A=np.matrix([[rho,0,0,0],[0.,1.-2.*H2211*C1122,0.,-4.*H2212*C1122],[0,0.,rho,0.],[0.,-2.*H2211*C1222,0.,1-4.*H2212*C1222]])
#     B=-np.array([[0.,1.,0.,0.],[C1111,0.,C1112,0.],[0,0.,0.,1.],[C1112,0.,C1212,0.]])
#     Ainv=np.linalg.inv(A)
#     Jac=np.dot(Ainv,B)
        
#     characteristiStructure = np.linalg.eig(Jac)
#     eigenvalues= characteristiStructure[0]
#     eigenvectors= characteristiStructure[1]
        
#     # Find fast and slow celerities
#     cf=max(eigenvalues);cs=np.min(np.abs(eigenvalues))
    
#     slow=np.where(np.abs(eigenvalues-cs)<1.e-2)[0][0]
#     slow2=np.where(np.abs(eigenvalues+cs)<1.e-2)[0][0]
#     if (eigenvectors[:,slow2]>0.).all():eigenvectors[:,slow2]*=-1.
        
#     dq=eigenvectors[:,slow]-eigenvectors[:,slow2]
#     Adsig=dq[1] ; Bdtau=dq[3]
#     return -Adsig/Bdtau


def computeEigenStresses(sig):
    #    | sig11 sig12   0   |
    #sig=| sig12 sig22   0   |
    #    |   0     0   sig33 |
    s3=sig[2,2]
    delta=(sig[0,0]-sig[1,1])**2+4.*sig[0,1]**2
    s1=0.5*(sig[0,0]+sig[1,1]-np.sqrt(delta))
    s2=0.5*(sig[0,0]+sig[1,1]+np.sqrt(delta))
    return np.array([s1,s2,s3])

def computeTangentFastSig(sig11,sig12,lamb,mu,beta):
    # sig12 driven
    n1=1.;n2=0.
    sig22=0.;sig33=0.
    H=tangentModulus(np.array([sig11,sig12,sig22,sig33]),lamb,mu,beta,'thinWalled')
    C=acousticTensor(H,np.array([n1,n2]))
    eigenf,eigens=acousticEigenStructure(C)
    alpha11= (H[0,1]*n1+H[1,1]*n2)*(H[1,2]*n1+H[2,2]*n2) - (H[0,2]*n1+H[1,2]*n2)*(H[1,1]*n1+H[1,2]*n2)
    alpha12=((H[0,1]*n1+H[0,2]*n2)*(H[0,2]*n1+H[1,2]*n2) - (H[0,0]*n1+H[0,1]*n2)*(H[1,2]*n1+H[2,2]*n2))/2.
    alpha22= (H[0,0]*n1+H[0,1]*n2)*(H[1,1]*n1+H[1,2]*n2) - (H[0,1]*n1+H[0,2]*n2)*(H[0,1]*n1+H[1,1]*n2)
    w1=eigens[1][0];w2=eigens[1][1]
    psi11=-w1/w2
    return psi11

def computeTangentFastTau(sig12,sig11,lamb,mu,beta):
    # sig12 driven
    n1=1.;n2=0.
    sig22=0.;sig33=0.
    H=tangentModulus(np.array([sig11,sig12,sig22,sig33]),lamb,mu,beta,'thinWalled')
    C=acousticTensor(H,np.array([n1,n2]))
    eigenf,eigens=acousticEigenStructure(C)
    alpha11= (H[0,1]*n1+H[1,1]*n2)*(H[1,2]*n1+H[2,2]*n2) - (H[0,2]*n1+H[1,2]*n2)*(H[1,1]*n1+H[1,2]*n2)
    alpha12=((H[0,1]*n1+H[0,2]*n2)*(H[0,2]*n1+H[1,2]*n2) - (H[0,0]*n1+H[0,1]*n2)*(H[1,2]*n1+H[2,2]*n2))/2.
    alpha22= (H[0,0]*n1+H[0,1]*n2)*(H[1,1]*n1+H[1,2]*n2) - (H[0,1]*n1+H[0,2]*n2)*(H[0,1]*n1+H[1,1]*n2)
    w1=eigens[1][0];w2=eigens[1][1]
    psi11=-w2/(1.*w1)
    return psi11

def integrateODE_tauDriven(dtau,sig0,tau0,lamb,mu,beta):
    sigma=sig0
    # computeTangentSlowTau(sig12,sig11,lamb,mu,beta)
    # subdivision of time step
    sub_steps = 1
    dTAU = dtau/sub_steps
    theta = 1.0
    for i in range(sub_steps):
        ## Nonlinear solution procedure
        ## R = s^{n+1} - s^{n} - RHS
        R=lambda x: x - sigma - dTAU*(theta*computeTangentFastTau(tau0+dTAU,x,lamb,mu,beta)+(1.0-theta)*computeTangentFastTau(tau0,sigma,lamb,mu,beta))
        solution = scipy.optimize.fsolve(R,sigma)
        sigma = solution
    return solution

def integrateODE_sigDriven(dsig,sig0,tau0,lamb,mu,beta):
    sigma=tau0
    # computeTangentSlowTau(sig12,sig11,lamb,mu,beta)
    # subdivision of time step
    sub_steps = 1
    dSIG = dsig/sub_steps
    theta = 1.0
    for i in range(sub_steps):
        ## Nonlinear solution procedure
        ## R = s^{n+1} - s^{n} - RHS
        R=lambda x: x - sigma - dSIG*(theta*computeTangentFastSig(sig0+dSIG,x,lamb,mu,beta)+(1.0-theta)*computeTangentFastSig(sig0,sigma,lamb,mu,beta))
        #solution = scipy.optimize.fsolve(R,sigma)
        solution = scipy.optimize.root(R,sigma).x[0]
        sigma = solution
    return solution

rho = 7800.
E = 2.e11
nu = 0.3
mu = 0.5*E/(1.+nu)
kappa = E/(3.*(1.-2.*nu))
lamb = kappa-2.*mu/3.
sigy = 100.0e6        
H = 100.08e6
beta=(6.*mu**2)/(3.*mu+H)
Niter=10000

## Initial yield surface
Samples=100
s11=0.*sigy
s12=np.linspace(0.,np.sqrt(1./3.)*sigy ,Samples)
s22=np.zeros(len(s12))
s22=np.zeros(Samples)
for i in range(Samples-1):
    fvm = lambda x : criterion(s12[i],x,0.,sigy,0.)
    s22[i] = scipy.optimize.brentq(fvm,0.*sigy,2.5*sigy)
fileName=path+'TWfast_yield0.pgf'
export2pgfPlotFile(fileName,np.array([s12,s22]),'sigma_12','sigma_11')
cylindre=vonMisesYieldSurface(sigy)
fileName=path+'TWCylindreDevPlane.pgf'
export2pgfPlot3D(str(fileName),cylindre[0,:],cylindre[1,:],cylindre[2,:])



parameter="tau"
parameter="sig"




TAU=np.zeros(Niter)    
SIG=np.zeros(Niter)    
TAUC=np.zeros(Niter)    
SIGC=np.zeros(Niter)    

## analysis in deviator plane
sigdev1=np.zeros(np.shape(TAU))
sigdev2=np.zeros(np.shape(TAU))
sigdev3=np.zeros(np.shape(TAU))
sigdev1C=np.zeros(np.shape(TAU))
sigdev2C=np.zeros(np.shape(TAU))
sigdev3C=np.zeros(np.shape(TAU))





   
if parameter=='tau':
    sig0=0.*sigy
    tau0=np.sqrt((-sig0**2+ sigy**2)/3.)
    tauEnd=0.0*sigy


    dtau=(tauEnd-tau0)/Niter
    TAU=np.linspace(tau0,tauEnd,Niter)
    SIG[0]=sig0
    r = ode(computeTangentFastTau).set_integrator('vode',method='bdf',order=5)
    r.set_initial_value(SIG[0],TAU[0]).set_f_params(lamb,mu,beta)
    
    for i in range(Niter-1):
        # r.set_f_params(SIG[i],lamb,mu,beta)
        # r.integrate(r.t+dtau)
        
        # SIG[i+1]=r.y
    
        #SIG[i+1]=integrateODE_tauDriven(dtau,SIG[i],TAU[i],lamb,mu,beta)

        sigDev=computeDeviatoricPart(np.array([SIG[i+1],TAU[i+1],0.,0.]))
        sigma = np.matrix([[sigDev[0],sigDev[1]/np.sqrt(2.),0.],[sigDev[1]/np.sqrt(2.),sigDev[2],0.],[0.,0.,sigDev[3]]])
        eigSigDev=computeEigenStresses(sigma)

        sigdev1[i+1]=eigSigDev[0]
        sigdev2[i+1]=eigSigDev[1]
        sigdev3[i+1]=eigSigDev[2]

    # Clifton solution
    TAUC=np.linspace(tau0,tauEnd,Niter)
    SIGC[0]=sig0

    r = ode(computeCliftonTangentFastTau).set_integrator('vode',method='bdf',order=5)
    r.set_initial_value(SIGC[0],TAUC[0]).set_f_params(lamb,mu,H)

    for i in range(Niter-1):
        r.set_f_params(SIGC[i],lamb,mu,H)
        r.integrate(r.t+dtau)
        
        SIGC[i+1]=r.y
    
        sigDev=computeDeviatoricPart(np.array([SIGC[i+1],TAUC[i+1],0.,0.]))
        sigma = np.matrix([[sigDev[0],sigDev[1]/np.sqrt(2.),0.],[sigDev[1]/np.sqrt(2.),sigDev[2],0.],[0.,0.,sigDev[3]]])
        eigSigDev=computeEigenStresses(sigma)
        
        sigdev1C[i+1]=eigSigDev[0]
        sigdev2C[i+1]=eigSigDev[1]
        sigdev3C[i+1]=eigSigDev[2]
elif parameter=='sig':

    sig0=0.000*sigy
    tau0=np.sqrt((-sig0**2+ sigy**2)/3.)
    sigEnd=1.*sigy

    sigDev=computeDeviatoricPart(np.array([sig0,tau0,0.,0.]))
    sigma = np.matrix([[sigDev[0],sigDev[1]/np.sqrt(2.),0.],[sigDev[1]/np.sqrt(2.),sigDev[2],0.],[0.,0.,sigDev[3]]])
    eigSigDev=computeEigenStresses(sigma)

    sigdev1[0]=eigSigDev[0]
    sigdev2[0]=eigSigDev[1]
    sigdev3[0]=eigSigDev[2]
    sigdev1C[0]=eigSigDev[0]
    sigdev2C[0]=eigSigDev[1]
    sigdev3C[0]=eigSigDev[2]
    
    dsig=(sigEnd-sig0)/Niter
    SIG=np.linspace(sig0,sigEnd,Niter)
    TAU[0]=tau0
    
    r = ode(computeTangentFastSig).set_integrator('vode',method='bdf',order=5)
    r.set_initial_value(TAU[0],SIG[0]).set_f_params(lamb,mu,beta)

    for i in range(Niter-1):
        # r.set_f_params(TAU[i],lamb,mu,beta)
        # r.integrate(r.t+dsig)
        
        # TAU[i+1]=r.y
    
        TAU[i+1]=integrateODE_sigDriven(dsig,SIG[i],TAU[i],lamb,mu,beta)

        sigDev=computeDeviatoricPart(np.array([SIG[i+1],TAU[i+1],0.,0.]))
        sigma = np.matrix([[sigDev[0],sigDev[1]/np.sqrt(2.),0.],[sigDev[1]/np.sqrt(2.),sigDev[2],0.],[0.,0.,sigDev[3]]])
        eigSigDev=computeEigenStresses(sigma)

        sigdev1[i+1]=eigSigDev[0]
        sigdev2[i+1]=eigSigDev[1]
        sigdev3[i+1]=eigSigDev[2]

    # Clifton solution
    SIGC=np.linspace(sig0,sigEnd,Niter)
    TAUC[0]=tau0

    r = ode(computeCliftonTangentFastSig).set_integrator('vode',method='bdf',order=5)
    r.set_initial_value(TAUC[0],SIGC[0]).set_f_params(lamb,mu,H)

    for i in range(Niter-1):
        r.set_f_params(TAUC[i],lamb,mu,H)
        r.integrate(r.t+dsig)
        
        TAUC[i+1]=r.y
    
        sigDev=computeDeviatoricPart(np.array([SIGC[i+1],TAUC[i+1],0.,0.]))
        sigma = np.matrix([[sigDev[0],sigDev[1]/np.sqrt(2.),0.],[sigDev[1]/np.sqrt(2.),sigDev[2],0.],[0.,0.,sigDev[3]]])
        eigSigDev=computeEigenStresses(sigma)
        
        sigdev1C[i+1]=eigSigDev[0]
        sigdev2C[i+1]=eigSigDev[1]
        sigdev3C[i+1]=eigSigDev[2]



    
fileName=path+'fastStressPlane_Stress.pgf'
export2pgfPlotFile(fileName,np.array([TAU[0:-1:Niter/100],SIG[0:-1:Niter/100]]),'sigma_12','sigma_11')

fileName=path+'TWfastStressPlane_Stress.pgf'
export2pgfPlotFile(fileName,np.array([TAUC[0:-1:Niter/100],SIGC[0:-1:Niter/100]]),'sigma_12','sigma_11')

fileName=path+'fastDevPlane_Stress.pgf'
dico={"xlabel":r'$s_1$',"ylabel":r'$s_2$',"zlabel":r'$s_3$'}
export2pgfPlot3D(fileName,sigdev1[0:-1:Niter/100],sigdev2[0:-1:Niter/100],sigdev3[0:-1:Niter/100],dico)

fileName=path+'TWfastDevPlane_Stress.pgf'
dico={"xlabel":r'$s_1$',"ylabel":r'$s_2$',"zlabel":r'$s_3$'}
export2pgfPlot3D(fileName,sigdev1C[0:-1:Niter/100],sigdev2C[0:-1:Niter/100],sigdev3C[0:-1:Niter/100],dico)

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
ax2d.set_ylim([0,np.max(s12)])
ax2d.set_xlim([0.,1.1*sigEnd])

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

ax2d.plot(SIG,TAU,lw=2.5,linestyle="--",label="Fast wave ")
ax2d.plot(SIGC,TAUC,lw=1.,label="Fast wave (Clifton)")
ax.plot(sigdev1,sigdev2,sigdev3,lw=2.5,linestyle="--")
ax.plot(sigdev1C,sigdev2C,sigdev3C,lw=1.)

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
plt.show()

