# !\usr\bin\python
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import scipy.optimize
from matplotlib import animation
from scipy.integrate import ode
import pdb
from matplotlib import rcParams
from buildTeXFiles import *
import os
import sys

directory=os.path.basename(__file__)[:20]
if not os.path.exists('pgf_'+str(directory)+'/'):
    os.system('mkdir pgf_'+str(directory)+'/')
path='pgf_'+str(directory)+'/'

rcParams['axes.labelsize'] = 16
rcParams['xtick.labelsize'] = 16
rcParams['ytick.labelsize'] = 16
rcParams['legend.fontsize'] = 16

# Material parameters
rho = 7800.
E = 2.e11
nu = 0.3
mu = 0.5*E/(1.+nu)
kappa = E/(3.*(1.-2.*nu))
lamb = kappa-2.*mu/3.
sigy = 100.0e6        
H = 100.0e6
#H = 100.0e8
beta=(6.*mu**2)/(3.*mu+H)


def export2pgfPlot2D(fileName,field1,field2,dico={"xlabel":'x',"ylabel":'y'}):
    #pdb.set_trace()
    dataFile=open(fileName,"w")
    xlabel=dico["xlabel"]
    ylabel=dico["ylabel"]
    dataFile.write('# Curve ('+str(xlabel)+';'+str(ylabel)+') '+str(len(field1))+' points.\n')
    for i,x in enumerate(field1):
        dataFile.write(str(x)+' '+str(field2[i])+' i\n')
    dataFile.close()

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
    sigdnorm2=np.dot(sigDev,sigDev)
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
    w1=np.array([-C12,C11-omega1])
    w2=np.array([-C12,C11-omega2])
    return [omega1,w1],[omega2,w2]

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

def computeDeviatoricPart(T):
    # T = [T11 T21 T22 T33]
    Pdev=np.array([[1.-1/3.,0.,-1./3.,-1./3.],[0.,1.,0.,0.],[-1./3.,0.,1.-1./3.,-1./3.],[-1./3.,0.,-1./3.,1.-1./3.]])
    Tdev=np.dot(Pdev,T)
    return np.array([Tdev[0],np.sqrt(2.)*Tdev[1],Tdev[2],Tdev[3]])

def computeCriterion(sig11,sig22,sig12,sig33,sigy):
    # deviatoric stress
    sDev=computeDeviatoricPart(np.array([sig11,sig12,sig22,sig33]))
    normSDev=np.sqrt(np.dot(sDev,sDev))
    f=np.sqrt(3./2.)*normSDev - sigy
    return f

def computePsiFast(sig12,sigma,sig33,lamb,mu,beta,tangent):
    # sig12 driven
    n1=1.;n2=0.
    sig11=sigma[0];sig22=sigma[1]
    H=tangentModulus(np.array([sig11,sig12,sig22,sig33]),lamb,mu,beta,tangent)
    C=acousticTensor(H,np.array([n1,n2]))
    eigenf,eigens=acousticEigenStructure(C)
    alpha11= (H[0,1]*n1+H[1,1]*n2)*(H[1,2]*n1+H[2,2]*n2)-(H[0,2]*n1+H[1,2]*n2)*(H[1,1]*n1+H[1,2]*n2)
    alpha12=((H[0,1]*n1+H[0,2]*n2)*(H[0,2]*n1+H[1,2]*n2)-(H[0,0]*n1+H[0,1]*n2)*(H[1,2]*n1+H[2,2]*n2))/2.
    alpha22= (H[0,0]*n1+H[0,1]*n2)*(H[1,1]*n1+H[1,2]*n2)-(H[0,1]*n1+H[0,2]*n2)*(H[0,1]*n1+H[1,1]*n2)
    w1=eigens[1][0];w2=eigens[1][1]
    psi11=-w2/(1.*w1)
    psi22=(w2*alpha11/(1.*w1)-alpha12)/alpha22
    return np.array([psi11,psi22])

def computePsiFast(sig12,sigma,sig33,lamb,mu,beta,tangent):
    # sig12 driven
    n1=1.;n2=0.
    sig11=sigma[0];sig22=sigma[1]
    H=tangentModulus(np.array([sig11,sig12,sig22,sig33]),lamb,mu,beta,tangent)
    C=acousticTensor(H,np.array([n1,n2]))
    eigenf,eigens=acousticEigenStructure(C)
    alpha11= (H[0,1]*n1+H[1,1]*n2)*(H[1,2]*n1+H[2,2]*n2)-(H[0,2]*n1+H[1,2]*n2)*(H[1,1]*n1+H[1,2]*n2)
    alpha12=((H[0,1]*n1+H[0,2]*n2)*(H[0,2]*n1+H[1,2]*n2)-(H[0,0]*n1+H[0,1]*n2)*(H[1,2]*n1+H[2,2]*n2))/2.
    alpha22= (H[0,0]*n1+H[0,1]*n2)*(H[1,1]*n1+H[1,2]*n2)-(H[0,1]*n1+H[0,2]*n2)*(H[0,1]*n1+H[1,1]*n2)
    w1=eigens[1][0];w2=eigens[1][1]
    psi11=-w2/(1.*w1)
    psi22=(w2*alpha11/(1.*w1)-alpha12)/alpha22
    return np.array([psi11,psi22])
def computePsiFast(sig12,sigma,sig33,lamb,mu,beta,tangent):
    # sig12 driven
    n1=1.;n2=0.
    sig11=sigma[0];sig22=sigma[1]
    H=tangentModulus(np.array([sig11,sig12,sig22,sig33]),lamb,mu,beta,tangent)
    C=acousticTensor(H,np.array([n1,n2]))
    eigenf,eigens=acousticEigenStructure(C)
    alpha11= (H[0,1]*n1+H[1,1]*n2)*(H[1,2]*n1+H[2,2]*n2)-(H[0,2]*n1+H[1,2]*n2)*(H[1,1]*n1+H[1,2]*n2)
    alpha12=((H[0,1]*n1+H[0,2]*n2)*(H[0,2]*n1+H[1,2]*n2)-(H[0,0]*n1+H[0,1]*n2)*(H[1,2]*n1+H[2,2]*n2))/2.
    alpha22= (H[0,0]*n1+H[0,1]*n2)*(H[1,1]*n1+H[1,2]*n2)-(H[0,1]*n1+H[0,2]*n2)*(H[0,1]*n1+H[1,1]*n2)
    w1=eigens[1][0];w2=eigens[1][1]
    psi11=-w2/(1.*w1)
    psi22=(w2*alpha11/(1.*w1)-alpha12)/alpha22
    return np.array([psi11,psi22])

def computeSpeed(sigma,lamb,mu,beta,tangent):
    # sig12 driven
    n1=1.;n2=0.
    sig11=sigma[0,0];sig22=sigma[1,1];sig12=sigma[0,1];sig33=sigma[2,2]
    H=tangentModulus(np.array([sig11,sig12,sig22,sig33]),lamb,mu,beta,tangent)
    C=acousticTensor(H,np.array([n1,n2]))
    eigenf,eigens=acousticEigenStructure(C)
    return eigens[0]


def loadingPathSig(dsig,sigma,sigman,E,nu,h,lamb,mu,beta,tangent):
    # sig11 driven
    n1=1.;n2=0.
    # Stress part
    sig11=sigman[0];sig12=sigman[1];sig22=sigman[2];sig33=sigman[3]
    H=tangentModulus(np.array([sig11,sig12,sig22,sig33]),lamb,mu,beta,tangent)
    C=acousticTensor(H,np.array([n1,n2]))
    eigenf,eigens=acousticEigenStructure(C)
    alpha11= (H[0,1]*n1+H[1,1]*n2)*(H[1,2]*n1+H[2,2]*n2)-(H[0,2]*n1+H[1,2]*n2)*(H[1,1]*n1+H[1,2]*n2)
    alpha12=((H[0,1]*n1+H[0,2]*n2)*(H[0,2]*n1+H[1,2]*n2)-(H[0,0]*n1+H[0,1]*n2)*(H[1,2]*n1+H[2,2]*n2))/2.
    alpha22= (H[0,0]*n1+H[0,1]*n2)*(H[1,1]*n1+H[1,2]*n2)-(H[0,1]*n1+H[0,2]*n2)*(H[0,1]*n1+H[1,1]*n2)
    w1=eigenf[1][0];w2=eigenf[1][1]
    psi12=-w1/w2
    psi22=(w1*alpha12/w2-alpha11)/alpha22
    ## Plastic strain part
    sigDev=computeDeviatoricPart(np.array([sigma[0],sigma[1],sigma[2],sigma[3]]))
    sigDevn=computeDeviatoricPart(np.array([sigman[0],sigman[1],sigman[2],sigman[3]]))
    norm=np.sqrt(np.dot(sigDevn,sigDevn))
    flow=sigDevn/norm
    dSig=sigman-sigma
    dSig[1]*=np.sqrt(2.)
    dp=(1./h)*np.sqrt(3./2.)*np.dot(flow,dSig)
    depsp33=np.sqrt(3./2.)*flow[3]*dp
    return np.array([psi12*dsig,psi22*dsig,nu*(1.+psi22)*dsig-E*depsp33,depsp33])

def integrateODESig(dsig,sig0,tau0,sig22_0,sig33,epsp33,nu,E,H,lamb,mu,beta,tangent):
    sigma=np.array([tau0,sig22_0,sig33,epsp33])
    # subdivision of time step
    sub_steps = 1
    dSIG = dsig/sub_steps
    theta = 1.0
    for i in range(sub_steps):
        ## Nonlinear solution procedure
        ### R = s^{n+1} - s^{n} - RHS
        R=lambda x : x - sigma - theta*loadingPathSig(dSIG,np.array([sig0,sigma[0],sigma[1],sigma[2]]),np.array([sig0+dSIG,x[0],x[1],x[2]]),E,nu,H,lamb,mu,beta,tangent)
        solution = scipy.optimize.root(R,np.array([sigma[0],sigma[1],sigma[2],sigma[3]])).x
        sigma = solution
    return solution[0],solution[1],solution[3],solution[2]

def loadingPath(dtau,sigma,sigman,E,nu,h,lamb,mu,beta,tangent):
    # Stress part
    n1=1.;n2=0.
    sig11=sigman[0];sig12=sigman[1];sig22=sigman[2];sig33=sigman[3]
    H=tangentModulus(np.array([sig11,sig12,sig22,sig33]),lamb,mu,beta,tangent)
    C=acousticTensor(H,np.array([n1,n2]))
    eigenf,eigens=acousticEigenStructure(C)
    alpha11= (H[0,1]*n1+H[1,1]*n2)*(H[1,2]*n1+H[2,2]*n2)-(H[0,2]*n1+H[1,2]*n2)*(H[1,1]*n1+H[1,2]*n2)
    alpha12=((H[0,1]*n1+H[0,2]*n2)*(H[0,2]*n1+H[1,2]*n2)-(H[0,0]*n1+H[0,1]*n2)*(H[1,2]*n1+H[2,2]*n2))/2.
    alpha22= (H[0,0]*n1+H[0,1]*n2)*(H[1,1]*n1+H[1,2]*n2)-(H[0,1]*n1+H[0,2]*n2)*(H[0,1]*n1+H[1,1]*n2)
    w1=eigenf[1][0];w2=eigenf[1][1]
    psi11=-w2/(1.*w1)
    psi22=(w2*alpha11/(1.*w1)-alpha12)/alpha22
    ## Plastic strain part
    sigDev=computeDeviatoricPart(np.array([sigma[0],sigma[1],sigma[2],sigma[3]]))
    sigDevn=computeDeviatoricPart(np.array([sigman[0],sigman[1],sigman[2],sigman[3]]))
    norm=np.sqrt(np.dot(sigDevn,sigDevn))
    flow=sigDevn/norm
    dSig=sigman-sigma
    dSig[1]*=np.sqrt(2.)
    #dSig=sigDevn-sigDev
    dp=(1./h)*np.sqrt(3./2.)*np.dot(flow,dSig)
    depsp33=np.sqrt(3./2.)*flow[3]*dp
    return np.array([psi11*dtau,psi22*dtau,nu*(psi11+psi22)*dtau-E*depsp33,depsp33])

def integrateODE(dtau,sig0,tau0,sig22_0,sig33,epsp33,nu,E,H,lamb,mu,beta,tangent):
    sigma=np.array([sig0,sig22_0,sig33,epsp33])
    # subdivision of time step
    sub_steps = 1
    dTAU = dtau/sub_steps
    theta = 1.0
    for i in range(sub_steps):
        ## Nonlinear solution procedure
        ## R = s^{n+1} - s^{n} - RHS
        R=lambda x: x - sigma - theta*loadingPath(dTAU,np.array([sigma[0],tau0,sigma[1],sigma[2]]),np.array([x[0],tau0+dTAU,x[1],x[2]]),E,nu,H,lamb,mu,beta,tangent)
        solution = scipy.optimize.root(R,np.array([sigma[0],sigma[1],sigma[2],sigma[3]])).x
        sigma = solution
    
    return solution[0],solution[1],solution[2],solution[3]

def computeLodeAngle(sig11,sig22,sig12,sig33):
    # deviatoric stress
    sDev=computeDeviatoricPart(np.array([sig11,sig12,sig22,sig33]))
    s11=sDev[0];s12=sDev[1]/np.sqrt(2.);s22=sDev[2];s33=sDev[3]
    sig=computeEigenStresses(np.matrix([[s11,s12,0.],[s12,s22,0.],[0.,0.,s33]]))
    # deviator 2nd and 3rd invariants
    J3=s33*(s11*s22-s12**2) ; sqrtJ2=np.sqrt(0.5*np.dot(sDev,sDev))
    # tan=(1./np.sqrt(3.))*(2.*(sig[1]-sig[2])/(sig[0]-sig[2])-1.)
    # theta=-np.sign(tan)*np.arccos((3./2.)*np.sqrt(3.)*J3/(sqrtJ2**3))/3.
    # theta=theta*360./(2.*np.pi)
    Seq = np.sqrt(np.dot(sDev,sDev))
    X=(3./2.)*np.sqrt(3.)*J3/(sqrtJ2**3)
    theta = np.arccos(X)/3.
    theta=theta*360./(2.*np.pi)
    return X

def updateEquivalentPlasticStrain(sig,sign,H):
    # sig=[sig11^n , sqrt(2)*sig12^n , sig22 , sig33^n]
    # sign=[sig11^n+1 , sqrt(2)*sig12^n+1 , sig22 , sig33^n+1]
    sigDevn=computeDeviatoricPart(np.array([sign[0],sign[1]/np.sqrt(2.),sign[2],sign[3]]))
    sigDev=computeDeviatoricPart(np.array([sig[0],sig[1]/np.sqrt(2.),sig[2],sig[3]]))
    norm=np.sqrt(np.dot(sigDevn,sigDevn))
    flow=sigDevn/norm
    dSig=sign-sig
    #dSig=sigDevn-sigDev
    dp=(1./H)*np.sqrt(3./2.)*np.dot(flow,dSig)
    return dp

def computeEigenStresses(sig):
    #    | sig11 sig12   0   |
    #sig=| sig12 sig22   0   |
    #    |   0     0   sig33 |
    s3=sig[2,2]
    delta=(sig[0,0]-sig[1,1])**2+4.*sig[0,1]**2
    s1=0.5*(sig[0,0]+sig[1,1]-np.sqrt(delta))
    s2=0.5*(sig[0,0]+sig[1,1]+np.sqrt(delta))
    return np.array([s1,s2,s3])

def computePlasticResidual(epsp33,sig,epsp33n,sign,E,H,nu):
    # sig = [sig11 , sig12*sqrt(2) , sig22 , sig33] (previous time step)
    # sig33n = nu*(sig11+sig22)-E*epsp33n (updated time step)
    sign[3]=nu*(sign[0]+sign[2])-E*epsp33n
    sigDev=computeDeviatoricPart(np.array([sig[0],sig[1]/np.sqrt(2.),sig[2],sig[3]]))
    sigDevn=computeDeviatoricPart(np.array([sign[0],sign[1]/np.sqrt(2.),sign[2],sign[3]]))
    norm=np.sqrt(np.dot(sigDevn,sigDevn))
    flow=sigDevn/norm
    #dSig=sign-sig
    dSig=sigDevn-sigDev
    dp=(1./H)*np.sqrt(3./2.)*np.dot(flow,dSig)
    res=epsp33n-epsp33-np.sqrt(3./2.)*flow[3]*dp
    return res


from mpl_toolkits.mplot3d import proj3d
def orthogonal_proj(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    return np.array([[1,0,0,0],
                        [0,1,0,0],
                        [0,0,a,b],
                        [0,0,0,zback]])
proj3d.persp_transformation = orthogonal_proj

Samples=10

## Setting initial values
# admissible range for sigma11
sig=np.linspace(-sigy*np.sqrt(4*(nu**2-nu+1.))/np.sqrt(3.*(4.*nu**2-4.*nu+1.)),sigy*np.sqrt(4*(nu**2-nu+1.))/np.sqrt(3.*(4.*nu**2-4.*nu+1.)),Samples)
# corresponding values of sig22 through elastic characteristics
sig22=sig*lamb/(2.*mu+lamb)
# associated value of sigma12
delta=((nu-nu**2)*(sig+sig22)**2 -sig**2-sig22**2+sig*sig22 + sigy**2)/3.
tau=np.sqrt(delta)


Niter=500
TAU=np.zeros((Niter,len(sig22)))
SIG11=np.zeros((Niter,len(sig22)))
SIG22=np.zeros((Niter,len(sig22)))
SIG33=np.zeros((Niter,len(sig22)))
eigsigS=np.zeros((Niter,len(sig22),3))
criterionF=np.zeros((Niter,len(sig22)))
PsiS=np.zeros((Samples,len(sig22)))

plast_F=np.zeros((Niter,len(sig22)))
rcs2=np.zeros((Niter,len(sig22)))
Epsp33=np.zeros((Niter,len(sig22)))
LodeAngle_F=np.zeros((Niter,len(sig22)))
sigma_b=np.zeros((Niter,len(sig22)))
radius_F=np.zeros((len(sig22)))
# Boolean to plot the upadted yield surface
updated_criterion=False

## LOADING PATHS PLOTS
pgfFilesList=[]
yields11_s12=[]
yields22_s12=[]
deviatorPlots=[]

fig = plt.figure()
ax3=plt.subplot2grid((1,1),(0,0),projection='3d')

cylindre=vonMisesYieldSurface(sigy)
ax3.plot(cylindre[0,:],cylindre[1,:],cylindre[2,:], color="k")
elevation_Angle_radian=np.arctan(1./np.sqrt(2.0))
angle_degree= 180.*elevation_Angle_radian/np.pi
radius=1.*np.sqrt((2./3.)*sigy**2)
ax3.set_xlim(-1.*radius,1.*radius)
ax3.set_ylim(-1.*radius,1.*radius)
ax3.set_zlim(-1.*radius,1.*radius)
ax3.view_init(angle_degree,45.)
ax3.plot([0.,sigy],[0.,sigy],[0.,sigy],color="k")
ax3.set_xlabel(r'$\sigma_1$',size=24.)
ax3.set_ylabel(r'$\sigma_2$',size=24.)
ax3.set_zlabel(r'$\sigma_3$',size=24.)
    
fig = plt.figure()

ax1=plt.subplot2grid((1,2),(0,0))
ax2=plt.subplot2grid((1,2),(0,1))
ax1.set_xlabel(r'$\sigma_{11}$')
#ax1.set_ylabel(r'$\sigma_{12}$')
ax2.set_xlabel(r'$\sigma_{22}$')
#ax2.set_ylabel(r'$\sigma_{12}$')
ax1.grid()
ax2.grid()

fig = plt.figure()
ax4=plt.subplot2grid((1,1),(0,0),projection='3d')
ax4.set_xlabel(r'$\sigma_{11}$',size=24.)
ax4.set_ylabel(r'$\sigma_{22}$',size=24.)
ax4.set_zlabel(r'$\sigma_{12}$',size=24.)
ax4.grid()

#### PLOTTING THE INITIAL YIELD SURFACE IN STRESS SPACE
samp=25

# sig=np.linspace(-sigy*np.sqrt(4*(nu**2-nu+1.))/np.sqrt(3.*(4.*nu**2-4.*nu+1.)),sigy*np.sqrt(4*(nu**2-nu+1.))/np.sqrt(3.*(4.*nu**2-4.*nu+1.)),Samples)
# # corresponding values of sig22 through elastic characteristics
# sig22=sig*lamb/(2.*mu+lamb)
# # associated value of sigma12
# delta=((nu-nu**2)*(sig+sig22)**2 -sig**2-sig22**2+sig*sig22 + sigy**2)/3.
# tau=np.sqrt(delta)

# Sample constant stress component sig22
sigma2=np.linspace(-sigy*np.sqrt(4*(nu**2-nu+1.))/np.sqrt(3.*(4.*nu**2-4.*nu+1.)),sigy*np.sqrt(4*(nu**2-nu+1.))/np.sqrt(3.*(4.*nu**2-4.*nu+1.)),samp)
sigma1=np.zeros((50,50))
sigma3=np.zeros((50,50))
sigma_lim=np.zeros(25)
tau_lim=np.zeros(25)
TeXFile=open("yieldSurfaceSingularity.tex","w")
TeXFile.write(r'\begin{tikzpicture}');TeXFile.write('\n')
TeXFile.write(r'\begin{axis}[width=.75\textwidth,view={135}{35.2643},xlabel=$\sigma_{11}$,ylabel=$\sigma_{22}$,zlabel=$\sigma_{12}$,every axis y label/.style={at={(rel axis cs:0.,.5,-0.65)}, anchor=west}, every axis x label/.style={at={(rel axis cs:0.5,.,-0.65)}, anchor=east}, every axis z label/.style={at={(rel axis cs:0.,.0,.18)}, anchor=north}]');TeXFile.write('\n')

for k in range(len(sigma2)-1)[1:]:
    s22=sigma2[k]

    Delta=(4.*(nu**2-nu+1.)*sigy**2- 3.*(4.*nu**2-4.*nu+1.)*s22**2)
    sigMax=(s22*(1.+2.*nu-2.*nu**2)+np.sqrt(Delta))/(2.*(nu**2-nu+1.))
    sigMin=(s22*(1.+2.*nu-2.*nu**2)-np.sqrt(Delta))/(2.*(nu**2-nu+1.))

    sigma_lim[k]= 0.5*s22*(2.*nu**2-2.*nu-1.)/(nu-nu**2-1.)
    tau_lim[k]= np.sqrt(((nu-nu**2)*(sigma_lim[k]+s22)**2 -sigma_lim[k]**2-s22**2+sigma_lim[k]*s22 + sigy**2)/3.)
    # Sample stress component sig11
    sigma1[:,k]=np.linspace(sigMin,sigMax,50)
    TeXFile.write(r'\addplot3+[gray,dashed,thin,no markers] coordinates {')
                    
    # Compute shear stress satisfying the criterion given sig11 and sig22
    for i in range(50):
        s11=sigma1[i,k]
        delta=((nu-nu**2)*(s11+s22)**2 -s11**2-s22**2+s11*s22 + sigy**2)/3.
        if np.abs(delta)<10. : delta=np.abs(delta)
        sigma3[i,k]=np.sqrt(delta)
        TeXFile.write('('+str(s22)+','+str(s11)+','+str(sigma3[i,k])+') ')
    TeXFile.write('};\n')
    ax4.plot(np.ones(50)*sigma2[k],sigma1[:,k],sigma3[:,k],'k--',lw=0.5)


TeXFile.write(r'\addplot3+[black,very thick,no markers] coordinates {')
for k in range(len(sigma2)-1)[1:]:
    TeXFile.write('('+str(sigma2[k])+','+str(sigma_lim[k])+','+str(tau_lim[k])+') ')
TeXFile.write('};\n')
TeXFile.write(r'\end{axis}')
TeXFile.write('\n')
TeXFile.write('\end{tikzpicture}')
TeXFile.write('\n')
TeXFile.write('%%% Local Variables:')
TeXFile.write('\n')
TeXFile.write('%%% mode: latex')
TeXFile.write('\n')
TeXFile.write('%%% TeX-master: "../manuscript"')
TeXFile.write('\n')
TeXFile.write('%%% End:')
TeXFile.write('\n')
TeXFile.close()
              
sigma1= 0.5*sigma2*(2.*nu**2-2.*nu-1.)/(nu-nu**2-1.)
maX= np.sqrt(((nu-nu**2)*(sigma1+sigma2)**2 -sigma1**2-sigma2**2+sigma1*sigma2 + sigy**2)/3.)
ax4.plot(sigma2,sigma1,maX,'r')

##########################################################

tangent='planeStrain'
## LOADING PATHS PLOTS

chain_waves=True
chaining_parameter='sig'
bound = 2
for k in range(len(sig22))[bound:-bound]:
    #for k in range(len(sig22)):
    s22=sig22[k]
    sig0=sig[k]
    tau0=tau[k]

    print "sig22=",s22,", sig11=",sig0,", tau=",tau0
    print "initial yield function ", computeCriterion(sig0,s22,tau0,nu*(sig0+s22),sigy)

    tauM = 1.05*np.max(tau[bound:-bound])
    if H == 100.0e8 :     tauM = 50.*np.max(tau[bound:-bound])
    dtau=(tauM-tau0)/Niter
    TAU[:,k]=np.linspace(tau0,tauM,Niter)
    SIG11[0,k]=sig0
    
    SIG22[0,k]=s22
    
    sig33=nu*(SIG11[0,k]+SIG22[0,k])
    SIG33[0,k]=sig33


    sigma = np.matrix([[SIG11[0,k],TAU[0,k],0.],[TAU[0,k],SIG22[0,k],0.],[0.,0.,sig33]])
    rcs2[0,k]=np.sqrt(computeSpeed(sigma,lamb,mu,beta,tangent)/rho)
            
    sigDev=computeDeviatoricPart(np.array([SIG11[0,k],TAU[0,k],SIG22[0,k],SIG33[0,k]]))
    sigma = np.matrix([[sigDev[0],sigDev[1]/np.sqrt(2.),0.],[sigDev[1]/np.sqrt(2.),sigDev[2],0.],[0.,0.,sigDev[3]]])
            
    eigsigS[0,k,:]=computeEigenStresses(sigma)
    LodeAngle_F[0,k]=computeLodeAngle(sigma[0,0],SIG22[0,k],sigma[0,1],sig33)
    LodeAngle_F[0,k]=computeLodeAngle(sig0,s22,tau0,sig33)
    sigma_b[0,k]=SIG11[0,k]*(1.+2.*nu-2.*nu**2)/(2.*nu**2-2.*nu+2.)
    
    plast=0.
    epsp33=0.
    for j in range(Niter-1):
        SIG11[j+1,k],SIG22[j+1,k],SIG33[j+1,k],epsp33=integrateODE(dtau,SIG11[j,k],TAU[j,k],SIG22[j,k],SIG33[j,k],epsp33,nu,E,H,lamb,mu,beta,tangent)
        
        sigma = np.array([SIG11[j,k],np.sqrt(2.)*TAU[j,k],SIG22[j,k],SIG33[j,k]])
        sigman=np.array([SIG11[j+1,k],TAU[j+1,k]*np.sqrt(2.),SIG22[j+1,k],SIG33[j+1,k]])
        # Plastic update
        dp=updateEquivalentPlasticStrain(sigma,sigman,H)
        plast+=dp
        sig33=SIG33[j+1,k]
        
        criterionF[j+1,k]=computeCriterion(SIG11[j+1,k],SIG22[j+1,k],TAU[j+1,k],sig33,sigy+H*plast)
        plast_F[j+1,k]=plast
        Epsp33[j+1,k]=epsp33
        LodeAngle_F[j+1,k]=computeLodeAngle(sigman[0],sigman[2],sigman[1]/np.sqrt(2.),sig33)
        sigma_b[j+1,k]=E*epsp33*(2.*nu-1)/(2.*nu**2-2.*nu+2.) +SIG11[j+1,k]*(1.+2.*nu-2.*nu**2)/(2.*nu**2-2.*nu+2.)
        
        # Eigenvalues of sigma (for deviatoric plane plots)
        sigma = np.matrix([[SIG11[j+1,k],TAU[j+1,k],0.],[TAU[j+1,k],SIG22[j+1,k],0.],[0.,0.,SIG33[j+1,k]]])
        rcs2[j+1,k]=np.sqrt(computeSpeed(sigma,lamb,mu,beta,tangent)/rho)
        
        # sigDev=computeDeviatoricPart(np.array([SIG11[j+1,k],TAU[j+1,k],SIG22[j+1,k],SIG33[j+1,k]]))
        # sigma = np.matrix([[sigDev[0],sigDev[1]/np.sqrt(2.),0.],[sigDev[1]/np.sqrt(2.),sigDev[2],0.],[0.,0.,sigDev[3]]])
        eigsigS[j+1,k,:]=computeEigenStresses(sigma)
        rcs2[j+1,k]=np.sqrt(computeSpeed(sigma,lamb,mu,beta,tangent)/rho)
        Nfinal=Niter
        if rcs2[j+1,k]>rcs2[j,k]:
            print "Simple wave condition violated. Determinant of sigma=",eigsigS[j,k,0]*eigsigS[j,k,1]*eigsigS[j,k,2]
            Nfinal=j
            break

    if chain_waves and Nfinal!=Niter:
        print "Trying to change the driving parameter (",chaining_parameter,")"
        sig0=SIG11[Nfinal,k]
        sigM=1.05*sig0
        dsig=(sigM-sig0)/(Niter-Nfinal-1)
        SIG11[Nfinal:,k]=np.linspace(sig0,sigM,Niter-Nfinal)

        for j in range(Niter-1)[Nfinal:]:
            TAU[j+1,k],SIG22[j+1,k],SIG33[j+1,k],epsp33=integrateODESig(dsig,SIG11[j,k],TAU[j,k],SIG22[j,k],SIG33[j,k],epsp33,nu,E,H,lamb,mu,beta,tangent)

            sigma = np.array([SIG11[j,k],np.sqrt(2.)*TAU[j,k],SIG22[j,k],SIG33[j,k]])
            sigman=np.array([SIG11[j+1,k],TAU[j+1,k]*np.sqrt(2.),SIG22[j+1,k],SIG33[j+1,k]])
            # Plastic update
            dp=updateEquivalentPlasticStrain(sigma,sigman,H)
            plast+=dp
            sig33=SIG33[j+1,k]
            
            criterionF[j+1,k]=computeCriterion(SIG11[j+1,k],SIG22[j+1,k],TAU[j+1,k],sig33,sigy+H*plast)
            plast_F[j+1,k]=plast
            Epsp33[j+1,k]=epsp33
            LodeAngle_F[j+1,k]=computeLodeAngle(sigman[0],sigman[2],sigman[1]/np.sqrt(2.),sig33)
            sigma_b[j+1,k]=E*epsp33*(2.*nu-1)/(2.*nu**2-2.*nu+2.) +SIG11[j+1,k]*(1.+2.*nu-2.*nu**2)/(2.*nu**2-2.*nu+2.)

            # Eigenvalues of sigma (for deviatoric plane plots)
            sigma = np.matrix([[SIG11[j+1,k],TAU[j+1,k],0.],[TAU[j+1,k],SIG22[j+1,k],0.],[0.,0.,SIG33[j+1,k]]])
            rcs2[j+1,k]=np.sqrt(computeSpeed(sigma,lamb,mu,beta,tangent)/rho)
            
            # sigDev=computeDeviatoricPart(np.array([SIG11[j+1,k],TAU[j+1,k],SIG22[j+1,k],SIG33[j+1,k]]))
            # sigma = np.matrix([[sigDev[0],sigDev[1]/np.sqrt(2.),0.],[sigDev[1]/np.sqrt(2.),sigDev[2],0.],[0.,0.,sigDev[3]]])
            eigsigS[j+1,k,:]=computeEigenStresses(sigma)
            rcs2[j+1,k]=np.sqrt(computeSpeed(sigma,lamb,mu,beta,tangent)/rho)
            Nfinal=Niter
            if rcs2[j+1,k]>rcs2[j,k]:
                print "Simple wave condition violated"
                Nfinal=j
                break

    print "Final equivalent plastic strain after slow wave : ",plast
    if Nfinal/100>=1:
        pas=Nfinal/100
        ranging=np.arange(0,Nfinal-1,pas)
    else:
        pas=1
        ranging=np.linspace(0,Nfinal-1,Nfinal,True,False,np.int)
    fileName=path+'DPslowStressPlane_Stress'+str(k)+'.pgf'
    pgfFilesList.append(fileName)
    export2pgfPlotFile(fileName,np.array([TAU[ranging,k],SIG11[ranging,k],SIG22[ranging,k],rcs2[ranging,k],LodeAngle_F[ranging,k],sigma_b[ranging,k],np.sqrt(eigsigS[ranging,k,0]**2 + eigsigS[ranging,k,1]**2 + eigsigS[ranging,k,2]**2)]),'sigma_12','sigma_11','sigma_22','p','Theta','sigmab','radius')
    fileName=path+'DPslowDevPlane_Stress'+str(k)+'.pgf'
    deviatorPlots.append(fileName)
    dico={"xlabel":r'$\sigma_1$',"ylabel":r'$\sigma_2$',"zlabel":r'$\sigma_3$'}
    export2pgfPlot3D(fileName,eigsigS[ranging,k,0],eigsigS[ranging,k,1],eigsigS[ranging,k,2],dico)
    
    fileName=path+'slowWave_DevPlane_StressValue'+str(k)+'.pgf'
    dico={"xlabel":r'$\sigma_1$',"ylabel":r'$\sigma_2$',"zlabel":r'$\sigma_3$'}
        

    fileName=path+'CylindreDevPlane.pgf'
    export2pgfPlot3D(str(fileName),cylindre[0,:],cylindre[1,:],cylindre[2,:])
    deviatorPlots.append(fileName)

    ### SUBPLOTS SETTINGS
    ax3.plot(eigsigS[ranging,k,0],eigsigS[ranging,k,1],eigsigS[ranging,k,2],lw=1.5)
    ax3.plot([-sigy,sigy],[0.,0.],[0.,0.],color="k",linestyle="--",lw=1.)
    ax3.plot([0.,0.],[-sigy,sigy],[0.,0.],color="k",linestyle="--",lw=1.)
    ax3.plot([-radius,radius],[radius,-radius],[0.,0.],color="k",linestyle="--",lw=1.)
    plt.tight_layout()
    #plt.show()

    ax4.plot(SIG11[ranging,k],SIG22[ranging,k],TAU[ranging,k])
    plt.tight_layout()
    
    # ax1.plot(SIG11[ranging,k],TAU[ranging,k])
    # #ax2.plot(SIG22[:,k],TAU[:,k])
    # ax2.plot(LodeAngle_F[ranging,k],SIG22[ranging,k]/sigma_b[ranging,k])

    # Look at triaxiality and invariants of sigma
    
    ax1.set_xlabel(r'$p$',size=28.)
    ax2.set_xlabel(r'$\eta$',size=28.)
    norm_vm = criterionF[:,k]+sigy+H*plast_F[:,k]
    press = (1./3.)*(eigsigS[:,k,0]+eigsigS[:,k,1]+eigsigS[:,k,2])
    triax = press/norm_vm
    det=eigsigS[:,k,0]*eigsigS[:,k,1]*eigsigS[:,k,2]
    # ax1.plot(press,2.*SIG22[:,k]/SIG11[:,k])
    ax1.plot(press[ranging]/max(abs(det)),SIG22[ranging,k]/sigma_b[ranging,k],linestyle='-.')
    #ax1.plot(np.arange(0,len(det),1),det,linestyle='-.')
    ax2.plot(triax[ranging]/max(abs(triax)),SIG22[ranging,k]/sigma_b[ranging,k],linestyle='-.')
    #ax2.plot(np.arange(0,len(triax),1),triax,linestyle='-.')
    #ax4.plot(SIG22[:,k],LodeAngle_F[:,k],linestyle='--')
    

    # Look at the principal stress with respect to the maximum-shear condition
    """
    ax1.set_xlabel(r'$\sigma_1$',size=28.)
    ax2.set_xlabel(r'$\sigma_2$',size=28.)
    ax1.plot(eigsigS[ranging,k,0],SIG22[ranging,k]/sigma_b[ranging,k],linestyle='-.')
    ax2.plot(eigsigS[ranging,k,2],SIG22[ranging,k]/sigma_b[ranging,k],linestyle='-.')
    """
    
    ## sig22 value will change here
    subtitle=[r'(a) ($\sigma_{11},\sigma_{12}$) plane',r'(b) ($\sigma_{22},\sigma_{12}$) plane',r'(c) Deviatoric plane']

    srcX=['sigma_11','sigma_22']
    srcY=['sigma_12','sigma_12']

    name1='DPslowWaves_sig11_tau'+str(k)+'.tex'
    name2='DPslowWaves_sig22_tau'+str(k)+'.tex'
    name3='DPslowWaves_deviator'+str(k)+'.tex'
    names=[[name1,name2],name3]
    
    files1=np.concatenate([pgfFilesList,yields11_s12])
    files2=pgfFilesList
    pgfFiles=[[files1,files2],deviatorPlots]
    xlabels=[['$\sigma_{11} (Pa)$','$\sigma_{22}  (Pa)$'],'$s_1 $'] #size=number of .tex files
    ylabels=[['$\sigma_{12}  (Pa)$','$\sigma_{12}  (Pa)$'],'$s_2 $'] #size=number of .tex files
    zlabels=[['',''],'$s_3$'] #size=number of .tex files
    
    TauMax=1.1*np.max(TAU[0:-1:Niter/100,k])
    #buildTeXFiles2(names,pgfFiles,xlabels,ylabels,zlabels,srcX,srcY,TauMax)
    
    pgfFilesList=[];yields11_s12=[];
plt.show()
