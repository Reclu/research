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
    w1=np.array([C22-omega1,-C12])
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
"""
def loadingPath(dsig,sigma,sigman,h,lamb,mu,beta,tangent):
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
    w1=eigens[1][0];w2=eigens[1][1]
    psi12=-w1/(1.*w2)
    psi22=(w1*alpha12/(1.*w2)-alpha11)/alpha22
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
    # return dtau,dsig22,dsig33,depsp33
    return np.array([psi12*dsig,psi22*dsig,depsp33])

def integrateODE(dsig,sig0,tau0,sig22_0,sig33,epsp33,nu,E,H,lamb,mu,beta,tangent):
    sigma=np.array([tau0,sig22_0,epsp33])
    # subdivision of time step
    sub_steps = 1
    dSIG = dsig/sub_steps
    theta = 1.0
    for i in range(sub_steps):
        ## Nonlinear solution procedure
        ## R = s^{n+1} - s^{n} - RHS
        #R=lambda x: x - sigma - theta*np.array([dSIG*computePsiFast(sig0+dSIG,np.array([x[0],x[1]]),nu*(sig0+dSIG+x[1])-E*x[2],lamb,mu,beta,tangent)[0],dSIG*computePsiFast(sig0+dSIG,np.array([x[0],x[1]]),nu*(sig0+dSIG+x[1])-E*x[2],lamb,mu,beta,tangent)[1],computePlasticResidual2(np.array([sig0,sigma[0]*np.sqrt(2.),sigma[1],nu*(sig0+sigma[1])-E*sigma[2]]),np.array([sig0+dSIG,x[0]*np.sqrt(2.),x[1],nu*(sig0+dSIG+x[1])-E*x[2]]),H)])
        R=lambda x : x - sigma - theta*loadingPath(dSIG,np.array([sig0,sigma[0],sigma[1],nu*(sig0+sigma[1])-E*sigma[2]]),np.array([sig0+dSIG,x[0],x[1],nu*(sig0+dSIG+x[1])-E*x[2]]),H,lamb,mu,beta,tangent)
        solution = scipy.optimize.fsolve(R,np.array([sigma[0],sigma[1],sigma[2]]))
        
        sigma = solution
    return solution[0],solution[1],solution[2],nu*(dsig+sig0+solution[1])-E*solution[2]

"""
def loadingPath(dsig,sigma,sigman,E,nu,h,lamb,mu,beta,tangent):
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
    w1=eigens[1][0];w2=eigens[1][1]
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
    # return dtau,dsig22,dsig33,depsp33
    ## Singularities of the problem
    #root1: sig12=0
    #root2: beta*s11**2=lamb+mu
    #root3: s11=0
    # if sig12==0.:
    #     if abs(beta*sigDev[0]**2-lamb-mu)==0.:
    #         print "Loss of hyperbolicity"
    #         #pdb.set_trace()
    #     elif  abs(beta*sigDev[0]**2-lamb-2.*mu)<=0.:
    #         print "Loss of hyperbolicity"
    #         #pdb.set_trace()
    #     else:
    #         print "Horizontal loading path"
    #         #pdb.set_trace()
    # elif sigDev[0]==0.:
    #     print "Horizontal loading path"
    #     #pdb.set_trace()
    return np.array([psi12*dsig,psi22*dsig,nu*(1.+psi22)*dsig-E*depsp33,depsp33])

def integrateODE(dsig,sig0,tau0,sig22_0,sig33,epsp33,nu,E,H,lamb,mu,beta,tangent):
    sigma=np.array([tau0,sig22_0,sig33,epsp33])
    # subdivision of time step
    sub_steps = 1
    dSIG = dsig/sub_steps
    theta = 1.0
    for i in range(sub_steps):
        ## Nonlinear solution procedure
        ### R = s^{n+1} - s^{n} - RHS
        R=lambda x : x - sigma - theta*loadingPath(dSIG,np.array([sig0,sigma[0],sigma[1],sigma[2]]),np.array([sig0+dSIG,x[0],x[1],x[2]]),E,nu,H,lamb,mu,beta,tangent)
        solution = scipy.optimize.root(R,np.array([sigma[0],sigma[1],sigma[2],sigma[3]])).x
        sigma = solution
    return solution[0],solution[1],solution[3],solution[2]

def computeSpeed(sigma,lamb,mu,beta,tangent):
    # sig12 driven
    n1=1.;n2=0.
    sig11=sigma[0,0];sig22=sigma[1,1];sig12=sigma[0,1];sig33=sigma[2,2]
    H=tangentModulus(np.array([sig11,sig12,sig22,sig33]),lamb,mu,beta,tangent)
    C=acousticTensor(H,np.array([n1,n2]))
    eigenf,eigens=acousticEigenStructure(C)
    return eigenf[0]


def computeLodeAngle(sig11,sig22,sig12,sig33):
    # deviatoric stress
    sDev=computeDeviatoricPart(np.array([sig11,sig12,sig22,sig33]))
    s11=sDev[0];s12=sDev[1]/np.sqrt(2.);s22=sDev[2];s33=sDev[3]
    sig=computeEigenStresses(np.matrix([[s11,s12,0.],[s12,s22,0.],[0.,0.,s33]]))
    # deviator 2nd and 3rd invariants
    J3=s33*(s11*s22-s12**2) ; sqrtJ2=np.sqrt(0.5*np.dot(sDev,sDev))
    tan=(1./np.sqrt(3.))*(2.*(sig[1]-sig[2])/(sig[0]-sig[2])-1.)
    theta=-np.sign(tan)*np.arccos((3./2.)*np.sqrt(3.)*J3/(sqrtJ2**3))/3.
    theta=theta*360./(2.*np.pi)
    return theta

def updateEquivalentPlasticStrain(sig,sign,H):
    # sig=[sig11^n , sqrt(2)*sig12^n , sig22^n , sig33^n]
    # sign=[sig11^n+1 , sqrt(2)*sig12^n+1 , sig22^n+1 , sig33^n+1]
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


def updateEpsilon(Eps,sigma,sigman,E,nu,h):
    sigDevn=computeDeviatoricPart(np.array([sigman[0],sigman[1],sigman[2],sigman[3]]))
    norm=np.sqrt(np.dot(sigDevn,sigDevn))
    flow=sigDevn/norm
    dSig=sigman-sigma
    sigmaDot=np.array([dSig[0],dSig[1]*np.sqrt(2.),dSig[2],dSig[3]])
    deps11=(1.+nu)*dSig[0]/E-nu*(dSig[0]+dSig[2]+dSig[3])/E +(3./(2.*h))*flow[0]*np.dot(flow,sigmaDot)
    deps22=(1.+nu)*dSig[2]/E-nu*(dSig[0]+dSig[2]+dSig[3])/E +(3./(2.*h))*flow[2]*np.dot(flow,sigmaDot)
    deps33=(1.+nu)*dSig[3]/E-nu*(dSig[0]+dSig[2]+dSig[3])/E +(3./(2.*h))*flow[3]*np.dot(flow,sigmaDot)
    deps12=(1.+nu)*dSig[1]/E +(3./(2.*h))*flow[1]*np.dot(flow,sigmaDot)/np.sqrt(2.)
    Eps+=np.array([deps11,deps12,deps22,deps33])
    return Eps

def updatePlasticStrain(sigma,sigman,h):
    sigDev=computeDeviatoricPart(np.array([sigma[0],sigma[1],sigma[2],sigma[3]]))
    sigDevn=computeDeviatoricPart(np.array([sigman[0],sigman[1],sigman[2],sigman[3]]))
    norm=np.sqrt(np.dot(sigDevn,sigDevn))
    flow=sigDevn/norm
    dSig=sigman-sigma
    dSig[1]*=np.sqrt(2.)
    #dSig=sigDevn-sigDev
    dp=(1./h)*np.sqrt(3./2.)*np.dot(flow,dSig)
    dEpsp=np.sqrt(3./2.)*flow*dp
    return dEpsp

from mpl_toolkits.mplot3d import proj3d
def orthogonal_proj(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    return np.array([[1,0,0,0],
                        [0,1,0,0],
                        [0,0,a,b],
                        [0,0,0,zback]])
proj3d.persp_transformation = orthogonal_proj

Samples=6

# Sample constant stress component sig22
sig22=np.linspace(0.,sigy*np.sqrt(4*(nu**2-nu+1.))/np.sqrt(3.*(4.*nu**2-4.*nu+1.)),Samples)
sig22=np.linspace(-sigy*np.sqrt(4*(nu**2-nu+1.))/np.sqrt(3.*(4.*nu**2-4.*nu+1.)),sigy*np.sqrt(4*(nu**2-nu+1.))/np.sqrt(3.*(4.*nu**2-4.*nu+1.)),Samples)

sig22=[0.]

Samples*=10
sig=np.zeros((Samples,Samples))
tau=np.zeros((Samples,Samples))


frames=[0,10,29,31,49,Samples-1]
frames=[0]

col=["r","g","b","y","c","m","k","p"]

# purple to red
col=['#781C81','#3F60AE','#539EB6','#6DB388','#CAB843','#E78532','#D92120']
tauM=1.5*sigy/np.sqrt(3.)
sigM=5.*sigy/np.sqrt(1-nu+nu**2)
tauM=sigM
Niter=5000
TAU=np.zeros((Niter,len(frames),len(sig22)))
SIG11=np.zeros((Niter,len(frames),len(sig22)))
SIG22=np.zeros((Niter,len(frames),len(sig22)))
SIG33=np.zeros((Niter,len(frames),len(sig22)))
eigsigS=np.zeros((Niter,len(frames),len(sig22),3))
criterionF=np.zeros((Niter,len(frames),len(sig22)))
PsiS=np.zeros((Samples,len(sig22)))

plast_F=np.zeros((Niter,len(frames),len(sig22)))
Epsp33=np.zeros((Niter,len(frames),len(sig22)))
Eps=np.zeros((4,Niter,len(frames),len(sig22)))
LodeAngle_F=np.zeros((Niter,len(frames),len(sig22)))
radius_F=np.zeros((len(frames),len(sig22)))
rcf2=np.zeros((Niter,len(frames),len(sig22)))
# Boolean to plot the upadted yield surface
updated_criterion=True

## LOADING PATHS PLOTS
pgfFilesList=[]
yields11_s12=[]
yields22_s12=[]
deviatorPlots=[]

#for k in range(len(sig22)-1)[1:]:
for k in range(len(sig22)):
    s22=sig22[k]
    
    Delta=(4.*(nu**2-nu+1.)*sigy**2- 3.*(4.*nu**2-4.*nu+1.)*s22**2)
    sigMax=(s22*(1.+2.*nu-2.*nu**2)+np.sqrt(Delta))/(2.*(nu**2-nu+1.))
    sigMin=(s22*(1.+2.*nu-2.*nu**2)-np.sqrt(Delta))/(2.*(nu**2-nu+1.))
    
    # Sample stress component sig11
    sig[:,k]=np.linspace(sigMin,sigMax,Samples)
    #sig[:,k]=np.linspace(0.,sigMax,Samples)
    
    # Compute shear stress satisfying the criterion given sig11 and sig22
    for i in range(Samples):
        s11=sig[i,k]
        delta=((nu-nu**2)*(s11+s22)**2 -s11**2-s22**2+s11*s22 + sigy**2)/3.
        if np.abs(delta)<10. : delta=np.abs(delta)
        tau[i,k]=np.sqrt(delta)
        if np.isnan(tau[i,k]):
            print "Nan ",delta,s11,s22

exportPgf = True
tangent='planeStrain'
## LOADING PATHS PLOTS
for k in range(len(sig22)):
    s22=sig22[k]
    
    approx=np.zeros((len(frames),2))
    ordonnees=np.zeros((len(frames),Samples))
    abscisses=np.zeros((len(frames),Samples))
    for s,i in enumerate(frames):
        #pdb.set_trace()
        
        # if i==Samples-1:
        #     sig0=-sig[-1-i+1,k]
        #     tau0=tau[-1-i+1,k]
        #     sigMax=-1.*sigM
        # elif i==Samples-2:
        #     sig0=sig[-1-i+1,k]
        #     tau0=tau[-1-i+1,k]
        #     sigMax=1.*sigM
        sig0=sig[-1-i,k]
        tau0=tau[-1-i,k]
        
        if i>Samples/2:
            sigMax=-2.*sigM
        else:
            sigMax=5.*sigM
        
        # maxCrit=0.5*(s22*(2.*nu**2-2.*nu-1.))/(nu-nu**2-1.)
        # if sig0<maxCrit :
        #     sigMax=-2.*sigM
        # else:
        #     sigMax=2.*sigM
        # print "Maximum stress ",sigMax
        
        dsig=(sigMax-sig0)/Niter
        
        SIG11[:,s,k]=np.linspace(sig0,sigMax,Niter)
        
        TAU[0,s,k]=tau0
        SIG22[0,s,k]=s22
        
        sig33=nu*(SIG11[0,s,k]+SIG22[0,s,k])
        SIG33[0,s,k]=sig33
        
        
        sigma = np.matrix([[SIG11[0,s,k],TAU[0,s,k],0.],[TAU[0,s,k],SIG22[0,s,k],0.],[0.,0.,sig33]])

        rcf2[0,s,k] = np.sqrt(computeSpeed(sigma,lamb,mu,beta,tangent)/rho)
        
        # sigDev=computeDeviatoricPart(np.array([SIG11[0,s,k],TAU[0,s,k],SIG22[0,s,k],SIG33[0,s,k]]))
        # sigma = np.matrix([[sigDev[0],sigDev[1]/np.sqrt(2.),0.],[sigDev[1]/np.sqrt(2.),sigDev[2],0.],[0.,0.,sigDev[3]]])
            
        eigsigS[0,s,k,:]=computeEigenStresses(sigma)
        LodeAngle_F[0,s,k]=computeLodeAngle(sigma[0,0],SIG22[0,s,k],sigma[0,1],sig33)
            
        plast=0.
        epsp33=0.
        epsp=np.zeros(4)
        for j in range(Niter-1):
            
            TAU[j+1,s,k],SIG22[j+1,s,k],epsp33,SIG33[j+1,s,k]=integrateODE(dsig,SIG11[j,s,k],TAU[j,s,k],SIG22[j,s,k],SIG33[j,s,k],epsp33,nu,E,H,lamb,mu,beta,tangent)
            sigma = np.array([SIG11[j,s,k],np.sqrt(2.)*TAU[j,s,k],SIG22[j,s,k],SIG33[j,s,k]])
            # Plastic update
            sig33=SIG33[j+1,s,k]
            sigman=np.array([SIG11[j+1,s,k],TAU[j+1,s,k]*np.sqrt(2.),SIG22[j+1,s,k],sig33])
            Eps[:,j+1,s,k]=updateEpsilon(Eps[:,j,s,k],np.array([SIG11[j,s,k],TAU[j,s,k],SIG22[j,s,k],SIG33[j,s,k]]),np.array([SIG11[j+1,s,k],TAU[j+1,s,k],SIG22[j+1,s,k],SIG33[j+1,s,k]]),E,nu,H)
            dp=updateEquivalentPlasticStrain(sigma,sigman,H)
            # Eigenvalues of sigma (for deviatoric plane plots)
            sigma = np.matrix([[SIG11[j+1,s,k],TAU[j+1,s,k],0.],[TAU[j+1,s,k],SIG22[j+1,s,k],0.],[0.,0.,SIG33[j+1,s,k]]])
            rcf2[j+1,s,k] = np.sqrt(computeSpeed(sigma,lamb,mu,beta,tangent)/rho)
            # if dp<-5.e-6 :
            #     print "equivalent plastic strain increment negative. Increment=",dp
            #     pdb.set_trace()
            plast+=dp
            
            criterionF[j+1,s,k]=computeCriterion(SIG11[j+1,s,k],SIG22[j+1,s,k],TAU[j+1,s,k],sig33,sigy+H*plast)
            plast_F[j+1,s,k]=plast
            Epsp33[j+1,s,k]=epsp33
            LodeAngle_F[j+1,s,k]=computeLodeAngle(sigman[0],sigman[2],sigman[1]/np.sqrt(2.),sig33)
            epsp+=updatePlasticStrain(np.array([SIG11[j,s,k],TAU[j,s,k],SIG22[j,s,k],SIG33[j,s,k]]),np.array([SIG11[j+1,s,k],TAU[j+1,s,k],SIG22[j+1,s,k],SIG33[j+1,s,k]]),H)
            tr=epsp[0]+epsp[2]+epsp[3]
            if abs(tr) >1.e-15:
                print "Caution, incompressibility not ensured"
            
            sigDev=computeDeviatoricPart(np.array([SIG11[j+1,s,k],TAU[j+1,s,k],SIG22[j+1,s,k],SIG33[j+1,s,k]]))
            sigma = np.matrix([[sigDev[0],sigDev[1]/np.sqrt(2.),0.],[sigDev[1]/np.sqrt(2.),sigDev[2],0.],[0.,0.,sigDev[3]]])
            sigma = np.matrix([[SIG11[j+1,s,k],TAU[j+1,s,k],0.],[TAU[j+1,s,k],SIG22[j+1,s,k],0.],[0.,0.,SIG33[j+1,s,k]]])
            eigsigS[j+1,s,k,:]=computeEigenStresses(sigma)

            Nfinal=Niter
            if rcf2[j+1,s,k]>rcf2[j,s,k]:
                print "Simple wave condition violated"
                Nfinal=j
                break
            
        # plt.plot(time[1:],Eps[0,1:,s,k],label='eps11')
        # plt.plot(time[1:],Eps[1,1:,s,k],label='eps12')
        # plt.plot(time[1:],Eps[2,1:,s,k],label='eps22')
        # plt.plot(time[1:],Eps[3,1:,s,k],label='eps33')
        # plt.legend()
        # plt.grid()
        # plt.show()
        # plt.plot(time[1:],SIG11[1:,s,k],label='sig11')
        # plt.plot(time[1:],SIG22[1:,s,k],label='sig22')
        # plt.plot(time[1:],SIG33[1:,s,k],label='sig33')
        # plt.legend()
        # plt.grid()
        # plt.show()
        print "Final equivalent plastic strain after fast wave : ",plast
        fileName=path+'DPfastStressPlane_frame'+str(s)+'_Stress'+str(k)+'.pgf'
        pgfFilesList.append(fileName)
        if Nfinal/100>1:
            pas=Nfinal/100
            ranging=np.linspace(0,Nfinal-1,Nfinal/100,True,False,np.int)
        else:
            pas=1
            ranging=np.linspace(0,Nfinal-1,Nfinal,True,False,np.int)
        if exportPgf: export2pgfPlotFile(fileName,np.array([TAU[ranging,s,k],SIG11[ranging,s,k],SIG22[ranging,s,k],rcf2[ranging,s,k],LodeAngle_F[ranging,s,k]]),'sigma_12','sigma_11','sigma_22','p','Theta')
        fileName=path+'DPfastDevPlane_frame'+str(s)+'_Stress'+str(k)+'.pgf'
        deviatorPlots.append(fileName)
        dico={"xlabel":r'$\sigma_1$',"ylabel":r'$\sigma_2$',"zlabel":r'$\sigma_3$'}
        if exportPgf: export2pgfPlot3D(fileName,eigsigS[ranging,s,k,0],eigsigS[ranging,s,k,1],eigsigS[ranging,s,k,2],dico)

        radius_F[s,k]=sigy+H*plast
    
    cylindre=vonMisesYieldSurface(sigy)
    fileName=path+'CylindreDevPlane.pgf'
    deviatorPlots.append(fileName)
    if exportPgf: export2pgfPlot3D(str(fileName),cylindre[0,:],cylindre[1,:],cylindre[2,:])
    
    TAU_MAX_F=np.max(ordonnees)
    SIG_MAX_F=np.max(abscisses)
    
    plot_path=True
   
            
    ### SUBPLOTS SETTINGS
    if plot_path :
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
    
        ax1=plt.subplot2grid((2,2),(0,0))
        ax2=plt.subplot2grid((2,2),(0,1))
        ax4=plt.subplot2grid((2,2),(1,0),colspan=2)
        ax1.set_xlabel(r'$\sigma_{11}$')
        ax1.set_ylabel(r'$\sigma_{12}$')
        ax2.set_xlabel(r'$\sigma_{22}$')
        ax2.set_ylabel(r'$\sigma_{12}$')
        ax4.set_xlabel(r'$\Theta$')
        ax4.set_ylabel('p')
        ax1.grid()
        ax2.grid()
        ax4.grid()
        ax1.plot(sig[:,k],tau[:,k],'k')

        Delta=(4.*(nu**2-nu+1.)*sigy**2- 3.*(4.*nu**2-4.*nu+1.)*s22**2)
        sigMax=(s22*(1.+2.*nu-2.*nu**2)+np.sqrt(Delta))/(2.*(nu**2-nu+1.))
        sigMin=(s22*(1.+2.*nu-2.*nu**2)-np.sqrt(Delta))/(2.*(nu**2-nu+1.))
    
        # Sample stress component sig11
        sig[:,k]=np.linspace(sigMin,sigMax,Samples)
        # Compute shear stress satisfying the criterion given sig11 and sig22
        for i in range(Samples):
            s11=sig[i,k]
            delta=((nu-nu**2)*(s11+s22)**2 -s11**2-s22**2+s11*s22 + sigy**2)/3.
            if np.abs(delta)<10. : delta=np.abs(delta)
            tau[i,k]=np.sqrt(delta)
            if np.isnan(tau[i,k]):
                print "Nan ",delta,s11,s22

        fileName=path+'DPfast_yield0_s11s12_Stress'+str(k)+'.pgf'
        yields11_s12.append(fileName)
        if exportPgf: export2pgfPlotFile(fileName,np.array([tau[:,k],sig[:,k]]),'sigma_12','sigma_11')
        for p,i in enumerate(frames):
            sig0=sig[-1-i,k]
            Delta=(4.*(nu**2-nu+1.)*sigy**2- 3.*(4.*nu**2-4.*nu+1.)*sig0**2)
            s22max=(sig0*(1.+2.*nu-2.*nu**2)+np.sqrt(Delta))/(2.*(nu**2-nu+1.))
            s22min=(sig0*(1.+2.*nu-2.*nu**2)-np.sqrt(Delta))/(2.*(nu**2-nu+1.))
            s22=np.linspace(s22min,s22max,Samples)
            delta=((nu-nu**2)*(sig0+s22)**2 -s22**2-sig0**2+s22*sig0 + sigy**2)/3.
            if (np.abs(delta)<10.).any() : delta=np.abs(delta)
            s12=np.sqrt(delta)
        
            ## export to pgf file
            fileName=path+'DPfast_yield0_s22s12_frame'+str(p)+'_Stress'+str(k)+'.pgf'
            #yields22_s12.append(fileName)
            if exportPgf: export2pgfPlotFile(fileName,np.array([s12,s22]),'sigma_12','sigma_22')
            ax2.plot(s22,s12,'k')
            if updated_criterion :
                # plot updated yield surface that takes into account hardening in both planes
                sig0=SIG11[-1,p,k]
                maxCrit=0.5*(sig0*(2.*nu**2-2.*nu-1.)+E*Epsp33[-1,p,k]*(1.-2.*nu))/(nu-nu**2-1.)
                plast=plast_F[-1,p,k]
                Delta=(4.*(nu**2-nu+1.)*(sigy+H*plast)**2- 3.*(E*Epsp33[-1,p,k]+(1.-2.*nu)*sig0)**2)
                s22max=(sig0*(1.+2.*nu-2.*nu**2) +E*Epsp33[-1,p,k]*(2.*nu-1.) +np.sqrt(Delta))/(2.*(nu**2-nu+1.))
                s22min=(sig0*(1.+2.*nu-2.*nu**2)  +E*Epsp33[-1,p,k]*(2.*nu-1.) -np.sqrt(Delta))/(2.*(nu**2-nu+1.))
                s22=np.linspace(s22min,s22max,Samples)
                delta=(-(E*Epsp33[-1,p,k])**2 +E*Epsp33[-1,p,k]*(sig0+s22)*(2.*nu-1) + sig0*s22*(2.*nu+1.-2.*nu**2) + (sig0**2+s22**2)*(nu-nu**2-1.)+ (sigy+H*plast)**2)/3.
                if (np.abs(delta)<100.).any() : delta=np.abs(delta)
                s12=np.sqrt(delta)
                ax2.plot(s22,s12,color=col[p],linestyle='--')
                ax2.plot([maxCrit,maxCrit],[0.,tauM],color=col[p],linestyle='-.')
                ## export to pgf file
                fileName=path+'DPfast_yieldfin_s22s12_frame'+str(p)+'_Stress'+str(k)+'.pgf'
                if exportPgf: export2pgfPlotFile(fileName,np.array([s12,s22]),'sigma_12','sigma_22')
                ## The same on ax1
                s22=SIG22[-1,p,k]
                Delta=(4.*(nu**2-nu+1.)*(sigy+H*plast)**2- 3.*(E*Epsp33[-1,p,k]+(1.-2.*nu)*s22)**2)
                sigMax=(s22*(1.+2.*nu-2.*nu**2) +E*Epsp33[-1,p,k]*(2.*nu-1.) +np.sqrt(Delta))/(2.*(nu**2-nu+1.))
                sigMin=(s22*(1.+2.*nu-2.*nu**2) +E*Epsp33[-1,p,k]*(2.*nu-1.) -np.sqrt(Delta))/(2.*(nu**2-nu+1.))
                s11=np.linspace(sigMin,sigMax,Samples)
                maxCrit=0.5*(s22*(2.*nu**2-2.*nu-1.)+E*Epsp33[-1,p,k]*(1.-2.*nu))/(nu-nu**2-1.)
                delta=(-(E*Epsp33[-1,p,k])**2 +E*Epsp33[-1,p,k]*(s11+s22)*(2.*nu-1) + s11*s22*(2.*nu+1.-2.*nu**2) + (s11**2+s22**2)*(nu-nu**2-1.)+ (sigy+H*plast)**2)/3.
                if (np.abs(delta)<10.).any() : delta=np.abs(delta)
                s12=np.sqrt(delta)
                ## export to pgf file
                fileName=path+'DPfast_yieldfin_s11s12_frame'+str(p)+'_Stress'+str(k)+'.pgf'
                if exportPgf: export2pgfPlotFile(fileName,np.array([s12,s11]),'sigma_12','sigma_11')
                ax1.plot(s11,s12,color=col[p],linestyle='--')
                ax1.plot([maxCrit,maxCrit],[0.,tauM],color=col[p],linestyle='-.')
            ax1.plot(SIG11[:,p,k],TAU[:,p,k],color=col[p],lw=2.5)
            ax2.plot(SIG22[:,p,k],TAU[:,p,k],color=col[p],lw=2.5)
            ax4.plot(LodeAngle_F[:,p,k],plast_F[:,p,k],color=col[p],lw=2.5)
            ax3.plot(eigsigS[:,p,k,0],eigsigS[:,p,k,1],eigsigS[:,p,k,2],color=col[p],lw=2.5)
            # plot s11=0 point in the space of stresses
            cylindre=vonMisesYieldSurface(radius_F[p,k])
            ax3.plot(cylindre[0,:],cylindre[1,:],cylindre[2,:], color=col[p],linestyle='--')
        ax3.plot([-sigy,sigy],[0.,0.],[0.,0.],color="k",linestyle="--",lw=1.)
        ax3.plot([0.,0.],[-sigy,sigy],[0.,0.],color="k",linestyle="--",lw=1.)
        ax3.plot([-radius,radius],[radius,-radius],[0.,0.],color="k",linestyle="--",lw=1.)
        plt.suptitle(r'Loading paths through fast waves for $\sigma_{22}$ ='+'{:.2e}'.format(sig22[k])+'Pa.', fontsize=24.)
        plt.show()

    else:
        fig = plt.figure()
        ax=plt.subplot2grid((1,1),(0,0),projection='3d')
        ax.set_xlabel(r'$\sigma_{11}$',size=24.)
        ax.set_ylabel(r'$\sigma_{12}$',size=24.)
        ax.set_zlabel(r'$\sigma_{22}$',size=24.)
        for l in range(len(sig22)):
            ax.plot(sig[:,l],tau[:,l],sig22[l],color="k")
        for p,i in enumerate(frames):
            ax.plot(SIG11[:,p,k],TAU[:,p,k],SIG22[:,p,k],color=col[p],lw=2.5)
        #plt.show()

    ## sig22 value will change here
    subtitle=[r'(a) ($\sigma_{11},\sigma_{12}$) plane',r'(b) ($\sigma_{22},\sigma_{12}$) plane',r'(c) Deviatoric plane']

    srcX=['sigma_11','sigma_22']
    srcY=['sigma_12','sigma_12']

    name1='DPfastWaves_sig11_tau'+str(k)+'.tex'
    name2='DPfastWaves_sig22_tau'+str(k)+'.tex'
    name3='DPfastWaves_deviator'+str(k)+'.tex'
    names=[[name1,name2],name3]
    
    files1=np.concatenate([pgfFilesList,yields11_s12])
    files2=np.concatenate([pgfFilesList,yields22_s12]) 
    pgfFiles=[[files1,files2],deviatorPlots]
    xlabels=[['$\sigma_{11} (Pa)$','$\sigma_{22}  (Pa)$'],'$s_1 $'] #size=number of .tex files
    ylabels=[['$\sigma_{12}  (Pa)$','$\sigma_{12}  (Pa)$'],'$s_2 $'] #size=number of .tex files
    zlabels=[['',''],'$s_3$'] #size=number of .tex files
    
    TauMax=1.*sigy
    buildTeXFiles2(names,pgfFiles,xlabels,ylabels,zlabels,srcX,srcY,TauMax)
    
    pgfFilesList=[];yields11_s12=[];deviatorPlots=[];yields22_s12=[]
