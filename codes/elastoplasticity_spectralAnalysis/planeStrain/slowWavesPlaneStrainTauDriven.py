# !/usr/bin/python
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

def export2DTeXFile(fileName,xFields,xlabel,ylabel1,ylabel2,yfields1,yfields2,*kwargs):
    TeXFile=open(fileName,"w")
    n_fields = np.shape(yfields1)[0]
    n_labels = np.shape(kwargs)[0]
    # Define Paul Tol's colors (purple to red)
    marker=['none','none','none','none','none','x','pentagone*','none','triangle*']
    style=['solid','solid','solid','solid','solid','only marks','solid','solid']
    thickness=['very thick','very thick','thick','very thick','very thick','thick','thin','thin','thick']
    couleur=['Red','Blue','Duck','Purple','Green','Duck','Yellow']
    TeXFile.write(r'\begin{tikzpicture}[scale=.8]');TeXFile.write('\n')
    TeXFile.write(r'\begin{axis}[xlabel='+str(xlabel)+',ylabel='+str(ylabel1)+',ymajorgrids=true,xmajorgrids=true,legend pos= outer north east]');TeXFile.write('\n')
    legend=''
    for i in range(n_fields):
        if i==0:
            legend=legend+kwargs[0][i]
        else:
            legend=legend+','+kwargs[0][i]
        TeXFile.write(r'\addplot['+str(couleur[i])+','+str(thickness[i])+',mark='+str(marker[i])+','+str(style[i])+'] coordinates {')
        for j in range(np.shape(yfields1[i])[0]):
            TeXFile.write('('+str(xFields[i][j])+','+str(yfields1[i][j])+') ')
        TeXFile.write('};\n')
    TeXFile.write(r'\legend{'+str(legend)+'}')
    TeXFile.write('\n')    
    TeXFile.write(r'\end{axis}')
    TeXFile.write('\n')
    ## Seond axis
    n_fields = np.shape(yfields2)[0]
    n_labels = np.shape(kwargs)[0]
    # Define Paul Tol's colors (purple to red)
    #marker=['','none','none','none','none','x','pentagone*','none','triangle*']
    style=['dashed','dashed','dashed','dashed','dashed','only marks','dashed','dashed']
    thickness=['very thick','very thick','thick','very thick','very thick','thick','thin','thin','thick']
    TeXFile.write(r'\begin{axis}[xlabel='+str(xlabel)+',ylabel='+str(ylabel2)+',ymajorgrids=true,xmajorgrids=true,axis y line*=right]');TeXFile.write('\n')
    legend=''
    for i in range(n_fields):
        TeXFile.write(r'\addplot['+str(couleur[i])+','+str(thickness[i])+',mark='+str(marker[i])+','+str(style[i])+'] coordinates {')
        for j in range(np.shape(yfields2[i])[0]):
            TeXFile.write('('+str(xFields[i][j])+','+str(yfields2[i][j])+') ')
        TeXFile.write('};\n')
    TeXFile.write(r'\legend{'+str(legend)+'}')
    TeXFile.write('\n')    
    TeXFile.write(r'\end{axis}')
    TeXFile.write('\n')

    TeXFile.write('\end{tikzpicture}')
    TeXFile.write('\n')
    TeXFile.write('%%% Local Variables:')
    TeXFile.write('\n')
    TeXFile.write('%%% mode: latex')
    TeXFile.write('\n')
    TeXFile.write('%%% TeX-master: "../../mainManuscript"')
    TeXFile.write('\n')
    TeXFile.write('%%% End:')
    TeXFile.write('\n')
    TeXFile.close()


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

def vonMisesYieldSurface(sigma):
    ## Version one : build a cylinder in sigma eigenspace and turn it around axis
    
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

def loadingPath(dtau,sigma,sigman,h,lamb,mu,beta,tangent):
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
    # return dsig11,dsig22,depsp33
    return np.array([psi11*dtau,psi22*dtau,depsp33])

def integrateODE(dtau,sig0,tau0,sig22_0,sig33,epsp33,nu,E,H,lamb,mu,beta,tangent):
    sigma=np.array([sig0,sig22_0,epsp33])
    # computePsiSlow(sig12,sigma,sig33,lamb,mu,beta,tangent)
    # subdivision of time step
    sub_steps = 1
    dTAU = dtau/sub_steps
    theta = 1.0
    for i in range(sub_steps):
        ## Nonlinear solution procedure
        ## R = s^{n+1} - s^{n} - RHS
        
        R=lambda x: x - sigma - theta*loadingPath(dTAU,np.array([sigma[0],tau0,sigma[1],nu*(sigma[0]+sigma[1])-E*sigma[2]]),np.array([x[0],tau0+dTAU,x[1],nu*(x[0]+x[1])-E*x[2]]),H,lamb,mu,beta,tangent)
        #pdb.set_trace()
        
        solution = scipy.optimize.fsolve(R,np.array([sigma[0],sigma[1],epsp33]))
        sigma = solution
    return solution[0],solution[1],solution[2],nu*(solution[0]+solution[1])-E*solution[2]

def computeSpeed(sigma,lamb,mu,beta,tangent):
    # sig12 driven
    n1=1.;n2=0.
    sig11=sigma[0,0];sig22=sigma[1,1];sig12=sigma[0,1];sig33=sigma[2,2]
    H=tangentModulus(np.array([sig11,sig12,sig22,sig33]),lamb,mu,beta,tangent)
    C=acousticTensor(H,np.array([n1,n2]))
    eigenf,eigens=acousticEigenStructure(C)
    return eigens[0]


def computeSpeedSlow(sig12,sigma,sig33,lamb,mu,beta,tangent):
    # sig12 driven
    n1=1.;n2=0.
    sig11=sigma[0];sig22=sigma[1]
    H=tangentModulus(np.array([sig11,sig12,sig22,sig33]),lamb,mu,beta,tangent)
    C=acousticTensor(H,np.array([n1,n2]))
    eigenf,eigens=acousticEigenStructure(C)
    rcs2=eigens[0]
    return rcs2

def computeLodeAngle(sig11,sig22,sig12,sig33):
    # deviatoric stress
    sDev=computeDeviatoricPart(np.array([sig11,sig12,sig22,sig33]))
    s11=sDev[0];s12=sDev[1]/np.sqrt(2.);s22=sDev[2];s33=sDev[3]
    sig=computeEigenStresses(np.matrix([[s11,s12,0.],[s12,s22,0.],[0.,0.,s33]]))
    #sig=np.linalg.eig(np.matrix([[s11,s12,0.],[s12,s22,0.],[0.,0.,s33]]))[0]
    J3=s33*(s11*s22-s12**2) ;J2=0.5*np.dot(sDev,sDev) ; sqrtJ2=np.sqrt(0.5*np.dot(sDev,sDev))
    sin3Theta=0.5*J3*(3./J2)**(3./2.)
    theta=np.arcsin(sin3Theta)/3.
    tan=(1./np.sqrt(3.))*(2.*(sig[1]-sig[2])/(sig[0]-sig[2])-1.)
    theta=-np.sign(tan)*np.arccos((1./2.)*np.sqrt(3.)*J3/(sqrtJ2**3))/3.
    theta=theta*360./(2.*np.pi)
    cos3Theta=(27./2.)*J3/(np.sqrt(3.*J2)**3)
    theta=(1.-(2/np.pi)*np.arccos(cos3Theta))*360./(2.*np.pi)
    return cos3Theta

def computeLodeAngle2(sig11,sig22,sig12,sig33):
    # deviatoric stress
    sDev=computeDeviatoricPart(np.array([sig11,sig12,sig22,sig33]))
    s11=sDev[0];s12=sDev[1]/np.sqrt(2.);s22=sDev[2];s33=sDev[3]
    sig=computeEigenStresses(np.matrix([[s11,s12,0.],[s12,s22,0.],[0.,0.,s33]]))
    #sig=np.linalg.eig(np.matrix([[s11,s12,0.],[s12,s22,0.],[0.,0.,s33]]))[0]
    J3=s33*(s11*s22-s12**2) ;J2=0.5*np.dot(sDev,sDev) ; sqrtJ2=np.sqrt(0.5*np.dot(sDev,sDev))
    sin3Theta=0.5*J3*(3./J2)**(3./2.)
    theta=np.arcsin(sin3Theta)/3.
    tan=(1./np.sqrt(3.))*(2.*(sig[1]-sig[2])/(sig[0]-sig[2])-1.)
    theta=-np.sign(tan)*np.arccos((1./2.)*np.sqrt(3.)*J3/(sqrtJ2**3))/3.
    theta=theta*360./(2.*np.pi)
    cos3Theta=(27./2.)*J3/(np.sqrt(3.*J2)**3)
    theta=(1.-(2/np.pi)*np.arccos(cos3Theta))*360./(2.*np.pi)
    sin3Theta=np.sign(tan)*np.sign(cos3Theta)*np.sqrt(1.-cos3Theta**2)
    return tan

def updateEquivalentPlasticStrain(sig,sign,H):
    # sig=[sig11^n , sqrt(2)*sig12^n , sig22 , sig33^n]
    # sign=[sig11^n+1 , sqrt(2)*sig12^n+1 , sig22 , sig33^n+1]
    sigDevn=computeDeviatoricPart(np.array([sign[0],sign[1]/np.sqrt(2.),sign[2],sign[3]]))
    sigDev=computeDeviatoricPart(np.array([sig[0],sig[1]/np.sqrt(2.),sig[2],sig[3]]))
    norm=np.sqrt(np.dot(sigDevn,sigDevn))
    flow=sigDev/norm
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
    # eigenvectors
    v1=np.array([sig[0,0]-sig[1,1]-np.sqrt(delta),sig[0,1],0])
    v2=np.array([sig[0,0]-sig[1,1]+np.sqrt(delta),sig[0,1],0])
    return np.array([s1,s2,s3])

def eigenVects(sig):
    #    | sig11 sig12   0   |
    #sig=| sig12 sig22   0   |
    #    |   0     0   sig33 |
    s3=sig[2,2]
    delta=(sig[0,0]-sig[1,1])**2+4.*sig[0,1]**2
    s1=0.5*(sig[0,0]+sig[1,1]-np.sqrt(delta))
    s2=0.5*(sig[0,0]+sig[1,1]+np.sqrt(delta))
    # eigenvectors
    v1=np.array([sig[0,0]-sig[1,1]-np.sqrt(delta),sig[0,1],0])
    v2=np.array([sig[0,0]-sig[1,1]+np.sqrt(delta),sig[0,1],0])
    v3=np.array([0.,0.,1.])
    return v1,v2,v3


def computePlasticResidual(epsp33,sig,epsp33n,sign,E,H,nu):
    # sig = [sig11 , sig12*sqrt(2) , sig22 , sig33] (previous time step)
    # sig33n = nu*(sig11+sig22)-E*epsp33n (updated time step)
    sign[3]=nu*(sign[0]+sign[2])-E*epsp33n
    sigDev=computeDeviatoricPart(np.array([sig[0],sig[1]/np.sqrt(2.),sig[2],sig[3]]))
    sigDevn=computeDeviatoricPart(np.array([sign[0],sign[1]/np.sqrt(2.),sign[2],sign[3]]))
    norm=np.sqrt(np.dot(sigDevn,sigDevn))
    flow=sigDevn/norm
    dSig=sign-sig
    dSig[1]*=np.sqrt(2.)
    #dSig=sigDevn-sigDev
    dp=(1./H)*np.sqrt(3./2.)*np.dot(flow,dSig)
    res=epsp33n-epsp33-np.sqrt(3./2.)*flow[3]*dp
    return res

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
            

from mpl_toolkits.mplot3d import proj3d
def orthogonal_proj(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    return np.array([[1,0,0,0],
                        [0,1,0,0],
                        [0,0,a,b],
                        [0,0,0,zback]])
proj3d.persp_transformation = orthogonal_proj

Samples=5#

# Sample constant stress component sig22
sig22=np.linspace(0.,sigy*np.sqrt(4*(nu**2-nu+1.))/np.sqrt(3.*(4.*nu**2-4.*nu+1.)),Samples)
sig22=np.linspace(-sigy*np.sqrt(4*(nu**2-nu+1.))/np.sqrt(3.*(4.*nu**2-4.*nu+1.)),sigy*np.sqrt(4*(nu**2-nu+1.))/np.sqrt(3.*(4.*nu**2-4.*nu+1.)),Samples)

#sig22=[0.]
Samples*=10
sig=np.zeros((Samples,Samples/10))
tau=np.zeros((Samples,Samples/10))

frames=[5,20,30,40]
#frames=[10,15,20,25,30,35]
#frames=[5]
col=["r","g","b","y","c","m","k","p"]
# purple to red
col=['#781C81','#3F60AE','#539EB6','#6DB388','#CAB843','#E78532','#D92120']
tauM=1.5*sigy/np.sqrt(3.)
sigM=1.5*sigy/np.sqrt(1-nu+nu**2)
tauM=sigM
Niter=5000
TAU=np.zeros((Niter,len(frames),len(sig22)))
SIG11=np.zeros((Niter,len(frames),len(sig22)))
SIG22=np.zeros((Niter,len(frames),len(sig22)))
SIG33=np.zeros((Niter,len(frames),len(sig22)))
MaxCrit=np.zeros((Niter,len(frames),len(sig22)))
eigsigS=np.zeros((Niter,len(frames),len(sig22),3))
eigsigDevS=np.zeros((Niter,len(frames),len(sig22),3))
criterionS=np.zeros((Niter,len(frames),len(sig22)))
PsiS=np.zeros((Samples,len(sig22)))

plast_S=np.zeros((Niter,len(frames),len(sig22)))
rcs2=np.zeros((Niter,len(frames),len(sig22)))
Epsp33=np.zeros((Niter,len(frames),len(sig22)))
Eps=np.zeros((4,Niter,len(frames),len(sig22)))
LodeAngle_S=np.zeros((Niter,len(frames),len(sig22)))
LodeAngle_S2=np.zeros((Niter,len(frames),len(sig22)))
speed_S=np.zeros((Niter,len(frames),len(sig22)))

# Boolean to plot the upadted yield surface
updated_criterion=True
for k in range(len(sig22)):
    s22=sig22[k]
    
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
        if np.isnan(tau[i,k]): print "Nan ",delta,s11,s22,i,k


tangent='planeStrain'
## LOADING PATHS PLOTS
pgfFilesList=[]
yields11_s12=[]
yields22_s12=[]
deviatorPlots=[]
for k in range(len(sig22)-1)[1:]:
#for k in range(len(sig22)):
    # if k==1 or k==3:
    #     continue
    s22=sig22[k]
    print "sigma22=",s22
    sigM=1.25*np.max(sig[:,k])
    tauM=1.5*np.max(tau[:,k])
    ## For each value of sig22 trace the loading paths given by psis from yield surface to an arbitrary shear stress level
    approx=np.zeros((len(frames),2))
    ordonnees=np.zeros((len(frames),Samples))
    abscisses=np.zeros((len(frames),Samples))
    radius_S=np.zeros(len(frames))
    for s,i in enumerate(frames):
        # if s22==0:
        #     continue
        sig0=sig[-1-i,k]
        tau0=tau[-1-i,k]
        
        dtau=(tauM-tau0)/Niter
        
        TAU[:,s,k]=np.linspace(tau0,tauM,Niter)
        
        SIG11[0,s,k]=sig0
        SIG22[0,s,k]=s22
        
        sig33=nu*(SIG11[0,s,k]+SIG22[0,s,k])
        SIG33[0,s,k]=sig33
        
        sig0=SIG11[0,s,k];
        MaxCrit[0,s,k]=0.5*(sig0*(2.*nu**2-2.*nu-1.))/(nu-nu**2-1.)

        sigma = np.matrix([[SIG11[0,s,k],TAU[0,s,k],0.],[TAU[0,s,k],SIG22[0,s,k],0.],[0.,0.,sig33]])
        rcs2[0,s,k]=np.sqrt(computeSpeed(sigma,lamb,mu,beta,tangent)/rho)
                    
        eigsigS[0,s,k,:]=computeEigenStresses(sigma)
        sigDev=computeDeviatoricPart(np.array([SIG11[0,s,k],TAU[0,s,k],SIG22[0,s,k],SIG33[0,s,k]]))
        sigma = np.matrix([[sigDev[0],sigDev[1]/np.sqrt(2.),0.],[sigDev[1]/np.sqrt(2.),sigDev[2],0.],[0.,0.,sigDev[3]]])
        eigsig=computeEigenStresses(sigma)
        eigsigDevS[0,s,k,:]=eigsig
        LodeAngle_S[0,s,k]=computeLodeAngle(sigma[0,0],SIG22[0,s,k],sigma[0,1],sig33)
        LodeAngle_S2[0,s,k]=computeLodeAngle2(sigma[0,0],SIG22[0,s,k],sigma[0,1],sig33)
        
        speed_S[0,s,k]=computeSpeedSlow(TAU[0,s,k],[SIG11[0,s,k],SIG22[0,s,k]],sig33,lamb,mu,beta,tangent)
            
        plast=0.
        epsp33=0.
        for j in range(Niter-1):
            SIG11[j+1,s,k],SIG22[j+1,s,k],epsp33,SIG33[j+1,s,k]=integrateODE(dtau,SIG11[j,s,k],TAU[j,s,k],SIG22[j,s,k],SIG33[j,s,k],epsp33,nu,E,H,lamb,mu,beta,tangent)
            sig0=SIG11[j+1,s,k];
            MaxCrit[j+1,s,k]=0.5*(sig0*(2.*nu**2-2.*nu-1.)+E*epsp33*(1.-2.*nu))/(nu-nu**2-1.)
                
            sigma = np.array([SIG11[j,s,k],np.sqrt(2.)*TAU[j,s,k],SIG22[j,s,k],SIG33[j,s,k]])
            
            sig33=SIG33[j+1,s,k]
            
            sigman=np.array([SIG11[j+1,s,k],TAU[j+1,s,k]*np.sqrt(2.),SIG22[j+1,s,k],sig33])
            speed_S[j+1,s,k]=computeSpeedSlow(TAU[j+1,s,k],[SIG11[j+1,s,k],SIG22[j+1,s,k]],sig33,lamb,mu,beta,tangent)

            if speed_S[j+1,s,k]>speed_S[j,s,k]:
                print "Simple wave condition violated"
                break
            dp=updateEquivalentPlasticStrain(sigma,sigman,H)
            plast+=dp

            criterionS[j+1,s,k]=computeCriterion(SIG11[j+1,s,k],SIG22[j+1,s,k],TAU[j+1,s,k],sig33,sigy+H*plast)
            plast_S[j+1,s,k]=plast
            
            Epsp33[j+1,s,k]=epsp33
            LodeAngle_S[j+1,s,k]=computeLodeAngle(sigman[0],sigman[2],sigman[1]/np.sqrt(2.),sig33)
            LodeAngle_S2[j+1,s,k]=computeLodeAngle2(sigman[0],sigman[2],sigman[1]/np.sqrt(2.),sig33)
            # Eigenvalues of sigma (for deviatoric plane plots)
            sigma = np.matrix([[SIG11[j+1,s,k],TAU[j+1,s,k],0.],[TAU[j+1,s,k],SIG22[j+1,s,k],0.],[0.,0.,SIG33[j+1,s,k]]])
            rcs2[j+1,s,k]=np.sqrt(computeSpeed(sigma,lamb,mu,beta,tangent)/rho)

            #Eps=np.zeros((4,Niter,len(frames),len(sig22)))
            Eps[:,j+1,s,k]=updateEpsilon(Eps[:,j,s,k],np.array([SIG11[j,s,k],TAU[j,s,k],SIG22[j,s,k],SIG33[j,s,k]]),np.array([SIG11[j+1,s,k],TAU[j+1,s,k],SIG22[j+1,s,k],SIG33[j+1,s,k]]),E,nu,H)
            
            #pdb.set_trace()
            eigsigS[j+1,s,k,:]=computeEigenStresses(sigma)
            # sigDev=computeDeviatoricPart(np.array([SIG11[j+1,s,k],TAU[j+1,s,k],SIG22[j+1,s,k],SIG33[j+1,s,k]]))
            # sigma = np.matrix([[sigDev[0],sigDev[1]/np.sqrt(2.),0.],[sigDev[1]/np.sqrt(2.),sigDev[2],0.],[0.,0.,sigDev[3]]])
            eigsigDevS[j+1,s,k,:]=computeEigenStresses(sigma)
            radi=np.sqrt(np.dot(eigsigDevS[j,s,k,:],eigsigDevS[j,s,k,:]))
            radi2=np.sqrt(2./3.)*sigy
            #pdb.set_trace()
        
        time=np.linspace(0,j+1,Niter)
        # plt.plot(time[1:],Eps[0,1:,s,k],label='eps11')
        # plt.plot(time[1:],Eps[1,1:,s,k],label='eps12')
        # plt.plot(time[1:],Eps[2,1:,s,k],label='eps22')
        # plt.plot(time[1:],Eps[3,1:,s,k],label='eps33')
        # plt.legend()
        # plt.grid()
        # plt.show()
        print "Final equivalent plastic strain after slow wave : ",plast
        fileName=path+'DPslowStressPlane_frame'+str(s)+'_Stress'+str(k)+'.pgf'
        ## color bar of p
        #export2pgfPlotFile(fileName,np.array([TAU[0:-1:Niter/100,s,k],SIG11[0:-1:Niter/100,s,k],SIG22[0:-1:Niter/100,s,k],plast_S[0:-1:Niter/100,s,k],LodeAngle_S[0:-1:Niter/100,s,k]]),'sigma_12','sigma_11','sigma_22','p','Theta')
        ## color bar of rcs2
        print np.min(LodeAngle_S[0:-1:Niter/100,s,k]),np.max(LodeAngle_S[0:-1:Niter/100,s,k])
        export2pgfPlotFile(fileName,np.array([TAU[0:-1:Niter/100,s,k],SIG11[0:-1:Niter/100,s,k],SIG22[0:-1:Niter/100,s,k],rcs2[0:-1:Niter/100,s,k],LodeAngle_S[0:-1:Niter/100,s,k],eigsigS[0:-1:Niter/100,s,k,2]+eigsigS[0:-1:Niter/100,s,k,1]+eigsigS[0:-1:Niter/100,s,k,0],eigsigDevS[0:-1:Niter/100,s,k,2],MaxCrit[0:-1:Niter/100,s,k]]),'sigma_12','sigma_11','sigma_22','p','Theta','hydro','s3','maxCrit')
        pgfFilesList.append(fileName)
        fileName=path+'DPslowDevPlane_frame'+str(s)+'_Stress'+str(k)+'.pgf'
        dico={"xlabel":r'$s_1$',"ylabel":r'$s_2$',"zlabel":r'$s_3$'}
        export2pgfPlot3D(fileName,eigsigS[0:-1:Niter/100,s,k,0],eigsigS[0:-1:Niter/100,s,k,1],eigsigS[0:-1:Niter/100,s,k,2],dico)
        deviatorPlots.append(fileName)
        """
        if k==1:
            fileName='slowWave_DevPlane_frame'+str(s)+'_Stress'+str(k)+'.pgf'
            dico={"xlabel":r'$\sigma_1$',"ylabel":r'$\sigma_2$',"zlabel":r'$\sigma_3$'}
            export2pgfPlot3D(fileName,eigsigDevS[0:-1:Niter/100,s,k,0],eigsigDevS[0:-1:Niter/100,s,k,1],eigsigDevS[0:-1:Niter/100,s,k,2],dico)
        """
        radius_S[s]=sigy+H*plast
    # plt.plot(TAU[:,0,k],maxCrit[:,0,k],'r')
    # plt.plot(TAU[:,0,k],SIG22[:,0,k],'r--')
    # plt.plot(TAU[:,1,k],maxCrit[:,1,k],'g')
    # plt.plot(TAU[:,1,k],SIG22[:,1,k],'g--')
    # plt.plot(TAU[:,2,k],maxCrit[:,2,k],'b')
    # plt.plot(TAU[:,2,k],SIG22[:,2,k],'b--')
    # plt.plot(TAU[:,3,k],maxCrit[:,3,k],'k')
    # plt.plot(TAU[:,3,k],SIG22[:,3,k],'k--')
    # plt.plot(TAU[:,0,k],SIG11[:,0,k]+SIG22[:,0,k]+SIG33[:,0,k],'r')
    # plt.grid()
    # plt.show()
        
    if k==1 or k==2 or k==3:
        # ran=Niter
        # plt.plot(eigsigDevS[0:ran,0,k,2],SIG22[0:ran,0,k],'r')
        # #plt.plot(SIG22[0:ran,0,k],TAU[0:ran,0,k],'r--')
        # plt.plot(eigsigDevS[0:ran,3,k,2],SIG22[0:ran,3,k],'b')
        # #plt.plot(SIG22[0:ran,3,k],TAU[0:ran,3,k],'b--')
        # plt.grid()
        # plt.show()
        #pdb.set_trace()

        legend=['loading path 1','loading path 2','loading path 3','loading path 4']
        #export2DTeXFile('maxCrit.tex.tex',np.array([LodeAngle_S[0:-1:Niter/100,0,k]]),r'$\Theta ()$',r'$\sigma_{11} \: (Pa)$','',np.array([SIG11[0:-1:Niter/100,0,k]]),legend)

        print "export Lode Angle"
        ## export Graph of Lode Angle 0:-1:Niter/100
        legend=['loading path 1','loading path 2','loading path 3','loading path 4']
        #export2DTeXFile('LodeAngle_s11.tex',np.array([LodeAngle_S[0:-1:Niter/100,0,k],LodeAngle_S[0:-1:Niter/100,1,k],LodeAngle_S[0:-1:Niter/100,2,k],LodeAngle_S[0:-1:Niter/100,3,k]]),r'$\Theta ()$',r'$\sigma_{11} \: (Pa)$','',np.array([SIG11[0:-1:Niter/100,0,k]/max(SIG11[0:-1:Niter/100,0,k]),SIG11[0:-1:Niter/100,1,k]/max(SIG11[0:-1:Niter/100,1,k]),SIG11[0:-1:Niter/100,2,k]/max(SIG11[0:-1:Niter/100,2,k]),SIG11[0:-1:Niter/100,3,k]/max(SIG11[0:-1:Niter/100,3,k])]),legend)
        #export2DTeXFile('LodeAngle_s22.tex',np.array([LodeAngle_S[0:-1:Niter/100,0,k],LodeAngle_S[0:-1:Niter/100,1,k],LodeAngle_S[0:-1:Niter/100,2,k],LodeAngle_S[0:-1:Niter/100,3,k]]),r'$\Theta ()$',r'$\sigma_{22} \: (Pa)$','',np.array([SIG22[0:-1:Niter/100,0,k]/np.max(SIG22[0:-1:Niter/100,0,k]),SIG22[0:-1:Niter/100,1,k]/np.max(SIG22[0:-1:Niter/100,1,k]),SIG22[0:-1:Niter/100,2,k]/np.max(SIG22[0:-1:Niter/100,2,k]),SIG22[0:-1:Niter/100,3,k]/np.max(SIG22[0:-1:Niter/100,3,k])]),legend)
        legend=['loading path 1','loading path 2','loading path 3','loading path 4']
        #export2DTeXFile('LodeAngle_s11.tex',np.array([LodeAngle_S[0:-1:Niter/100,0,k]]),r'$\Theta ()$',r'$\sigma_{11} \: (Pa)$','',np.array([SIG11[0:-1:Niter/100,0,k]]),legend)

        ## fileName,xFields,xlabel,ylabel1,ylabel2,yfields1,yfields2,*kwargs):
        #export2DTeXFile('LodeAngle.tex',np.array([TAU[0:-1:Niter/200,0,k],TAU[0:-1:Niter/200,1,k],TAU[0:-1:Niter/200,2,k],TAU[0:-1:Niter/200,3,k]]),r'$\sigma_{12} \: (Pa)$',r'$\sigma_{22} \: (Pa)$',r'$\Theta ()$',np.array([SIG22[0:-1:Niter/200,0,k],SIG22[0:-1:Niter/200,1,k],SIG22[0:-1:Niter/200,2,k],SIG22[0:-1:Niter/200,3,k]]),np.array([LodeAngle_S[0:-1:Niter/200,0,k],LodeAngle_S[0:-1:Niter/200,1,k] ,LodeAngle_S[0:-1:Niter/200,2,k],LodeAngle_S[0:-1:Niter/200,3,k]]),legend)
        # export2DTeXFile('LodeAngle_s11.tex',np.array([LodeAngle_S[0:-1:Niter/100,0,k],LodeAngle_S[0:-1:Niter/100,1,k],LodeAngle_S[0:-1:Niter/100,2,k],LodeAngle_S[0:-1:Niter/100,3,k]]),r'$\Theta ()$',r'$\sigma_{11} \: (Pa)$','',np.array([SIG11[0:-1:Niter/100,0,k],SIG11[0:-1:Niter/100,1,k],SIG11[0:-1:Niter/100,2,k],SIG11[0:-1:Niter/100,3,k]]),legend)
        # export2DTeXFile('LodeAngle_s22.tex',np.array([LodeAngle_S[0:-1:Niter/100,0,k],LodeAngle_S[0:-1:Niter/100,1,k],LodeAngle_S[0:-1:Niter/100,2,k],LodeAngle_S[0:-1:Niter/100,3,k]]),r'$\Theta ()$',r'$\sigma_{11} \: (Pa)$','',np.array([SIG22[0:-1:Niter/100,0,k],SIG22[0:-1:Niter/100,1,k],SIG22[0:-1:Niter/100,2,k],SIG22[0:-1:Niter/100,3,k]]),legend)

    cylindre=vonMisesYieldSurface(sigy)
    fileName=path+'CylindreDevPlane.pgf'
    export2pgfPlot3D(str(fileName),cylindre[0,:],cylindre[1,:],cylindre[2,:])
    deviatorPlots.append(fileName)
    TAU_MAX_S=np.max(ordonnees)
    SIG_MAX_S=np.max(abscisses)


    plot_path=True
    ### SUBPLOTS SETTINGS
    if plot_path :
        fig = plt.figure(figsize=(10,10))
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
        ax3.scatter(sigy,0.,0.);ax3.scatter(-sigy,0.,0.)
        ax3.scatter(0.,sigy,0.);ax3.scatter(0.,-sigy,0.)
        ax3.scatter(0.,0.,sigy);ax3.scatter(0.,0.,-sigy)
        ax3.plot([0.,sigy],[0.,sigy],[0.,sigy],color="k")
        ax3.set_xlabel(r'$s_1$',size=24.)
        ax3.set_ylabel(r'$s_2$',size=24.)
        ax3.set_zlabel(r'$s_3$',size=24.)
        
        fig = plt.figure()
        ax4=plt.subplot2grid((2,2),(1,0),colspan=2)
        ax1=plt.subplot2grid((2,2),(0,0))
        ax2=plt.subplot2grid((2,2),(0,1))
        ax1.set_xlabel(r'$\sigma_{11}$')
        ax1.set_ylabel(r'$\sigma_{12}$')
        ax2.set_xlabel(r'$\sigma_{22}$')
        ax2.set_ylabel(r'$\sigma_{12}$')
        ax4.set_xlabel(r'$\tau$')
        ax4.set_ylabel('p')
        ax1.grid()
        ax2.grid()
        ax4.grid()

        ax1.plot(sig[:,k],tau[:,k],'k')
        ## export to pgf file
        fileName=path+'DPslow_yield0_s11s12_Stress'+str(k)+'.pgf'
        export2pgfPlotFile(fileName,np.array([tau[:,k],sig[:,k]]),'sigma_12','sigma_11')
        yields11_s12.append(fileName)
        for p,i in enumerate(frames):
            plast=plast_S[-1,p,k]
            sig0=sig[-1-i,k]
            Delta=(4.*(nu**2-nu+1.)*sigy**2- 3.*(4.*nu**2-4.*nu+1.)*sig0**2)
            s22max=(sig0*(1.+2.*nu-2.*nu**2)+np.sqrt(Delta))/(2.*(nu**2-nu+1.))
            s22min=(sig0*(1.+2.*nu-2.*nu**2)-np.sqrt(Delta))/(2.*(nu**2-nu+1.))
            s22=np.linspace(s22min,s22max,Samples)
            delta=((nu-nu**2)*(sig0+s22)**2 -s22**2-sig0**2+s22*sig0 + sigy**2)/3.
            if (np.abs(delta)<10.).any() : delta=np.abs(delta)
            s12=np.sqrt(delta)
            ## export to pgf file
            fileName=path+'DPslow_yield0_s22s12_frame'+str(p)+'_Stress'+str(k)+'.pgf'
            #export2pgfPlotFile(fileName,np.array([s12,s22]),'sigma_12','sigma_22')
            #yields22_s12.append(fileName)
            ax2.plot(s22,s12,'k')
            if updated_criterion :
                # plot updated yield surface that takes into account hardening in both planes
                sig0=SIG11[-1,p,k]
                maxCrit=0.5*(sig0*(2.*nu**2-2.*nu-1.)+E*Epsp33[-1,p,k]*(1.-2.*nu))/(nu-nu**2-1.)
                Delta=(4.*(nu**2-nu+1.)*(sigy+H*plast)**2- 3.*(E*Epsp33[-1,p,k]+(1.-2.*nu)*sig0)**2)
                s22max=(sig0*(1.+2.*nu-2.*nu**2) +E*Epsp33[-1,p,k]*(2.*nu-1.) +np.sqrt(Delta))/(2.*(nu**2-nu+1.))
                s22min=(sig0*(1.+2.*nu-2.*nu**2)  +E*Epsp33[-1,p,k]*(2.*nu-1.) -np.sqrt(Delta))/(2.*(nu**2-nu+1.))
                s22=np.linspace(s22min,s22max,Samples)
                delta=(-(E*Epsp33[-1,p,k])**2 +E*Epsp33[-1,p,k]*(sig0+s22)*(2.*nu-1) + sig0*s22*(2.*nu+1.-2.*nu**2) + (sig0**2+s22**2)*(nu-nu**2-1.)+ (sigy+H*plast)**2)/3.
                if (np.abs(delta)<10.).any() : delta=np.abs(delta)
                s12=np.sqrt(delta)
                ax2.plot(s22,s12,color=col[p],linestyle='--')
                ax2.plot([maxCrit,maxCrit],[0.,tauM],color=col[p],linestyle='-.')
                limitX=np.array([maxCrit,maxCrit])
                limitY=np.array([0.,tauM])
                ## export to pgf file
                fileName=path+'DPslow_yieldfin_s22s12_frame'+str(p)+'_Stress'+str(k)+'.pgf'
                #export2pgfPlotFile(fileName,np.array([s12,s22]),'sigma_12','sigma_22')
                export2pgfPlotFile(fileName,np.array([limitY,limitX]),'sigma_12','sigma_22')
                yields22_s12.append(fileName)
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
                fileName=path+'DPslow_yieldfin_s11s12_frame'+str(p)+'_Stress'+str(k)+'.pgf'
                #export2pgfPlotFile(fileName,np.array([s12,s11]),'sigma_12','sigma_11')
        
                ax1.plot(s11,s12,color=col[p],linestyle='--')
                ax1.plot([maxCrit,maxCrit],[0.,tauM],color=col[p],linestyle='-.')
                limitX=np.array([maxCrit,maxCrit])
                limitY=np.array([0.,tauM])
                export2pgfPlotFile(fileName,np.array([limitY,limitX]),'sigma_12','sigma_11')
                yields11_s12.append(fileName)
                
            ax1.plot(SIG11[:,p,k],TAU[:,p,k],color=col[p],lw=2.5)
            ax2.plot(SIG22[:,p,k],TAU[:,p,k],color=col[p],lw=2.5)
            ax4.semilogy(LodeAngle_S[:,p,k],plast_S[:,p,k],color=col[p],lw=2.5)
            #ax4.plot(TAU[:,p,k],plast_S[:,p,k],color=col[p],lw=2.5)
            #ax3.plot(eigsigDevS[:,p,k,0],eigsigDevS[:,p,k,1],eigsigDevS[:,p,k,2],color=col[p],lw=2.5)
            ax3.plot(eigsigS[:,p,k,0],eigsigS[:,p,k,1],eigsigS[:,p,k,2],color=col[p],lw=2.5)
        ax3.plot([-sigy,sigy],[0.,0.],[0.,0.],color="k",linestyle="--",lw=1.)
        ax3.plot([0.,0.],[-sigy,sigy],[0.,0.],color="k",linestyle="--",lw=1.)
        ax3.plot([-radius,radius],[radius,-radius],[0.,0.],color="k",linestyle="--",lw=1.)
        plt.suptitle(r'Loading paths through slow waves for $\sigma_{22}$ ='+'{:.2e}'.format(sig22[k])+'Pa.', fontsize=24.)
        plt.show()

    ## sig22 value will change here
    xlabels=['$\sigma_{11} $','$\sigma_{22} $','$s_1 $'] #size=number of .tex files
    ylabels=['$\sigma_{12} $','$\sigma_{12} $','$s_2 $'] #size=number of .tex files
    zlabels=['','','$s_3$'] #size=number of .tex files


    subtitle=[r'(a) ($\sigma_{11},\sigma_{12}$) plane',r'(b) ($\sigma_{22},\sigma_{12}$) plane',r'(c) Deviatoric plane']

    srcX=['sigma_11','sigma_22']
    srcY=['sigma_12','sigma_12']

    name1='DPslowWaves_sig11_tau'+str(k)+'.tex'
    name2='DPslowWaves_sig22_tau'+str(k)+'.tex'
    name3='DPslowWaves_deviator'+str(k)+'.tex'
    names=[name1,name2,name3]
    
    files1=np.concatenate([pgfFilesList,yields11_s12])
    #files2=np.concatenate([pgfFilesList,yields22_s12])
    files2=pgfFilesList
    files2=np.concatenate([pgfFilesList,yields22_s12]) 
    pgfFiles=[files1,files2,deviatorPlots]
    #buildTeXFiles(names,pgfFiles,xlabels,ylabels,zlabels,subtitle,srcX,srcY)
    names=[[name1,name2],name3]
    pgfFiles=[[files1,files2],deviatorPlots]
    xlabels=[['$\sigma_{11} (Pa)$','$\sigma_{22}  (Pa)$'],'$s_1 $'] #size=number of .tex files
    ylabels=[['$\sigma_{12}  (Pa)$','$\sigma_{12}  (Pa)$'],'$s_2 $'] #size=number of .tex files
    zlabels=[['',''],'$s_3$'] #size=number of .tex files

    TauMax=1.1*np.max(TAU[0:-1:Niter/100,:,k])
    buildTeXFiles2(names,pgfFiles,xlabels,ylabels,zlabels,srcX,srcY,TauMax)
    
    pgfFilesList=[];yields11_s12=[];deviatorPlots=[];yields22_s12=[]

Samples=50#
# Sample constant stress component sig22
sig22=np.linspace(0.,sigy*np.sqrt(4*(nu**2-nu+1.))/np.sqrt(3.*(4.*nu**2-4.*nu+1.)),Samples)
sig22=np.linspace(-sigy*np.sqrt(4*(nu**2-nu+1.))/np.sqrt(3.*(4.*nu**2-4.*nu+1.)),sigy*np.sqrt(4*(nu**2-nu+1.))/np.sqrt(3.*(4.*nu**2-4.*nu+1.)),Samples)

#sig22=[0.]
sig=np.zeros((Samples,Samples))
tau=np.zeros((Samples,Samples))

fig = plt.figure(figsize=(10,10))
ax=plt.subplot2grid((1,1),(0,0),projection='3d')
ax.grid()

for k in range(len(sig22)):
    s22=sig22[k]
    
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
        if np.isnan(tau[i,k]): print "Nan ",delta,s11,s22,i,k
    ax.plot(sig[:,k],tau[:,k],s22,'k-.',lw=0.5)

ax.view_init(0.,0.)
#ax.set_xlabel(r'$\sigma_{11}$',size=24.)
ax.set_ylabel(r'$\sigma_{12}$',size=24.)
ax.set_zlabel(r'$\sigma_{22}$',size=24.)
ax.plot(SIG11[:,0,2],TAU[:,0,2],SIG22[:,0,2])                
# for k in range(len(sig[])):
#     s22=sig22[k]
# sig0=SIG11[-1,p,k]
# maxCrit=0.5*(sig0*(2.*nu**2-2.*nu-1.)+E*Epsp33[-1,p,k]*(1.-2.*nu))/(nu-nu**2-1.)
# Delta=(4.*(nu**2-nu+1.)*(sigy+H*plast)**2- 3.*(E*Epsp33[-1,p,k]+(1.-2.*nu)*sig0)**2)
# s22max=(sig0*(1.+2.*nu-2.*nu**2) +E*Epsp33[-1,p,k]*(2.*nu-1.) +np.sqrt(Delta))/(2.*(nu**2-nu+1.))
# s22min=(sig0*(1.+2.*nu-2.*nu**2)  +E*Epsp33[-1,p,k]*(2.*nu-1.) -np.sqrt(Delta))/(2.*(nu**2-nu+1.))
# s22=np.linspace(s22min,s22max,Samples)
# delta=(-(E*Epsp33[-1,p,k])**2 +E*Epsp33[-1,p,k]*(sig0+s22)*(2.*nu-1) + sig0*s22*(2.*nu+1.-2.*nu**2) + (sig0**2+s22**2)*(nu-nu**2-1.)+ (sigy+H*plast)**2)/3.
# if (np.abs(delta)<10.).any() : delta=np.abs(delta)
# s12=np.sqrt(delta)
# ax2.plot(s22,s12,color=col[p],linestyle='--')
                
# for i in range(SIG22.shape[2])[1:-1]:
#     for k in range(SIG22.shape[1]):
#         ax.plot(SIG11[:,k,i],TAU[:,k,i],SIG22[:,k,i],'r')

plt.show()

