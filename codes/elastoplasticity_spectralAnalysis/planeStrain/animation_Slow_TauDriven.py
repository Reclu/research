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
H = 100.0e8
# sigy = 400.0e6        
# H = 100.0e8                     
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

############ SIG DRIVEN
def loadingPath_Sig(dtau,sigma,sigman,E,nu,h,lamb,mu,beta,tangent):
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
    psi12=-w1/w2
    psi22=(w1*alpha12/w2-alpha11)/alpha22
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
    #return np.array([psi11*dtau,psi22*dtau,depsp33])
    return np.array([psi12*dsig,psi22*dsig,nu*(1.+psi22)*dsig-E*depsp33,depsp33])
    
def integrateODE_Sig(dsig,sig0,tau0,sig22_0,sig33,epsp33,nu,E,H,lamb,mu,beta,tangent):
    sigma=np.array([tau0,sig22_0,sig33,epsp33])
    # computePsiSlow(sig12,sigma,sig33,lamb,mu,beta,tangent)
    # subdivision of time step
    sub_steps = 1
    dSIG = dsig/sub_steps
    theta = 1.0
    for i in range(sub_steps):
        ## Nonlinear solution procedure
        ## R = s^{n+1} - s^{n} - RHS
        
        R=lambda x: x - sigma - theta*loadingPath_Sig(dSIG,np.array([sig0,sigma[0],sigma[1],sigma[2]]),np.array([sig0+dSIG,x[0],x[1],x[2]]),E,nu,H,lamb,mu,beta,tangent)
        #pdb.set_trace()
        
        #solution = scipy.optimize.fsolve(R,np.array([sigma[0],sigma[1],sigma[2],sigma[3]]))
        solution = scipy.optimize.root(R,np.array([sigma[0],sigma[1],sigma[2],sigma[3]])).x
        sigma = solution
    return solution[0],solution[1],solution[2],solution[3]
#######
####### TAU DRIVEN
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
    # computePsiSlow(sig12,sigma,sig33,lamb,mu,beta,tangent)
    # subdivision of time step
    sub_steps = 1
    dTAU = dtau/sub_steps
    theta = 1.0
    for i in range(sub_steps):
        ## Nonlinear solution procedure
        ## R = s^{n+1} - s^{n} - RHS
        
        R=lambda x: x - sigma - theta*loadingPath(dTAU,np.array([sigma[0],tau0,sigma[1],sigma[2]]),np.array([x[0],tau0+dTAU,x[1],x[2]]),E,nu,H,lamb,mu,beta,tangent)
        #pdb.set_trace()
        
        solution = scipy.optimize.fsolve(R,np.array([sigma[0],sigma[1],sigma[2],sigma[3]]))
        sigma = solution
    return solution[0],solution[1],solution[2],solution[3]
##################

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

frames=[0,5,20,30,40]
frames=[5]
col=["r","g","b","y","c","m","k","p"]
# purple to red
col=['#781C81','#3F60AE','#539EB6','#6DB388','#CAB843','#E78532','#D92120']
tauM=1.5*sigy/np.sqrt(3.)
sigM=1.5*sigy/np.sqrt(1-nu+nu**2)
tauM=sigM
Niter=500
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
iteration=np.zeros((len(frames),len(sig22)),dtype=int)
Ranges=np.zeros((len(frames),len(sig22)),dtype=int)
# Boolean to plot the upadted yield surface
updated_criterion=True
chain_waves=True
chaining_parameter='sig'
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
        tauM=1.5*np.max(tau[:,k])
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
            SIG11[j+1,s,k],SIG22[j+1,s,k],SIG33[j+1,s,k],epsp33=integrateODE(dtau,SIG11[j,s,k],TAU[j,s,k],SIG22[j,s,k],SIG33[j,s,k],epsp33,nu,E,H,lamb,mu,beta,tangent)
            sig0=SIG11[j+1,s,k];
            MaxCrit[j+1,s,k]=0.5*(sig0*(2.*nu**2-2.*nu-1.)+E*epsp33*(1.-2.*nu))/(nu-nu**2-1.)
                
            sigma = np.array([SIG11[j,s,k],np.sqrt(2.)*TAU[j,s,k],SIG22[j,s,k],SIG33[j,s,k]])
            
            sig33=SIG33[j+1,s,k]
            
            sigman=np.array([SIG11[j+1,s,k],TAU[j+1,s,k]*np.sqrt(2.),SIG22[j+1,s,k],sig33])
            
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
            
            Nfinal=Niter
            iteration[s,k]=Niter-1
            if rcs2[j+1,s,k]>rcs2[j,s,k]:
                print "Simple wave condition violated"
                iteration[s,k]=int(j)
                Nfinal=int(j)
                break
            #pdb.set_trace()
        if chain_waves and Nfinal!=Niter:
            print "Trying to change the driving parameter (",chaining_parameter,")"
            #pdb.set_trace()
            ## If numerical issues by driving with tau, try to carry on the integration by driving with sigma
            if chaining_parameter=='sig':
                sig0=SIG11[Nfinal,s,k]
                sigM=1.05*sig0
                dsig=(sigM-sig0)/(Niter-Nfinal-1)
                SIG11[Nfinal:,s,k]=np.linspace(sig0,sigM,Niter-Nfinal)
            elif chaining_parameter=='tau':
                tau0=TAU[Nfinal,s,k]
                tauM=0.8*tau0
                dtau=(tauM-tau0)/(Niter-Nfinal-1)
                TAU[Nfinal:,s,k]=np.linspace(tau0,tauM,Niter-Nfinal)
                
            for j in range(Niter-1)[Nfinal:]:
                if chaining_parameter=='sig':
                    TAU[j+1,s,k],SIG22[j+1,s,k],SIG33[j+1,s,k],epsp33=integrateODE_Sig(dsig,SIG11[j,s,k],TAU[j,s,k],SIG22[j,s,k],SIG33[j,s,k],epsp33,nu,E,H,lamb,mu,beta,tangent)
                elif chaining_parameter=='tau':
                    SIG11[j+1,s,k],SIG22[j+1,s,k],SIG33[j+1,s,k],epsp33=integrateODE(dtau,SIG11[j,s,k],TAU[j,s,k],SIG22[j,s,k],SIG33[j,s,k],epsp33,nu,E,H,lamb,mu,beta,tangent)
                sigma = np.array([SIG11[j,s,k],np.sqrt(2.)*TAU[j,s,k],SIG22[j,s,k],SIG33[j,s,k]])
            
                sig33=SIG33[j+1,s,k]
            
                sigman=np.array([SIG11[j+1,s,k],TAU[j+1,s,k]*np.sqrt(2.),SIG22[j+1,s,k],sig33])
                dp=updateEquivalentPlasticStrain(sigma,sigman,H)
                plast+=dp
                
                #criterionS[j+1,s,k]=computeCriterion(SIG11[j+1,s,k],SIG22[j+1,s,k],TAU[j+1,s,k],sig33,sigy+H*plast)
                plast_S[j+1,s,k]=plast
                
                Epsp33[j+1,s,k]=epsp33
                LodeAngle_S[j+1,s,k]=computeLodeAngle(sigman[0],sigman[2],sigman[1]/np.sqrt(2.),sig33)
                
                # Eigenvalues of sigma (for deviatoric plane plots)
                sigma = np.matrix([[SIG11[j+1,s,k],TAU[j+1,s,k],0.],[TAU[j+1,s,k],SIG22[j+1,s,k],0.],[0.,0.,SIG33[j+1,s,k]]])
                rcs2[j+1,s,k]=np.sqrt(computeSpeed(sigma,lamb,mu,beta,tangent)/rho)

                Eps[:,j+1,s,k]=updateEpsilon(Eps[:,j,s,k],np.array([SIG11[j,s,k],TAU[j,s,k],SIG22[j,s,k],SIG33[j,s,k]]),np.array([SIG11[j+1,s,k],TAU[j+1,s,k],SIG22[j+1,s,k],SIG33[j+1,s,k]]),E,nu,H)
            
                eigsigS[j+1,s,k,:]=computeEigenStresses(sigma)
                # sigDev=computeDeviatoricPart(np.array([SIG11[j+1,s,k],TAU[j+1,s,k],SIG22[j+1,s,k],SIG33[j+1,s,k]]))
                # sigma = np.matrix([[sigDev[0],sigDev[1]/np.sqrt(2.),0.],[sigDev[1]/np.sqrt(2.),sigDev[2],0.],[0.,0.,sigDev[3]]])
                eigsigDevS[j+1,s,k,:]=computeEigenStresses(sigma)
            
                Nfinal=Niter
                iteration[s,k]=Niter
                if rcs2[j+1,s,k]>rcs2[j,s,k]:
                    print "Simple wave condition violated"
                    Nfinal=int(j)
                    iteration[s,k]=int(j)
                    break
        
        time=np.linspace(0,j+1,Niter)
        print "Final equivalent plastic strain after slow wave : ",plast
        fileName=path+'DPslowStressPlane_frame'+str(s)+'_Stress'+str(k)+'.pgf'
        ## color bar of rcs2
        if Nfinal/100>1:
            pas=Nfinal/100
            ranging=np.linspace(0,Nfinal-1,Nfinal/100,True,False,np.int)
        else:
            pas=1
            ranging=np.linspace(0,Nfinal-1,Nfinal,True,False,np.int)
        #pdb.set_trace()
        
        radius_S[s]=sigy+H*plast
          
    
    cylindre=vonMisesYieldSurface(sigy)
    fileName=path+'CylindreDevPlane.pgf'
    export2pgfPlot3D(str(fileName),cylindre[0,:],cylindre[1,:],cylindre[2,:])
    deviatorPlots.append(fileName)
    TAU_MAX_S=np.max(ordonnees)
    SIG_MAX_S=np.max(abscisses)

    
# Bluid animated loading paths for one initial value of sig22
## Evolution of the driving quantity on the left
## Loading path in deviatoric plane on the right
fig, (ax1, ax2) = plt.subplots(2,1)
# fig = plt.figure(figsize=(10,10))
ax1=plt.subplot2grid((1,2),(0,0))
ax2=plt.subplot2grid((1,2),(0,1),projection='3d')
# fig, ax1 = plt.subplots(1,1)

cylindre=vonMisesYieldSurface(sigy)
ax2.plot(cylindre[0,:],cylindre[1,:],cylindre[2,:], color="k",lw=0.5)
elevation_Angle_radian=np.arctan(1./np.sqrt(2.0))
angle_degree= 180.*elevation_Angle_radian/np.pi
radius=1.*np.sqrt((2./3.)*sigy**2)
ax2.set_xlim(-1.*radius,1.*radius)
ax2.set_ylim(-1.*radius,1.*radius)
ax2.set_zlim(-1.*radius,1.*radius)
ax2.view_init(angle_degree,45.)
ax2.plot([-sigy,sigy],[0.,0.],[0.,0.],color="k",linestyle="--",lw=1.)
ax2.plot([0.,0.],[-sigy,sigy],[0.,0.],color="k",linestyle="--",lw=1.)
ax2.plot([-radius,radius],[radius,-radius],[0.,0.],color="k",linestyle="--",lw=1.)

sig22=2
#line1, = ax1.plot([], [],col[0] )
# line2, = ax1.plot([], [],col[1])
# line3, = ax1.plot([], [],col[2])
# line=[line1,line2,line3]
ax1.set_xlim(0,Niter)
ax1.set_ylim(0.,np.max(TAU[:,0,sig22]))
ax1.grid()

line3D=[]
line=[]
for i,s in enumerate(frames):
    line1D, = ax1.plot([], [],col[i] )
    line.append(line1D)
    line2D, = ax2.plot(eigsigS[0:1,0,i,0],eigsigS[0:1,0,i,1],eigsigS[0:1,0,i,2],col[i])
    line3D.append(line2D)

def init():
    for k,s in enumerate(frames):
        line[k].set_data([], [])
    return line

def animate(i):
    
    for k,s in enumerate(frames):
        line[k].set_data(np.arange(0,i,1),TAU[:i,k,sig22])
        line3D[k].set_data(eigsigS[0:i,k,sig22,0],eigsigS[0:i,k,sig22,1])
        line3D[k].set_3d_properties(eigsigS[0:i,k,sig22,2])
    return line,line3D

video_length=10#seconds
frames_per_second = int(Niter/video_length)
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=Niter, interval=100, blit=False)
Writer = animation.writers['ffmpeg']
writers = Writer(fps=frames_per_second, metadata=dict(artist='Me'), bitrate=-1)
# plt.rcParams['animation.ffmpeg_path']='C:/ffmpeg/bin/ffmpeg.exe'
# writers=animation.FFMpegWriter(bitrate=500)
anim.save('loadingpath_slow.mp4', writer=writers)
#anim.save('loadingpath_slow.mp4', extra_args=['-vcodec', 'libx264'])

plt.show()
