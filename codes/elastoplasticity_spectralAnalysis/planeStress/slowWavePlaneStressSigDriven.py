# !\usr\bin\python
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import scipy.optimize
from matplotlib import animation
from scipy.integrate import ode
import pdb

# Material parameters
rho = 7800.
E = 2.e11
nu = 0.3
mu = 0.5*E/(1.+nu)
kappa = E/(3.*(1.-2.*nu))
lamb = kappa-2.*mu/3.
sigy = 100.0e6        
H = 100.08e6
beta=(6.*mu**2)/(3.*mu+H)

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
        H[0,2]=lamb-BETA*s11*s12
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

def computePsiSlow(sig11,sigma,sig33,lamb,mu,beta,tangent,rho):
    # sig11 driven
    n1=1.;n2=0.
    sig12=sigma[0];sig22=sigma[1]
    H=tangentModulus(np.array([sig11,sig12,sig22,sig33]),lamb,mu,beta,tangent)
    C=acousticTensor(H,np.array([n1,n2]))
    
    
    eigenf,eigens=acousticEigenStructure(C)
    alpha11=H[0,1]*H[1,2]- H[1,1]*H[0,2]
    alpha12=-H[0,1]*H[0,2]-H[0,0]*H[2,1]
    alpha22=H[0,0]*H[1,1]-H[0,1]**2
    w1=eigenf[1][0];w2=eigenf[1][1]
    psi12=-2.*w1/w2
    psi22=(2.*w1*alpha12/w2-alpha11)/alpha22
    """
    n1=1.;n2=0.
    JN=-np.array([[0.,0.,n1/rho,n2/rho,0.],[0.,0.,0.,n1/rho,n2/rho],[H[0,0]*n1+H[0,1]*n2,H[0,1]*n1+H[0,2]*n2,0.,0.,0.],[H[0,1]*n1+H[1,1]*n2,H[1,1]*n1+H[1,2]*n2,0,0,0],[H[2,0]*n1+H[2,1]*n2,H[2,1]*n1+H[2,2]*n2,0,0,0]])
    eigenStructure=np.linalg.eig(JN.T)
    contact=np.where(eigenStructure[0]==0)[0][0]
    cfplus=np.where(eigenStructure[0]==np.max(eigenStructure[0]))[0][0]
    cfminus=np.where(eigenStructure[0]==np.min(eigenStructure[0]))[0][0]
    index=np.ones(5);index[[contact,cfminus,cfplus]]-=1
    cs=np.where(index!=0.)[0]
    csminus=np.where(eigenStructure[0]==np.min(eigenStructure[0][cs]))[0][0]
    csplus=np.where(eigenStructure[0]==np.max(eigenStructure[0][cs]))[0][0]
    lcfminus=eigenStructure[1][:,cfminus];lcfplus=eigenStructure[1][:,cfplus]
    lcontact=eigenStructure[1][:,contact]
    dl=lcfminus-lcfplus

    if not (dl[4]!=0. and dl[0]!=0. and dl[1]!=0.):
        psi12=-dl[2]/dl[3]

    if not (lcontact[0]>1.e-6 and lcontact[1]>1.e-6):
        psi22=(lcontact[3]*(dl[2]/dl[3])-lcontact[2])/lcontact[4]
    """
    return np.array([psi12,psi22])

def computeLodeAngle(sig11,sig22,sig12,sig33):
    # deviatoric stress
    sDev=computeDeviatoricPart(np.array([sig11,sig12,sig22,sig33]))
    s11=sDev[0];s12=sDev[1]/np.sqrt(2.);s22=sDev[2];s33=sDev[3]
    # deviator 2nd and 3rd invariants
    J3=s33*(s11*s22-s12**2) ; sqrtJ2=np.sqrt(0.5*np.dot(sDev,sDev))
    theta=np.arccos((3./2.)*np.sqrt(3.)*J3/(sqrtJ2**3))/3.
    theta=theta*360./(2.*np.pi)
    return theta

def updateEquivalentPlasticStrain(sig,sign,H):
    # sig=[sig11^n , sqrt(2)*sig12^n , sig22 , sig33^n]
    # sign=[sig11^n+1 , sqrt(2)*sig12^n+1 , sig22 , sig33^n+1]
    sigDev=computeDeviatoricPart(np.array([sign[0],sign[1]/np.sqrt(2.),sign[2],sign[3]]))
    norm=np.sqrt(np.dot(sigDev,sigDev))
    flow=sigDev/norm
    dSig=sign-sig
    dp=(1./H)*np.sqrt(3./2.)*np.dot(flow,dSig)
    return dp

def plasticResidual(sig,sign,p,pn,H):
    # sig=[sig11^n , sqrt(2)*sig12^n , sig22 , sig33^n]
    # sign=[sig11^n+1 , sqrt(2)*sig12^n+1 , sig22 , sig33^n+1]
    sigDev=computeDeviatoricPart(np.array([sign[0],sign[1]/np.sqrt(2.),sign[2],sign[3]]))
    norm=np.sqrt(np.dot(sigDev,sigDev))
    flow=sigDev/norm
    dSig=sign-sig
    dp=(1./H)*np.sqrt(3./2.)*np.dot(flow,dSig)
    res=pn-p-dp
    return res

def computeEigenStresses(sig):
    #    | sig11 sig12   0   |
    #sig=| sig12 sig22   0   |
    #    |   0     0   sig33 |
    s3=sig[2,2]
    delta=(sig[0,0]-sig[1,1])**2+4.*sig[0,1]**2
    s1=0.5*(sig[0,0]+sig[1,1]-np.sqrt(delta))
    s2=0.5*(sig[0,0]+sig[1,1]+np.sqrt(delta))
    return np.array([s1,s2,s3])

from mpl_toolkits.mplot3d import proj3d
def orthogonal_proj(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    return np.array([[1,0,0,0],
                        [0,1,0,0],
                        [0,0,a,b],
                        [0,0,0,zback]])
proj3d.persp_transformation = orthogonal_proj

Samples=5

# Sample constant stress component sig22
sig22=np.linspace(0.,sigy,Samples)
#sig22=np.linspace(-sigy/np.sqrt(1-nu+nu**2),sigy/np.sqrt(1-nu+nu**2),Samples)

Samples*=10
sig=np.zeros((Samples,Samples))
tau=np.zeros((Samples,Samples))

frames=[10,20,40]
frames=[5,10,15,20]
col=["r","g","b","y","c","m","k","p"]
tauM=1.5*sigy/np.sqrt(3.)
sigM=1.5*sigy/np.sqrt(1-nu+nu**2)
tauM=sigM
Niter=1000
TAU=np.zeros((Niter,len(frames),len(sig22)))
SIG11=np.zeros((Niter,len(frames),len(sig22)))
SIG22=np.zeros((Niter,len(frames),len(sig22)))
eigsigS=np.zeros((Niter,len(frames),len(sig22),3))
criterionS=np.zeros((Niter,len(frames)))
PsiS=np.zeros((Samples,len(sig22)))

plast_S=np.zeros((Niter,len(frames)))
LodeAngle_S=np.zeros((Niter,len(frames)))
# Boolean to plot the upadted yield surface
updated_criterion=False
for k in range(len(sig22)-1):
    s22=sig22[k]
    
    Delta=(4.*sigy**2- 3.*s22**2)
    sigMax=(s22+np.sqrt(Delta))/2.
    sigMin=(s22-np.sqrt(Delta))/2.

    # Sample stress component sig11
    sig[:,k]=np.linspace(sigMin,sigMax,Samples)
    sig[:,k]=np.linspace(0.,sigMax,Samples)
    
    # Compute shear stress satisfying the criterion given sig11 and sig22
    for i in range(Samples):
        s11=sig[i,k]
        delta=(s11*s22 -s11**2-s22**2 + sigy**2)/3.
        if np.abs(delta)<10. : delta=np.abs(delta)
        tauMax=np.sqrt(delta)
        f_vm=lambda x:computeCriterion(s11,s22,x,0.,sigy)
        tau[i,k]=np.sqrt(delta)
        
    

## LOADING PATHS PLOTS
for k in range(len(sig22)-1)[1:]:
    s22=sig22[k]
    sigM=1.25*np.max(sig[:,k])
    tauM=1.25*np.max(tau[:,k])
    ## For each value of sig22 trace the loading paths given by psis from yield surface to an arbitrary shear stress level
    approx=np.zeros((len(frames),2))
    ordonnees=np.zeros((len(frames),Samples))
    abscisses=np.zeros((len(frames),Samples))
    radius_S=np.zeros(len(frames))
    for s,i in enumerate(frames):
        if i==0:
            continue
        sig0=sig[-1-i,k]
        tau0=tau[-1-i,k]

        
        dsig=(sigM-sig0)/Niter
        
        SIG11[:,s,k]=np.linspace(sig0,sigM,Niter)
        
        TAU[0,s,k]=tau0
        SIG22[0,s,k]=s22

        
        #rSlow = ode(computePsiSlow).set_integrator('vode',method='bdf')
        rSlow = ode(computePsiSlow).set_integrator('vode',method='adams',order=12)
        rSlow.set_initial_value(np.array([TAU[0,s,k],SIG22[0,s,k]]),SIG11[0,s,k]).set_f_params(0.,lamb,mu,beta,'planeStress',rho)
        sigma = np.matrix([[SIG11[0,s,k],TAU[0,s,k],0.],[TAU[0,s,k],SIG22[0,s,k],0.],[0.,0.,0.]])
        eigsig=np.linalg.eig(sigma)[0]
        eigsigS[0,s,k,:]=eigsig
        LodeAngle_S[0,s]=computeLodeAngle(sigma[0,0],SIG22[0,s,k],sigma[0,1],0.)
            
        p=0.
        epsp33=0.
        for j in range(Niter-1):
            rSlow.set_f_params(np.array([TAU[j,s,k],SIG22[j,s,k]]),0.,lamb,mu,beta,'planeStress',rho)
            if not rSlow.successful():
                print "Integration issues in slow wave path"
                break
            rSlow.integrate(rSlow.t+dsig)
        
            TAU[j+1,s,k],SIG22[j+1,s,k]=rSlow.y
            
            sigma = np.array([SIG11[j,s,k],np.sqrt(2.)*TAU[j,s,k],SIG22[j,s,k],0.])
            sigman = np.array([SIG11[j+1,s,k],np.sqrt(2.)*TAU[j+1,s,k],SIG22[j+1,s,k],0.])
            f_vm=computeCriterion(SIG11[j+1,s,k],SIG22[j+1,s,k],TAU[j+1,s,k],0.,sigy+H*p)
            #if f_vm>0. : 
            #p+=updateEquivalentPlasticStrain(sigma,sigman,H)
            #residual=lambda x: plasticResidual(sigma,sigman,p,x,H)
            residual=lambda x: computeCriterion(SIG11[j+1,s,k],SIG22[j+1,s,k],TAU[j+1,s,k],0.,sigy+H*x)
            p=scipy.optimize.root(residual,p,method='hybr',options={'xtol':1.e-12}).x[0]
            criterionS[j+1,s]=computeCriterion(SIG11[j+1,s,k],SIG22[j+1,s,k],TAU[j+1,s,k],0.,sigy+H*p)
            plast_S[j+1,s]=p
            LodeAngle_S[j+1,s]=computeLodeAngle(sigman[0],sigman[2],sigman[1]/np.sqrt(2.),0.)
            
            # Eigenvalues of sigma (for deviatoric plane plots)
            sigma = np.matrix([[SIG11[j+1,s,k],TAU[j+1,s,k],0.],[TAU[j+1,s,k],SIG22[j+1,s,k],0.],[0.,0.,0.]])
            eigsigS[j+1,s,k,:]=computeEigenStresses(sigma)
            
        print "Final equivalent plastic strain after slow wave : ",p

        radius_S[s]=sigy+H*p
    TAU_MAX_S=np.max(ordonnees)
    SIG_MAX_S=np.max(abscisses)

    ### SUBPLOTS SETTINGS
    fig = plt.figure()
    ax2=plt.subplot2grid((1,2),(0,1),projection='3d')
    ax1d1=plt.subplot2grid((1,2),(0,0))
    ax1d1.grid()
    ax1d1.set_xlabel(r'$\Theta$', fontsize=24)
    ax1d1.set_ylabel('p', fontsize=24)
    fvm1=ax1d1.twinx()
    fvm1.set_ylabel('f',fontsize=18.)
    fvm1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    cylindre=vonMisesYieldSurface(sigy)
    ax2.plot_wireframe(cylindre[0,:],cylindre[1,:],cylindre[2,:], color="k")
    elevation_Angle_radian=np.arctan(1./np.sqrt(2.0))
    angle_degree= 180.*elevation_Angle_radian/np.pi
    radius=1.*np.sqrt((2./3.)*sigy**2)
    
    ax2.set_xlim(-1.*radius,1.*radius)
    ax2.set_ylim(-1.*radius,1.*radius)
    ax2.set_zlim(-1.*radius,1.*radius)
    ax2.view_init(angle_degree,45.)
    ax2.plot([0.,sigy],[0.,sigy],[0.,sigy],color="k")
    ax2.set_xlabel(r'$\sigma_1$',size=24.)
    ax2.set_ylabel(r'$\sigma_2$',size=24.)
    ax2.set_zlabel(r'$\sigma_3$',size=24.)
    
    
    for p in range(len(frames)):
        if updated_criterion :
            cylindre=vonMisesYieldSurface(radius_S[p])
            ax2.plot_wireframe(cylindre[0,:],cylindre[1,:],cylindre[2,:], color=col[p],linestyle='--')
        ## 2D plot of equivalent plastic strain evolution
        ax1d1.plot(LodeAngle_S[:Niter/5,p],plast_S[:Niter/5,p],col[p])
        #ax1d1_2.plot(LodeAngle_S[:Niter/5,p],SIG33_S[:Niter/5,p,k],col[p],marker='o')
        fvm1.plot(LodeAngle_S[:,p],criterionS[:,p],col[p],linestyle='--')
        ## 3D plots of loading paths (deviatoric plane)
        ax2.plot(eigsigS[:,p,k,0],eigsigS[:,p,k,1],eigsigS[:,p,k,2],color=col[p],marker="o")
    ax2.plot([-sigy,sigy],[0.,0.],[0.,0.],color="k",linestyle="--",lw=1.)
    ax2.plot([0.,0.],[-sigy,sigy],[0.,0.],color="k",linestyle="--",lw=1.)
    ax2.plot([-radius,radius],[radius,-radius],[0.,0.],color="k",linestyle="--",lw=1.)
    
    #plt.show()
    
    fig = plt.figure()
    ax1=plt.subplot2grid((1,2),(0,0))
    ax2=plt.subplot2grid((1,2),(0,1))
    ax1.set_xlabel(r'$\sigma_{11}$',size=28.)
    ax1.set_ylabel(r'$\sigma_{12}$',size=28.)
    #ax1.set_zlabel(r'$\sigma_{22}$',size=28.)
    ax2.set_xlabel(r'$\sigma_{22}$',size=28.)
    ax2.set_ylabel(r'$\sigma_{12}$',size=28.)
    #ax2.set_zlabel(r'$\sigma_{11}$',size=28.)
    ax1.grid()
    ax2.grid()
    #ax2.view_init(-90.,-0.)
    #ax1.view_init(-90.,0.)
    for s,i in enumerate(frames):
        sig0=sig[-1-i,k]
        s22max=(sig0+np.sqrt(4*sigy**2-3.*sig0**2))/2.
        s22min=(sig0-np.sqrt(4*sigy**2-3.*sig0**2))/2.
        s22=np.linspace(s22min,s22max,Samples)
        s12=np.sqrt((sigy**2- sig0**2-s22**2+sig0*s22)/3.)
        ax2.plot(s22,s12,color=col[s])
    ax1.plot(sig[:,k],tau[:,k],'k')
    #ax2.plot(sig[:,k],tau[:,k],sig22[k],'k')
    
    for p in range(len(frames)):
        ax1.plot(SIG11[:,p,k],TAU[:,p,k],color=col[p])
        ax2.plot(SIG22[:,p,k],TAU[:,p,k],color=col[p])
    plt.show()
