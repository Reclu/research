#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from scipy.optimize import curve_fit
import pdb
from matplotlib import rcParams

rcParams['axes.labelsize'] = 24
rcParams['xtick.labelsize'] = 24
rcParams['ytick.labelsize'] = 24
rcParams['legend.fontsize'] = 16



def computeLpNorm(Unum,Uexact,dx,p):
    return (((dx*((np.abs(Unum-Uexact))**p)).sum())**(1.0/p))

def computeRelativeError(Unum,Uexact,dx,p):
    return (computeLpNorm(Unum,Uexact,dx,p)/computeLpNorm(np.zeros(len(Unum)),Uexact,dx,p))

def export2DTeXFile(fileName,xField,xlabel,ylabel,fields,*kwargs):
    TeXFile=open(fileName,"w")
    n_fields = np.shape(fields)[0]
    n_labels = np.shape(kwargs)[0]
    # Define Paul Tol's colors (purple to red)
    color=['Purple','Blue','Duck','Green','Yellow','Orange','Red']
    col=['120,28,129','63,96,174','83,158,182','109,179,136','202,184,67','231,133,50','217,33,32']
    marker=['*','x','triangle*','square*','+','star','pentagone*']
    for i in range(len(col)):
        TeXFile.write(r'\definecolor{'+color[i]+'}{RGB}{'+col[i]+'}')
        TeXFile.write('\n')
    TeXFile.write(r'\begin{tikzpicture}[scale=0.7]');TeXFile.write('\n')
    TeXFile.write(r'\begin{loglogaxis}[xlabel='+str(xlabel)+',ylabel='+str(ylabel)+',ymajorgrids=true,xmajorgrids=true,legend pos=south east,title={() -- particles per cell},xmin=0.07,xmax=2.,ymin=0.001,ymax=0.2]');TeXFile.write('\n')
    legend=''
    TeXFile.write(r'\addplot[Blue,very thick,mark=*] coordinates {')
    for j in range(len(fields[0,:])):
        TeXFile.write('('+str(xField[j])+','+str(fields[0,j])+') ')
    TeXFile.write('};\n')
    TeXFile.write(r'\addplot[Purple,very thick,mark=triangle*] coordinates {')
    for j in range(len(fields[3,:])):
        TeXFile.write('('+str(xField[j])+','+str(fields[3,j])+') ')
    TeXFile.write('};\n')
    TeXFile.write(r'\addplot[Red,very thick,mark=+] coordinates {')
    for j in range(len(fields[1,:])):
        TeXFile.write('('+str(xField[j])+','+str(fields[1,j])+') ')
    TeXFile.write('};\n')
    TeXFile.write(r'\addplot[Orange,very thick,mark=square] coordinates {')
    for j in range(len(fields[2,:])):
        TeXFile.write('('+str(xField[j])+','+str(fields[2,j])+') ')
    TeXFile.write('};\n')
    
    
    TeXFile.write(r'\legend{dgmpm (Euler),dgmpm (RK2),usl,usf}')
    TeXFile.write('\n')    
    TeXFile.write(r'\draw (axis cs:0.2,0.1) -- (axis cs:0.2/1.4,0.1/1.4);')
    TeXFile.write('\n')    
    TeXFile.write(r'\draw (axis cs:0.2,0.1) -- (axis cs:0.2,0.1/1.4) node [midway,right] {\scriptsize 1};')
    TeXFile.write('\n')    
    TeXFile.write(r'\draw (axis cs:0.2,0.1/1.4) -- (axis cs:0.2/1.4,0.1/1.4) node [midway,below] {\scriptsize 1};')
    TeXFile.write('\n')    
    
    TeXFile.write(r'\end{loglogaxis}')
    TeXFile.write('\n')
    TeXFile.write('\end{tikzpicture}')
    TeXFile.write('\n')
    TeXFile.close()

###Opening the files and computation of the solution by each method
###Parameters####
MP = [20,40,60,80,100,120,140,160,200,300]
ppc=[2,3,4,8]


MP=[4,8,16,32,64]
ppc=[4,6,8]
# MP = [10,20,40,60,80]
# ppc=[2,3]

#MP = [120,160,200,240,300]
loading=-75.
dx = np.zeros((len(MP),len(ppc)))
DX = np.zeros(len(MP))

ErrS_dgmpmRK2 = np.zeros((len(MP),len(ppc)))
ErrS_dgmpm = np.zeros((len(MP),len(ppc)))
ErrS_dgmpm2 = np.zeros((len(MP),len(ppc)))
ErrS_usl = np.zeros((len(MP),len(ppc)))
ErrS_usf = np.zeros((len(MP),len(ppc)))

ErrV_dgmpmRK2 = np.zeros((len(MP),len(ppc)))
ErrV_dgmpm = np.zeros((len(MP),len(ppc)))
ErrV_dgmpm2 = np.zeros((len(MP),len(ppc)))
ErrV_usl = np.zeros((len(MP),len(ppc)))
ErrV_usf = np.zeros((len(MP),len(ppc)))

Error_order1=2
Error_order2=2
period=0.4
for p in range(len(ppc)):
    print "==============================",ppc[p],"PPC =============================="
    for i,Ne in enumerate(MP):
        parameters = {"Mp":Ne,"ppc":ppc[p],"loading":loading,"t_order":1,"period":period}
        ##DGMPM: Discontinuous Galerkin Material Point Method
        DGMPM = dict(parameters)
        execfile('dgmpm.py', DGMPM)
        parameters = {"Mp":Ne,"ppc":ppc[p],"loading":loading,"t_order":2,"period":period,"alg":'USL'}
        ##DGMPM: Discontinuous Galerkin Material Point Method with RK2 time integration
        DGMPMRK2 = dict(parameters)
        execfile('dgmpm.py', DGMPMRK2)
        USL = dict(parameters)
        execfile('mpm.py', USL)
        parameters = {"Mp":Ne,"ppc":ppc[p],"loading":loading,"t_order":2,"period":period,"alg":'USF'}
        USF = dict(parameters)
        execfile('mpm.py', USF)
        ########## Compute Errors #################
        dx[i,p]=DGMPM["dx"]
        n1= DGMPM["Increments"]
        n2= USL["Increments"]
        n1RK2= DGMPMRK2["Increments"]
        #Sigma
        ErrS_dgmpm[i,p] = computeRelativeError(DGMPM["Stress"][:,n1],DGMPM["Sth"][:,n1],dx[i,p],Error_order1)
        ErrS_dgmpmRK2[i,p] = computeRelativeError(DGMPMRK2["Stress"][:,n1RK2],DGMPMRK2["Sth"][:,n1RK2],dx[i,p],Error_order1)
        ErrS_usl[i,p] = computeRelativeError(USL["Stress"][:,n2],USL["Sth"][:,n2],dx[i,p],Error_order1)
        ErrS_usf[i,p] = computeRelativeError(USF["Stress"][:,n2],USF["Sth"][:,n2],dx[i,p],Error_order1)
        ErrS_dgmpm2[i,p] = computeRelativeError(DGMPM["Stress"][:,n1],DGMPM["Sth"][:,n1],dx[i,p],Error_order2)

        #Velocity
        ErrV_dgmpm[i,p] = computeRelativeError(DGMPM["Velocity"][:,n1],DGMPM["Vth"][:,n1],dx[i,p],Error_order1)
        ErrV_dgmpmRK2[i,p] = computeRelativeError(DGMPMRK2["Velocity"][:,n1RK2],DGMPMRK2["Vth"][:,n1RK2],dx[i,p],Error_order1)
        ErrV_usl[i,p] = computeRelativeError(USL["Velocity"][:,n2],USL["Vth"][:,n2],dx[i,p],Error_order1)
        ErrV_usf[i,p] = computeRelativeError(USF["Velocity"][:,n2],USF["Vth"][:,n2],dx[i,p],Error_order1)
        ErrV_dgmpm2[i,p] = computeRelativeError(DGMPM["Velocity"][:,n1],DGMPM["Vth"][:,n1],dx[i,p],Error_order2)
        #########################  Comparison  ######################################
        position=np.linspace(0.,1.,10)
        
        E=1.e7;rho=1000.;c=np.sqrt(E/rho);G=1.e-4
        reference=E*G*np.pi*np.sin(np.pi*c*0.02)*np.cos(np.pi*position)
        
        ##Stress
        l1,l2,l3,l4,l5=plt.plot(DGMPM["Pos"][:,n1],DGMPM["Stress"][:,n1],'r+', \
                                DGMPM["Pos"][:,n1],DGMPM["Sth"][:,n1],'k', \
                                DGMPMRK2["Pos"][:,n1RK2],DGMPMRK2["Stress"][:,n1RK2],'bo', \
                                USF["Pos"][:,n2],USF["Stress"][:,n2],'g^',\
                                position,reference,'r',linestyle='--')
        # if p==0 :
        #     plt.plot(FEM["y"],FEM["sigma"][:,n3],'b-o',label='fem')
        #     plt.plot(FEM["y"],FEM["anal"][:,n3],'k--')
        plt.grid()
        plt.xlabel('x (m)')
        plt.ylabel(r'$\sigma$')
        plt.show()
        # l1,l2,l3,l4=plt.plot(DGMPM["Pos"][:,n1],DGMPM["Stress"][:,n1],'r+', \
        #                      DGMPM["Pos"][:,n1],DGMPM["Sth"][:,n1],'k', \
        #                      DGMPMRK2["Pos"][:,n1RK2],DGMPMRK2["Stress"][:,n1RK2],'bo', \
        #                      USF["Pos"][:,n2],USF["Stress"][:,n2],'g^')
        # # if p==0 :
        # #     plt.plot(FEM["y"],FEM["sigma"][:,n3],'b-o',label='fem')
        # #     plt.plot(FEM["y"],FEM["anal"][:,n3],'k--')
        # plt.grid()
        # plt.xlabel('x (m)')
        # plt.ylabel(r'$\sigma$')
        # plt.show()
        #############################################################################
    export2DTeXFile('dgmpm_mpm_accuracyS_'+str(ppc[p])+'ppc.tex',dx[:,p],r'$\Delta X (m)$',r'$\epsilon_\sigma$',np.array([ErrS_dgmpm[:,p],ErrS_usl[:,p],ErrS_usf[:,p],ErrS_dgmpmRK2[:,p]]),['dgmpm (Euler)','usl','usf','dgmpm (RK2)'])
    export2DTeXFile('dgmpm_mpm_accuracyV_'+str(ppc[p])+'ppc.tex',dx[:,p],r'$\Delta X$',r'$\epsilon_v$',np.array([ErrV_dgmpm[:,p],ErrV_usl[:,p],ErrV_usf[:,p],ErrV_dgmpmRK2[:,p]]),['dgmpm (Euler)','usl','usf','dgmpm (RK2)'])
         
###########Plot convergence curves####################
#####Assess the constant and slopes of the
#####convergence curves by fitting a power law with curve_fit
def func(x,C,n):
    return (C*(x**n))


# errorList=[]
# for p in range(len(ppc)):
#     errorList.append(ErrS_dgmpm[:,p])
#     errorList.append(ErrS_usl[:,p])
#     errorList.append(ErrS_usf[:,p])
#     errorList.append(ErrS_dgmpmRK2[:,p])
#     errorList.append(ErrS_dgmpm2[:,p])

# popt = np.zeros((len(errorList),2))
# for i,convCurve in enumerate(errorList):
#     popt[i,:], pcov = curve_fit(func,dx[:,i/4],convCurve)


# print 'Taux de convergence dgmpm: '
# for i in range(len(ppc)):
#     print str(ppc[i]),"ppc : ",popt[2*i,1]
    
# print 'Taux de convergence mpm: '
# for i in range(len(ppc)):
#     print str(ppc[i]),"ppc : ",popt[2*i+1,1]


# errorList=[Err2_S_fem]
# popt2 = np.zeros((len(errorList),2))
# popt2, pcov2 = curve_fit(func,DX,Err2_S_fem)
# print 'Taux de convergence fem: ',popt2[1]

#pdb.set_trace()
#####Plot of convergence curves#########
fig1 = plt.figure()
#Stress
color_dgmpm=['r-+','y-+','c-+','g-+','b-+','m-+']
color_mpm=['y-o','g-o','c-o','r-o','b-o','m-o']
color_dgmpmRK2=['m-+','b-+','g-+','c-+','y-+','r-+']
color_dgmpm2=['m--','b--','g--','c--','y--','r--']
for p in range(len(ppc)):
    plt.loglog(dx[:,p],ErrS_dgmpm[:,p],color_dgmpm[p],lw=2.5,label='dgmpm')
    # plt.loglog(dx[:,p],func(dx[:,p],popt[4*p,0],popt[4*p,1]),'k-',label='dgmpm')

    plt.loglog(dx[:,p],ErrS_usl[:,p],color_mpm[p],lw=2.5,label='usl')
    # plt.loglog(dx[:,p],func(dx[:,p],popt[4*p+1,0],popt[4*p+1,1]),'k-')

    plt.loglog(dx[:,p],ErrS_usl[:,p],color_mpm[p],linestyle='-.',lw=2.5,label='usf')
    # plt.loglog(dx[:,p],func(dx[:,p],popt[4*p+1,0],popt[4*p+1,1]),'m-')

    plt.loglog(dx[:,p],ErrS_dgmpmRK2[:,p],color_dgmpmRK2[p],lw=2.5,label='dgmpm RK2')
    # plt.loglog(dx[:,p],func(dx[:,p],popt[4*p+2,0],popt[4*p+2,1]),'k-')


plt.loglog(dx[:,0],func(dx[:,0],1.,1.),'k-')
plt.loglog(dx[:,0],func(dx[:,0],1.,2.),'k--')

plt.grid(True,which="both")
plt.xlabel('grid size')
plt.legend(numpoints=1,loc='best')
plt.ylabel(r'$\frac{||\sigma- \sigma_{exact}||_2}{||\sigma_{exact}||_2}$')
plt.show()


#####Plot of convergence curves#########
fig1 = plt.figure()
#Velocity
color_dgmpm=['r-+','y-+','c-+','g-+','b-+','m-+']
color_mpm=['r--+','y--+','c--+','g--+','b--+','m--+']
color_dgmpmRK2=['m-+','b-+','g-+','c-+','y-+','r-+']
color_dgmpm2=['m--','b--','g--','c--','y--','r--']
for p in range(len(ppc)):
    plt.loglog(dx[:,p],ErrV_dgmpm[:,p],color_dgmpm[p],lw=2.5,label='dgmpm')
    # plt.loglog(dx[:,p],func(dx[:,p],popt[4*p,0],popt[4*p,1]),'k-',label='dgmpm')

    plt.loglog(dx[:,p],ErrV_usl[:,p],color_mpm[p],lw=2.5,label='usl')
    # plt.loglog(dx[:,p],func(dx[:,p],popt[4*p+1,0],popt[4*p+1,1]),'k-')

    plt.loglog(dx[:,p],ErrV_usl[:,p],color_mpm[p],linestyle='-.',lw=2.5,label='usf')
    # plt.loglog(dx[:,p],func(dx[:,p],popt[4*p+1,0],popt[4*p+1,1]),'m-')

    plt.loglog(dx[:,p],ErrV_dgmpmRK2[:,p],color_dgmpmRK2[p],lw=2.5,label='dgmpm RK2')
    # plt.loglog(dx[:,p],func(dx[:,p],popt[4*p+2,0],popt[4*p+2,1]),'k-')

plt.loglog(dx[:,0],func(dx[:,0],1.,1.),'k-')
plt.loglog(dx[:,0],func(dx[:,0],1.,2.),'k--')

plt.grid(True,which="both")
plt.xlabel('grid size')
plt.legend(numpoints=1,loc='best')
plt.ylabel(r'$\frac{||\sigma- \sigma_{exact}||_2}{||\sigma_{exact}||_2}$')
plt.show()
