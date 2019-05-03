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


rho=7800.
E=2.e11
c=np.sqrt(E/rho)

def computeLpNorm(Unum,Uexact,dx,p):
    return (((dx*((np.abs(Unum-Uexact))**p)).sum())**(1.0/p))

def computeRelativeError(Unum,Uexact,dx,p):
    return (computeLpNorm(Unum,Uexact,dx,p)/computeLpNorm(np.zeros(len(Unum)),Uexact,dx,p))

def export2DTeXFile(fileName,title,xField,xlabel,ylabel,fields,legend):
    TeXFile=open(fileName,"w")
    n_fields = np.shape(fields)[0]
    # Define Paul Tol's colors (purple to red)
    color=['Blue','Purple','Red','Orange','black','black','black']
    col=['120,28,129','63,96,174','83,158,182','109,179,136','202,184,67','231,133,50','217,33,32']
    marker=['*','triangle*','square','+','none','none','none']
    size=['very thick','very thick','very thick','very thick','thin','thin',]
    # for i in range(len(col)):
    #     TeXFile.write(r'\definecolor{'+color[i]+'}{RGB}{'+col[i]+'}')
    #     TeXFile.write('\n')
    TeXFile.write(r'\begin{tikzpicture}[scale=0.7]');TeXFile.write('\n')
    TeXFile.write(r'\begin{loglogaxis}[xlabel='+str(xlabel)+',ylabel='+str(ylabel)+',ymajorgrids=true,xmajorgrids=true,legend pos=south east,title={'+str(title)+'}]');TeXFile.write('\n')
    for i in range(np.shape(fields)[0]):
        TeXFile.write(r'\addplot['+str(color[i])+','+str(size[i])+',mark='+str(marker[i])+'] coordinates {')
        for j in range(len(fields[i,:])):
            TeXFile.write('('+str(xField[j])+','+str(fields[i,j])+') ')
        TeXFile.write('};\n')
        TeXFile.write(r'\addlegendentry{'+legend[i]+'}\n')
    
    TeXFile.write(r'\draw (axis cs:0.2,0.08) -- (axis cs:0.2/1.4,0.08/1.4);')
    TeXFile.write('\n')    
    TeXFile.write(r'\draw (axis cs:0.2,0.08) -- (axis cs:0.2,0.08/1.4) node [midway,right] {\scriptsize 1};')
    TeXFile.write('\n')    
    TeXFile.write(r'\draw (axis cs:0.2,0.08/1.4) -- (axis cs:0.2/1.4,0.08/1.4) node [midway,below] {\scriptsize 1};')
    TeXFile.write('\n')    
    TeXFile.write(r'\end{loglogaxis}')
    TeXFile.write('\n')
    TeXFile.write('\end{tikzpicture}')
    TeXFile.write('\n')
    TeXFile.close()

def func(x,C,n):
    return (C*(x**n))

titles=['a','b','c','d']
###Opening the files and computation of the solution by each method
###Parameters####
MP = [20,40,60,80,100,120,140,160,200,300]
ppc=[2]

# MP = [10,20,40,60,80,100]
# ppc=[2]
MP = [4,8,16,32,64,128]
MP = [4,8,16,32,64]
ppc=[2,3,4,8]
#ppc=[2,6,10,20]
ppc=[2]
#MP = [120,160,200,240,300]
loading=-75.
dx = np.zeros((len(MP),len(ppc)))
DX = np.zeros(len(MP))

ErrS_dgmpm = np.zeros((len(MP),len(ppc)))
ErrS_dgmpm_cons = np.zeros((len(MP),len(ppc)))

ErrV_dgmpm = np.zeros((len(MP),len(ppc)))
ErrV_dgmpm_cons = np.zeros((len(MP),len(ppc)))

Error_order1=2
Error_order2=2
period=0.4
DX=[]
compute_CFL=True
for p in range(len(ppc)):
    print "==============================",ppc[p],"PPC =============================="
    for i,Ne in enumerate(MP):
        parameters = {"Mp":Ne,"ppc":ppc[p],"loading":loading,"t_order":1,"period":period,"compute_CFL":compute_CFL,"lumping":True}
        ##DGMPM: Discontinuous Galerkin Material Point Method
        DGMPM = dict(parameters)
        execfile('dgmpm.py', DGMPM)

        parameters = {"Mp":Ne,"ppc":ppc[p],"loading":loading,"t_order":1,"period":period,"compute_CFL":compute_CFL,"lumping":False}
        ##DGMPM: Discontinuous Galerkin Material Point Method
        DGMPM2 = dict(parameters)
        execfile('dgmpm.py', DGMPM2)

        ########## Compute Errors #################
        dx[i,p]=DGMPM["dx"]
        n1= DGMPM["Increments"]
        n2= USL["Increments"]
        n1RK2= DGMPMRK2["Increments"]
        nFEM= FEM["Increments"]

        #Sigma
        ErrS_dgmpm[i,p] = computeRelativeError(DGMPM["Stress"][:,n1],DGMPM["Sth"][:,n1],dx[i,p],Error_order1)
        ErrS_dgmpm_cons[i,p] = computeRelativeError(DGMPM2["Stress"][:,n1],DGMPM["Sth"][:,n1],dx[i,p],Error_order1)
        
        #Velocity
        ErrV_dgmpm[i,p] = computeRelativeError(DGMPM["Velocity"][:,n1],-DGMPM["Sth"][:,n1]/(rho*c),dx[i,p],Error_order1)
        ErrV_dgmpm2[i,p] = computeRelativeError(DGMPM2["Velocity"][:,n1],-DGMPM2["Sth"][:,n1]/(rho*c),dx[i,p],Error_order1)
        
        
        #########################  Comparison  ######################################
        
        ##Stress
        l1,l2,l3=plt.plot(DGMPM["Pos"][:,n1],DGMPM["Velocity"][:,n1],'r+', \
                                DGMPM["Pos"][:,n1],-DGMPM["Sth"][:,n1]/(rho*c),'k', \
                                DGMPM2["Pos"][:,n1],DGMPM2["Velocity"][:,n1],'bo')
        
        plt.grid()
        plt.xlabel('x (m)')
        plt.ylabel(r'$\sigma$')
        plt.show()
        #############################################################################

    ## Fit curves
    errorList=[]
    errorList.append(ErrS_dgmpm[:,p])
    errorList.append(ErrS_dgmpm_cons[:,p])
    
    popt = np.zeros((len(errorList),2))
    for j,convCurve in enumerate(errorList):
        popt[j,:], pcov = curve_fit(func,dx[:,p],convCurve)

    Fitdgmpm=func(dx[:,p],popt[0,0],popt[0,1])
    Fitdgmpm_cons=func(dx[:,p],popt[1,0],popt[1,1])
    print "+++ Order of accuracy in stress :"
    print "dgmpm (lumped): ", popt[0,1]
    print "dgmpm (consistent): ", popt[1,1]
    
    
    soustitre='('+str(titles[p])+') '+str(ppc[p])+' particles per cell'
    title='convS_'+str(ppc[p])+'ppc.tex'

    ## Export Stress curves
    export2DTeXFile(title,soustitre,dx[:,p],r'$\Delta X (m)$',r'$\epsilon_\sigma$',np.array([ErrS_dgmpm[:,p],ErrS_dgmpm_cons[:,p]]),['dgmpm (lumped)','dgmpm (consistent)'])
    
    
    # plt.loglog(dx[:,p],ErrS_dgmpm[:,p],'r')
    # plt.loglog(dx[:,p],FitdgmpmEuler,'k')
    # plt.grid()
    # plt.show()
    
    errorList=[]
    errorList.append(ErrV_dgmpm[:,p])
    errorList.append(ErrV_dgmpm_cons[:,p])
    
    
    popt = np.zeros((len(errorList),2))
    for j,convCurve in enumerate(errorList):
        popt[j,:], pcov = curve_fit(func,dx[:,p],convCurve)
    Fitdgmpm=func(dx[:,p],popt[0,0],popt[0,1])
    Fitdgmpm_cons=func(dx[:,p],popt[1,0],popt[1,1])
    print "+++ Order of accuracy in velocity :"
    print "dgmpm (lumped): ", popt[0,1]
    print "dgmpm (consistent): ", popt[1,1]
    
    # if p==0:
    #     popt, pcov = curve_fit(func,DX,ErrV_fem)
    #     print "fem : ", popt[1]
    ## Export Velocity curves 
    #export2DTeXFile('convV_'+str(ppc[p])+'ppc.tex',soustitre,dx[:,p],r'$\Delta X$',r'$\epsilon_v$',np.array([ErrV_dgmpm[:,p],ErrV_dgmpmRK2[:,p],ErrV_usl[:,p],FitdgmpmEuler,FitdgmpmRK2,Fitmpm]),['dgmpm (Euler)','usl','dgmpm (RK2)'])
    title='convV_'+str(ppc[p])+'ppc.tex'
    export2DTeXFile(title,soustitre,dx[:,p],r'$\Delta X (m)$',r'$\epsilon_v$',np.array([ErrV_dgmpm[:,p],ErrV_dgmpm_cons[:,p]]),['dgmpm (lumped)','dgmpm (consistent)'])
        
###########Plot convergence curves####################
#####Assess the constant and slopes of the
#####convergence curves by fitting a power law with curve_fit


errorList=[]
for p in range(len(ppc)):
    errorList.append(ErrS_dgmpm[:,p])
    errorList.append(ErrS_dgmpm_cons[:,p])
    
popt = np.zeros((len(errorList),2))
for i,convCurve in enumerate(errorList):
    popt[i,:], pcov = curve_fit(func,dx[:,i/3],convCurve)


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
#Velocity
color_dgmpm=['r-+','y-+','c-+','g-+','b-+','m-+']
color_mpm=['r--+','y--+','c--+','g--+','b--+','m--+']
color_dgmpmRK2=['m-+','b-+','g-+','c-+','y-+','r-+']
color_dgmpm2=['m--','b--','g--','c--','y--','r--']
for p in range(len(ppc)):
    plt.loglog(dx[:,p],ErrS_dgmpm[:,p],color_dgmpm[p],lw=2.5,label='dgmpm lumped')
    
    plt.loglog(dx[:,p],ErrS_dgmpm_cons[:,p],color_dgmpmRK2[p],lw=2.5,label='dgmpm consistent')
    
    
    
plt.grid(True,which="both")
plt.xlabel('grid size')
plt.legend(numpoints=1,loc='best')
plt.ylabel(r'$\frac{||\sigma- \sigma_{exact}||_2}{||\sigma_{exact}||_2}$')
plt.show()

#####Plot of convergence curves#########
#Velocity
for p in range(len(ppc)):
    plt.loglog(dx[:,p],ErrV_dgmpm[:,p],color_dgmpm[p],lw=2.5,label='dgmpm lumped')
    
    plt.loglog(dx[:,p],ErrV_dgmpm_cons[:,p],color_dgmpmRK2[p],lw=2.5,label='dgmpm consistent')
    
    plt.loglog(DX,ErrV_fem,'y',lw=2.5,label='FEM') 
    
plt.grid(True,which="both")
plt.xlabel('grid size')
plt.legend(numpoints=1,loc='best')
plt.ylabel(r'$\frac{||\sigma- \sigma_{exact}||_2}{||\sigma_{exact}||_2}$')
plt.show()


