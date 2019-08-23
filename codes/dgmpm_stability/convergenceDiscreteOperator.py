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
    return (((((np.abs(Unum-Uexact))**p)).sum())**(1.0/p))

def computeRelativeError(Unum,Uexact,dx,p):
    return (dx**(1.0/p)*computeLpNorm(Unum,Uexact,dx,p)/computeLpNorm(np.zeros(len(Unum)),Uexact,dx,p))

def export2DTeXFile(fileName,title,xField,xlabel,ylabel,fields,legend):
    TeXFile=open(fileName,"w")
    n_fields = np.shape(fields)[0]
    # Define Paul Tol's colors (purple to red)
    color=['Blue','Blue','Red','Red','black','black','black']
    col=['120,28,129','63,96,174','83,158,182','109,179,136','202,184,67','231,133,50','217,33,32']
    marker=['none','triangle*','none','square*','none','none','none']
    size=['very thick,dashed','very thick,only marks','very thick,dotted','very thick,only marks','thin','thin',]
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
    
    # TeXFile.write(r'\draw (axis cs:0.2,0.08) -- (axis cs:0.2/1.4,0.08/1.4);')
    # TeXFile.write('\n')    
    # TeXFile.write(r'\draw (axis cs:0.2,0.08) -- (axis cs:0.2,0.08/1.4) node [midway,right] {\scriptsize 1};')
    # TeXFile.write('\n')    
    # TeXFile.write(r'\draw (axis cs:0.2,0.08/1.4) -- (axis cs:0.2/1.4,0.08/1.4) node [midway,below] {\scriptsize 1};')
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
MP = [4,8,16,32,64]
MP = [8,16,32,64,128]
#MP = [32]

ppc=[2,4]

dx = np.zeros((len(MP),len(ppc)))
DX = np.zeros(len(MP))

Err_dgmpm = np.zeros((len(MP),len(ppc)))
Err_dgmpm_H = np.zeros((len(MP),len(ppc)))

Err_dgmpm_rk2 = np.zeros((len(MP),len(ppc)))
Err_dgmpm_rk2_H = np.zeros((len(MP),len(ppc)))


Error_order1=2
Error_order2=2
DX=[]
for p in range(len(ppc)):
    print "==============================",ppc[p],"PPC =============================="
    for i,Ne in enumerate(MP):
        parameters = {"Mp":Ne,"ppc":ppc[p],"t_order":1}
        ##DGMPM: Discontinuous Galerkin Material Point Method
        DGMPM = dict(parameters)
        execfile('linearAdvection_convergence.py', DGMPM)

        parameters = {"Mp":Ne,"ppc":ppc[p],"t_order":2}
        ##DGMPM: Discontinuous Galerkin Material Point Method
        DGMPM2 = dict(parameters)
        execfile('linearAdvection_convergence.py', DGMPM2)

        ########## Compute Errors #################
        dx[i,p]=DGMPM["dx"]
        n1= DGMPM["Increments"]
        n1RK2= DGMPM2["Increments"]
        
        #Sigma
        Err_dgmpm[i,p] = computeRelativeError(DGMPM["Stress"][:,n1],DGMPM["Exact"][:,n1],dx[i,p],Error_order1)
        Err_dgmpm_H[i,p] = computeRelativeError(DGMPM["Stressh"][:,n1],DGMPM["Exact"][:,n1],dx[i,p],Error_order1)

        Err_dgmpm_rk2[i,p] = computeRelativeError(DGMPM2["Stress"][:,n1RK2],DGMPM2["Exact"][:,n1RK2],dx[i,p],Error_order1)
        Err_dgmpm_rk2_H[i,p] = computeRelativeError(DGMPM2["Stressh"][:,n1RK2],DGMPM2["Exact"][:,n1RK2],dx[i,p],Error_order1)
        
        # plt.plot(DGMPM2["xp"][:,0],DGMPM2["Exact"][:,n1RK2],'k--')
        # plt.plot(DGMPM2["xp"][:,0],DGMPM2["Stress"][:,n1RK2],'r')
        # plt.plot(DGMPM2["xp"][:,0],DGMPM2["Stressh"][:,n1RK2],'yo')
        # plt.grid()
        # plt.show()
    ## Fit curves
    errorList=[]
    errorList.append(Err_dgmpm[:,p])
    errorList.append(Err_dgmpm_H[:,p])
    errorList.append(Err_dgmpm_rk2[:,p])
    errorList.append(Err_dgmpm_rk2_H[:,p])
    
    popt = np.zeros((len(errorList),2))
    for j,convCurve in enumerate(errorList):
        popt[j,:], pcov = curve_fit(func,dx[:,p],convCurve)

    Fitdgmpm=func(dx[:,p],popt[0,0],popt[0,1])
    Fitdgmpm_H=func(dx[:,p],popt[1,0],popt[1,1])
    Fitdgmpm_rk2=func(dx[:,p],popt[2,0],popt[2,1])
    Fitdgmpm_rk2_H=func(dx[:,p],popt[3,0],popt[3,1])
    print "+++ Order of accuracy in stress :"
    print "dgmpm scheme: ", popt[0,1]
    print "dgmpm discrete: ", popt[1,1]
    print "dgmpm scheme + rk2: ", popt[2,1]
    print "dgmpm discrete + rk2: ", popt[3,1]
    
    
    soustitre=''
    title='conv1D_'+str(ppc[p])+'ppc.tex'

    ## Export Stress curves
    export2DTeXFile(title,soustitre,dx[:,p],r'$\Delta X (m)$',r'$\epsilon_q$',np.array([Err_dgmpm[:,p],Err_dgmpm_H[:,p],Err_dgmpm_rk2[:,p],Err_dgmpm_rk2_H[:,p]]),['dgmpm euler ','dgmpm euler (discrete)','dgmpm rk2 ','dgmpm rk2 (discrete)'])
    
    
    plt.loglog(dx[:,p],Err_dgmpm[:,p],'r--')
    plt.loglog(dx[:,p],Err_dgmpm_H[:,p],'r*')
    plt.loglog(dx[:,p],Err_dgmpm_rk2[:,p],'y--')
    plt.loglog(dx[:,p],Err_dgmpm_rk2_H[:,p],'yo')
    # plt.loglog(dx[:,p],FitdgmpmEuler,'k')
    plt.grid()
    plt.show()
    
        
###########Plot convergence curves####################
#####Assess the constant and slopes of the
#####convergence curves by fitting a power law with curve_fit

"""
errorList=[]
for p in range(len(ppc)):
    errorList.append(Err_dgmpm[:,p])
    errorList.append(Err_dgmpm_H[:,p])
    
popt = np.zeros((len(errorList),2))
for i,convCurve in enumerate(errorList):
    popt[i,:], pcov = curve_fit(func,dx[:,i/3],convCurve)

#pdb.set_trace()
#####Plot of convergence curves#########
fig1 = plt.figure()
#Velocity
color_dgmpm=['r-+','y-+','c-+','g-+','b-+','m-+']
color_mpm=['r--+','y--+','c--+','g--+','b--+','m--+']
color_dgmpmRK2=['m-+','b-+','g-+','c-+','y-+','r-+']
color_dgmpm2=['m--','b--','g--','c--','y--','r--']
for p in range(len(ppc)):
    plt.loglog(dx[:,p],Err_dgmpm[:,p],color_dgmpm[p],lw=2.5,label='dgmpm scheme')
    
    plt.loglog(dx[:,p],Err_dgmpm_H[:,p],color_dgmpmRK2[p],lw=2.5,label='dgmpm consistent')
    
    
    
plt.grid(True,which="both")
plt.xlabel('grid size')
plt.legend(numpoints=1,loc='best')
plt.ylabel(r'$\frac{||\sigma- \sigma_{exact}||_2}{||\sigma_{exact}||_2}$')
plt.show()


"""
