#!/usr/bin/python

import numpy as np
from scipy import optimize
from sympy import *
import matplotlib.pyplot as plt
import pdb

def symbolResidual(point,S,Sp):
    CFL = symbols('CFL')
    Res=0.
    if S.shape[0]==1:
        S1=[S[0,0]]
        S2=[S[0,1]]
        Sum1=np.sum(S1) ; Sum2=np.sum(S2)
        Nmp=1
    else:
        S1=np.asarray(S[0,:])[0]
        S2=np.asarray(S[1,:])[0]
        Sum1=np.sum(S1) ; Sum2=np.sum(S2)
        Nmp=len(S1) 
    if Sp.shape[0]==1:
        Sp1=[Sp[0,0]]
        Sp2=[Sp[0,0]]
        Sump1=np.sum(Sp1) ; Sump2=np.sum(Sp2)
        Nmpp=1
    else:
        Sp1=np.asarray(Sp[0,:])[0]
        Sp2=np.asarray(Sp[1,:])[0]
        Sump1=np.sum(Sp1) ; Sump2=np.sum(Sp2)
        Nmpp=len(Sp1) 
    # Sum over material points in curent cell
    #pdb.set_trace()
    for p in range(Nmp):
        ## First order contributions
        D_mu = S1[p]*S1[point]/Sum1 + S2[p]*S2[point]/Sum2 + CFL*( S2[point]/Sum2 - S1[point]/Sum1 -Nmp*S2[p]*S2[point]/(Sum2**2) )
        ## Second order contributions
        D_mu += 0.5*Nmp*(CFL**2)*((S2[p]/Sum2)*(S1[point]/Sum1-S2[point]/Sum2) + (S2[point]/Sum2)*(Nmp*S2[p]/Sum2-1.)/Sum2)
        Res = Res +np.abs(D_mu)
    # Sum over material points in previous cell
    #pdb.set_trace()
    for p in range(Nmpp):
        ## First order contributions
        D_mu = CFL*Nmp*Sp2[p]*S1[point]/(Sum1*Sump2)
        ## Second order contributions
        D_mu +=0.5*Nmp*(CFL**2)*( S1[point]/(Sum1*Sump2)*(1.-Nmpp*Sp2[p]/Sump2) -(Sp2[p]/Sump2)*(S1[point]/Sum1-S2[point]/Sum2) )
        Res=Res + np.abs(D_mu)    
    Residual = lambdify((CFL),Res-1.)
    #Residual = lambdify((CFL),np.abs(Res)-1.)
    return Residual

def computeResidual(CFL,point,S,Sp):
    Res=0.
    #pdb.set_trace()
    if S.shape[0]==1:
        S1=[S[0,0]]
        S2=[S[0,1]]
        Sum1=np.sum(S1) ; Sum2=np.sum(S2)
        Nmp=1
    else:
        S1=np.asarray(S[0,:])[0]
        S2=np.asarray(S[1,:])[0]
        Sum1=np.sum(S1) ; Sum2=np.sum(S2)
        Nmp=len(S1) 
    if Sp.shape[0]==1:
        Sp1=[Sp[0,0]]
        Sp2=[Sp[0,0]]
        Sump1=np.sum(Sp1) ; Sump2=np.sum(Sp2)
        Nmpp=1
    else:
        Sp1=np.asarray(Sp[0,:])[0]
        Sp2=np.asarray(Sp[1,:])[0]
        Sump1=np.sum(Sp1) ; Sump2=np.sum(Sp2)
        Nmpp=len(Sp1)
    # Sum over material points in curent cell
    for p in range(Nmp):
        ## First order contributions
        D_mu = S1[point]*S1[p]/Sum1 + S2[point]*S2[p]/Sum2 + CFL*( S2[point]/Sum2 - S1[point]/Sum1 -Nmp*S2[point]*S2[p]/(Sum2**2) )
        ## Second order contributions
        D_mu += 0.5*Nmp*(CFL**2)*((S2[p]/Sum2)*(S1[point]/Sum1-S2[point]/Sum2) + (S2[point]/(Sum2**2))*(Nmp*S2[p]/Sum2-1.) )
        Res = Res +np.abs(D_mu)
    # Sum over material points in previous cell
    #pdb.set_trace()
    for p in range(Nmpp):
        ## First order contributions
        D_mu = CFL*Nmp*Sp2[p]*S1[point]/(Sum1*Sump2)
        ## Second order contributions
        D_mu +=0.5*Nmp*(CFL**2)*( S1[point]/(Sum1*Sump2)*(1.-Nmpp*Sp2[p]/Sump2) -(Sp2[p]/Sump2)*(S1[point]/Sum1-S2[point]/Sum2) )
        Res=Res + np.abs(D_mu)  
    return Res

# Symbolic function to evaluate shape functions
shape_functions=lambda x: np.matrix([(1-x)/DX,x/DX])

xn = np.array([0.,1.])
DX = 1.


CFL=np.linspace(0.,1.,100.)

############### 1PPC
print "**************************************************************"
print "******************  1PPC discretization **********************"
print "**************************************************************"
solution=optimize.newton(symbolResidual(0,shape_functions(0.5),shape_functions(0.5)),1.)
print "Solution is: ",solution
############### 2PPC
print "**************************************************************"
print "******************  2PPC discretization **********************"
print "**************************************************************"
solution=optimize.newton(symbolResidual(0,shape_functions(np.array([0.25,0.75])),shape_functions(np.array([0.25,0.75]))),1.)
print "Solution is: ",solution
############### 3PPC
print "**************************************************************"
print "******************  3PPC discretization **********************"
print "**************************************************************"
solution=optimize.newton(symbolResidual(0,shape_functions(np.array([0.25,0.5,0.75])),shape_functions(np.array([0.25,0.5,0.75]))),1.)
print "Solution is: ",solution
############### 4PPC
print "**************************************************************"
print "******************  4PPC discretization **********************"
print "**************************************************************"
solution=optimize.newton(symbolResidual(0,shape_functions(np.array([0.2,0.4,0.6,0.8])),shape_functions(np.array([0.2,0.4,0.6,0.8]))),1.)
print "Solution is: ",solution


#f_2mp_centered1=computeResidual(CFL,0,shape_functions(np.array([0.])),shape_functions(np.array([0.,0.75])))

#f_1mp=computeResidual(CFL,0,shape_functions(0.5),shape_functions(0.5))
#f_1mp_shifted=computeResidual(CFL,0,shape_functions(0.75),shape_functions(0.75))

#f_2mp_centered1=computeResidual(CFL,0,shape_functions(np.array([0.25,0.75])),shape_functions(np.array([0.25,0.75])))
f_2mp_centered1=computeResidual(CFL,0,shape_functions(0.5),shape_functions(0.5))
f_2mp_centered2=computeResidual(CFL,1,shape_functions(np.array([0.25,0.75])),shape_functions(np.array([0.25,0.75])))
f_2mp_shifted1=computeResidual(CFL,0,shape_functions(np.array([0.25+0.2,0.75+0.2])),shape_functions(np.array([0.25+0.2,0.75+0.2])))
f_2mp_shifted2=computeResidual(CFL,1,shape_functions(np.array([0.25+0.2,0.75+0.2])),shape_functions(np.array([0.25+0.2,0.75+0.2])))




plt.plot(CFL,f_2mp_centered1,'gs',lw=4.5,ms=14.,label='Symmetrical configuration (point 1)')
plt.plot(CFL,f_2mp_centered2,'b*',lw=2.5,ms=14.,label='Symmetrical configuration (point 2)')
plt.plot(CFL,f_2mp_shifted1,'r--',lw=2.5,ms=14.,label='Shifted configuration (point 1)')
plt.plot(CFL,f_2mp_shifted2,'k',lw=2.5,ms=10.,label='Shifted configuration (point 2)')
plt.xlabel('CFL number',size=34.)
plt.ylabel('Amplification factor',size=34.)
plt.xticks(size=34.)
plt.yticks(size=34.)
plt.legend(numpoints=1,loc='best',prop={'size':34})
plt.grid()
plt.show()

"""
shift=0.1
f_3mp_centered1=computeResidual(CFL,0,shape_functions(np.array([0.5-1./3.,0.5,0.5+1./3.])),shape_functions(np.array([0.5-1./3.,0.5,0.5+1./3.])))
f_3mp_centered2=computeResidual(CFL,1,shape_functions(np.array([0.5-1./3.,0.5,0.5+1./3.])),shape_functions(np.array([0.5-1./3.,0.5,0.5+1./3.])))
f_3mp_centered3=computeResidual(CFL,2,shape_functions(np.array([0.5-1./3.,0.5,0.5+1./3.])),shape_functions(np.array([0.5-1./3.,0.5,0.5+1./3.])))
f_3mp_shifted1=computeResidual(CFL,0,shape_functions(np.array([0.5-1./3.+shift,0.5+shift,0.5+1./3.+shift])),shape_functions(np.array([0.5-1./3.+shift,0.5+shift,0.5+1./3.+shift])))
f_3mp_shifted2=computeResidual(CFL,1,shape_functions(np.array([0.5-1./3.+shift,0.5+shift,0.5+1./3.+shift])),shape_functions(np.array([0.5-1./3.+shift,0.5+shift,0.5+1./3.+shift])))
f_3mp_shifted3=computeResidual(CFL,2,shape_functions(np.array([0.5-1./3.+shift,0.5+shift,0.5+1./3.+shift])),shape_functions(np.array([0.5-1./3.+shift,0.5+shift,0.5+1./3.+shift])))

solution=[]
solution.append(optimize.newton(symbolResidual(0,shape_functions(np.array([0.5-1./3.,0.5,0.5+1./3.])),shape_functions(np.array([0.5-1./3.,0.5,0.5+1./3.]))),1.))
solution.append(optimize.newton(symbolResidual(1,shape_functions(np.array([0.5-1./3.,0.5,0.5+1./3.])),shape_functions(np.array([0.5-1./3.,0.5,0.5+1./3.]))),1.))
solution.append(optimize.newton(symbolResidual(2,shape_functions(np.array([0.5-1./3.,0.5,0.5+1./3.])),shape_functions(np.array([0.5-1./3.,0.5,0.5+1./3.]))),1.))
solution_s=[]

solution_s.append(optimize.newton(symbolResidual(0,shape_functions(np.array([0.5-1./3.+shift,0.5+shift,0.5+1./3.+shift])),shape_functions(np.array([0.5-1./3.+shift,0.5+shift,0.5+1./3.+shift]))),1.))
solution_s.append(optimize.newton(symbolResidual(1,shape_functions(np.array([0.5-1./3.+shift,0.5+shift,0.5+1./3.+shift])),shape_functions(np.array([0.5-1./3.+shift,0.5+shift,0.5+1./3.+shift]))),1.))
solution_s.append(optimize.newton(symbolResidual(2,shape_functions(np.array([0.5-1./3.+shift,0.5+shift,0.5+1./3.+shift])),shape_functions(np.array([0.5-1./3.+shift,0.5+shift,0.5+1./3.+shift]))),1.))

plt.plot(np.array([0.5-1./3.,0.5,0.5+1./3.]),np.zeros(3),'bo',ms=15)
plt.plot(np.array([0.5-1./3.+shift,0.5+shift,0.5+1./3.+shift]),np.zeros(3),'ro',ms=15)
plt.grid()
plt.show()
print "Symmetrical solution is: ",np.min(solution)," ; Shiffted solution is : ",np.min(solution_s)

plt.plot(CFL,f_3mp_centered1,'gs',lw=4.5,ms=14.,label='Symmetrical configuration (point 1)')
plt.plot(CFL,f_3mp_centered2,'b*',lw=2.5,ms=14.,label='Symmetrical configuration (point 2)')
plt.plot(CFL,f_3mp_centered3,'yo',lw=4.5,ms=14.,label='Symmetrical configuration (point 3)')
plt.plot(CFL,f_3mp_shifted1,'r--',lw=2.5,ms=14.,label='Shifted configuration (point 1)')
plt.plot(CFL,f_3mp_shifted2,'k',lw=2.5,ms=10.,label='Shifted configuration (point 2)')
plt.plot(CFL,f_3mp_shifted3,'m>',lw=4.5,ms=14.,label='Shifted configuration (point 3)')
plt.xlabel('CFL number',size=34.)
plt.ylabel('Amplification factor',size=34.)
plt.xticks(size=34.)
plt.yticks(size=34.)
plt.legend(numpoints=1,loc='best',prop={'size':34})
plt.grid()
plt.show()
"""
"""
f_2mp_centered1=computeResidual(CFL,0,shape_functions(np.array([0.25,0.75])),shape_functions(np.array([0.25,0.75])))
f_2mp_centered1_1p=computeResidual(CFL,0,shape_functions(np.array([0.25,0.75])),shape_functions(np.array([0.5])))
f_2mp_centered1_3p=computeResidual(CFL,0,shape_functions(np.array([0.25,0.75])),shape_functions(np.array([0.25,0.5,0.75])))
f_2mp_centered2_3p=computeResidual(CFL,1,shape_functions(np.array([0.25,0.75])),shape_functions(np.array([0.25,0.5,0.75])))

plt.plot(CFL,f_2mp_centered1,'gs',lw=4.5,ms=8.,label='Regular material point set')
plt.plot(CFL,f_2mp_centered1_1p,'b*',lw=2.5,ms=10.,label='One upwind material point')
plt.plot(CFL,f_2mp_centered1_3p,'r.-',lw=2.5,ms=14.,label='Three upwind material points')
#plt.plot(CFL,f_2mp_centered2_3p,'k',lw=2.5,ms=10.,label='Three upwind particles (2nd point condition)')
plt.xlabel('CFL number',size=28.)
plt.ylabel('Amplification factor',size=28.)
plt.xticks(size=28.)
plt.yticks(size=28.)
plt.legend(numpoints=1,loc='best',prop={'size':28})
plt.grid()
plt.show()
"""

"""
pos1 = (-np.sqrt(3.)/3.)/2. +0.5
pos2 = (np.sqrt(3.)/3.)/2. +0.5
f_exact_integ=computeResidual(CFL,0,shape_functions(np.array([pos1,pos2])),shape_functions(np.array([pos1,pos2])))
f_exact_integ2=computeResidual(CFL,1,shape_functions(np.array([pos1,pos2])),shape_functions(np.array([pos1,pos2])))

solution=[]
solution.append(optimize.newton(symbolResidual(0,shape_functions(np.array([pos1,0.,pos2])),shape_functions(np.array([pos1,0.,pos2]))),1.))
solution.append(optimize.newton(symbolResidual(1,shape_functions(np.array([pos1,0.,pos2])),shape_functions(np.array([pos1,0.,pos2]))),1.))

print "Symmetrical solution is: ",np.min(solution) 

plt.plot(CFL,f_exact_integ,'gs',lw=4.5,ms=8.,label='Exact integration (1st point condition)')
plt.plot(CFL,f_exact_integ2,'ro',lw=4.5,ms=8.,label='Exact integration (2nd point condition)')
plt.xlabel('CFL number',size=28.)
plt.ylabel('Amplification factor',size=28.)
plt.xticks(size=28.)
plt.yticks(size=28.)
plt.legend(numpoints=1,loc='best',prop={'size':28})
plt.grid()
plt.show()
"""

