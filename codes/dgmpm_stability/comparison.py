#!/usr/bin/python

import numpy as np
from scipy import optimize
from sympy import *
import matplotlib.pyplot as plt
import pdb

def residualRK2(point,S,Sp):
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
    for p in range(Nmp):
        ## First order contributions
        D_mu = S1[p]*S1[point]/Sum1 + S2[p]*S2[point]/Sum2 + CFL*( S2[point]/Sum2 - S1[point]/Sum1 -Nmp*S2[p]*S2[point]/(Sum2**2) )
        ## Second order contributions
        D_mu += 0.5*Nmp*(CFL**2)*((S2[p]/Sum2)*(S1[point]/Sum1-S2[point]/Sum2) + (S2[point]/Sum2)*(Nmp*S2[p]/Sum2-1.)/Sum2)
        Res = Res +np.abs(D_mu)
    # Sum over material points in previous cell
    for p in range(Nmpp):
        ## First order contributions
        D_mu = CFL*Nmp*Sp2[p]*S1[point]/(Sum1*Sump2)
        ## Second order contributions
        D_mu +=0.5*Nmp*(CFL**2)*( S1[point]/(Sum1*Sump2)*(1.-Nmpp*Sp2[p]/Sump2) -(Sp2[p]/Sump2)*(S1[point]/Sum1-S2[point]/Sum2) )
        Res=Res + np.abs(D_mu)    
    Residual = lambdify((CFL),Res-1.)
    return Residual

def residualEuler(point,S,Sp):
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
    for p in range(Nmp):
        D_ma = S1[point]*S1[p]/Sum1 + S2[point]*S2[p]/Sum2 + CFL*( S2[point]/Sum2 - S1[point]/Sum1 -Nmp*S2[point]*S2[p]/(Sum2**2) )
        Res = Res +np.abs(D_ma)
    for p in range(Nmpp):
        D_ma = CFL*Nmp*S1[point]*Sp2[p]/(Sum1*Sump2)
        Res=Res + np.abs(D_ma)
    Residual = lambdify((CFL),Res-1.)
    return Residual

# Symbolic function to evaluate shape functions
shape_functions=lambda x: np.matrix([(1-x)/DX,x/DX])

xn = np.array([0.,1.])
DX = 1.

## required for plotting residual
CFL=np.linspace(0.,1.,100.)

shift=0.1

# 1PPC
print "**************************************************************"
print "******************  1PPC discretization **********************"
print "**************************************************************"
print "   "
shapes=shape_functions(0.25)
eulerSolution=optimize.root(residualEuler(0,shapes,shapes),1.,method='hybr',options={'xtol':1.e-4}).x
rk2Solution=optimize.root(residualRK2(0,shapes,shapes),1.,method='hybr',options={'xtol':1.e-4}).x
print "Euler solution, CFL= ",eulerSolution
print "RK2 solution, CFL= ",rk2Solution

# 2PPC
print "**************************************************************"
print "******************  2PPC discretization **********************"
print "**************************************************************"
print "   "
shift=0.25
shapes=shape_functions(np.array([0.25,0.75]))
## Gauss-Legendre integration
shapes=shape_functions(0.5*np.array([1.-1./np.sqrt(3.),1.+1./np.sqrt(3.)]))

eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[0]):
    eulerSolution.append(optimize.root(residualEuler(i,shapes,shapes),1.,method='hybr',options={'xtol':1.e-4}).x[0])
    rk2Solution.append(optimize.root(residualRK2(i,shapes,shapes),1.,method='hybr',options={'xtol':1.e-4}).x[0])
print "Euler solution, CFL= ",min(eulerSolution)
print "RK2 solution, CFL= ",min(rk2Solution)
pdb.set_trace()
print "   "
print "Shifted"
shapes=shape_functions(np.array([0.25+shift,0.75+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[0]):
    eulerSolution.append(optimize.root(residualEuler(i,shapes,shapes),1.,method='hybr',options={'xtol':1.e-4}).x[0])
    rk2Solution.append(optimize.root(residualRK2(i,shapes,shapes),1.,method='hybr',options={'xtol':1.e-4}).x[0])
print "Euler solution, CFL= ",min(eulerSolution)
print "RK2 solution, CFL= ",min(rk2Solution)

# 3PPC
print "**************************************************************"
print "******************  3PPC discretization **********************"
print "**************************************************************"
print "   "
shapes=shape_functions(np.array([0.25,0.5,0.75]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[0]):
    eulerSolution.append(optimize.root(residualEuler(i,shapes,shapes),1.,method='hybr',options={'xtol':1.e-4}).x[0])
    rk2Solution.append(optimize.root(residualRK2(i,shapes,shapes),1.,method='hybr',options={'xtol':1.e-4}).x[0])
print "Euler solution, CFL= ",min(eulerSolution)
print "RK2 solution, CFL= ",min(rk2Solution)

print "   "
print "Shifted"
shapes=shape_functions(np.array([0.25+shift,0.5+shift,0.75+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[0]):
    eulerSolution.append(optimize.root(residualEuler(i,shapes,shapes),1.,method='hybr',options={'xtol':1.e-4}).x[0])
    rk2Solution.append(optimize.root(residualRK2(i,shapes,shapes),1.,method='hybr',options={'xtol':1.e-4}).x[0])
print "Euler solution, CFL= ",min(eulerSolution)
print "RK2 solution, CFL= ",min(rk2Solution)

# 4PPC
print "**************************************************************"
print "******************  4PPC discretization **********************"
print "**************************************************************"
print "   "
shapes=shape_functions(np.array([0.2,0.4,0.6,0.8]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[0]):
    eulerSolution.append(optimize.root(residualEuler(i,shapes,shapes),1.,method='hybr',options={'xtol':1.e-4}).x[0])
    rk2Solution.append(optimize.root(residualRK2(i,shapes,shapes),1.,method='hybr',options={'xtol':1.e-4}).x[0])
print "Euler solution, CFL= ",min(eulerSolution)
print "RK2 solution, CFL= ",min(rk2Solution)

print "   "
print "Shifted"
shapes=shape_functions(np.array([0.2+shift,0.4+shift,0.6+shift,0.8+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[0]):
    eulerSolution.append(optimize.root(residualEuler(i,shapes,shapes),1.,method='hybr',options={'xtol':1.e-4}).x[0])
    rk2Solution.append(optimize.root(residualRK2(i,shapes,shapes),1.,method='hybr',options={'xtol':1.e-4}).x[0])
print "Euler solution, CFL= ",min(eulerSolution)
print "RK2 solution, CFL= ",min(rk2Solution)
