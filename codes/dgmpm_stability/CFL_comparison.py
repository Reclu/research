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
    for k in range(Nmp):
        ## First order contributions
        D_mu = S1[k]*S1[point]/Sum1 + S2[k]*S2[point]/Sum2 + CFL*( S2[point]/Sum2 - S1[point]/Sum1 -Nmp*S2[k]*S2[point]/(Sum2**2) )
        ## Second order contributions
        D_mu += 0.5*Nmp*(CFL**2)*((S2[k]/Sum2)*(S1[point]/Sum1-S2[point]/Sum2) + (S2[point]/Sum2)*(Nmp*S2[k]/Sum2-1.)/Sum2)
        Res = Res +np.abs(D_mu)
    # Sum over material points in previous cell
    for k in range(Nmpp):
        ## First order contributions
        D_mu = CFL*Nmp*Sp2[k]*S1[point]/(Sum1*Sump2)
        ## Second order contributions
        D_mu +=0.5*Nmp*(CFL**2)*( S1[point]/(Sum1*Sump2)*(1.-Nmp*Sp2[k]/Sump2) -(Sp2[k]/Sump2)*(S1[point]/Sum1-S2[point]/Sum2) )
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
    for k in range(Nmp):
        D_ma = S1[point]*S1[k]/Sum1 + S2[point]*S2[k]/Sum2 + CFL*( S2[point]/Sum2 - S1[point]/Sum1 -Nmp*S2[point]*S2[k]/(Sum2**2) )
        Res = Res +np.abs(D_ma)
    for k in range(Nmpp):
        D_ma = CFL*Nmp*S1[point]*Sp2[k]/(Sum1*Sump2)
        Res=Res + np.abs(D_ma)
    Residual = lambdify((CFL),Res-1.)
    return Residual

# def gridSearch(function,tol=1.e-7):
#     samples=100000
#     # Find the bigest root of the residual by grid search algorithm
#     CFL=np.linspace(0.,1.,samples)
#     for i in CFL:
#         if i==CFL[samples-1]: return i
#         a0=function(i)
#         if a0<tol:
#             continue
#         else:
#             return i
def gridSearch(function,tol=1.e-7):
    samples=100000
    # Find the bigest root of the residual by grid search algorithm
    CFL=np.linspace(0.,1.,samples)
    for i in range(samples):
        value=CFL[-1-i]
        a0=function(value)
        if a0<tol:
            return value
        else:
            continue
    return 0.

# Symbolic function to evaluate shape functions
shape_functions=lambda x: np.matrix([(1-x)/DX,x/DX])

xn = np.array([0.,1.])
DX = 1.

## required for plotting residual
CFL=np.linspace(0.,1.,100.)

shift=0.25

# # 1PPC
# print "**************************************************************"
# print "******************  1PPC discretization **********************"
# print "**************************************************************"
# print "   "
# shapes=shape_functions(0.00000000000001)
# eulerSolution=optimize.root(residualEuler(0,shapes,shapes),1.,method='hybr',options={'xtol':1.e-4}).x
# rk2Solution=optimize.root(residualRK2(0,shapes,shapes),1.,method='hybr',options={'xtol':1.e-4}).x
# print "Euler solution, CFL= ",eulerSolution
# print "RK2 solution, CFL= ",rk2Solution

# 2PPC
print "**************************************************************"
print "******************  2PPC discretization **********************"
print "**************************************************************"
print "   "
shapes=shape_functions(np.array([0.25,0.75]))
#shapes=shape_functions(np.array([(1.-1./np.sqrt(3.))/2.,(1.+1./np.sqrt(3.))/2.]))
eulerSolution=[]
rk2Solution=[]
    
for i in range(np.shape(shapes)[0]):
    # eulerSolution.append(optimize.root(residualEuler(i,shapes,shapes),1.,method='hybr',options={'xtol':1.e-12}).x[0])
    # rk2Solution.append(optimize.root(residualRK2(i,shapes,shapes),1.,method='hybr',options={'xtol':1.e-12}).x[0])
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ", (eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)

print "   "
print "Shifted"
shift=DX/10.
shapes=shape_functions(np.array([0.25+shift,0.75+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[1]):
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ", (eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)

print "   "
print "Shifted left on nodes"
shift=-0.25
shapes=shape_functions(np.array([0.25+shift,0.75+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[1]):
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ", (eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)

print "   "
print "Shifted right on nodes"
shift=+0.25
shapes=shape_functions(np.array([0.25+shift,0.75+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[1]):
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ", (eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)

print "   "
print "Shifted symmetrically"
shift=+DX/10.
shapes=shape_functions(np.array([0.25-shift,0.75+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[1]):
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ", (eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)

print "   "
print "Shifted symmetrically to nodes"
shift=0.25
shapes=shape_functions(np.array([0.25-shift,0.75+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[1]):
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ", (eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)

# 3PPC
print "**************************************************************"
print "******************  3PPC discretization **********************"
print "**************************************************************"
print "   "
shift=0.
shapes=shape_functions(np.array([DX/6.,0.5,1-DX/6.]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[1]):
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ", (eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)

print "   "
print "Shifted"
shift=0.1
shapes=shape_functions(np.array([0.5*DX/3.+shift,0.5+shift,1.-0.5*DX/3.+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[1]):
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ",(eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)

print "   "
print "Shifted left on nodes"
shift=-DX/6.
shapes=shape_functions(np.array([0.5*DX/3.+shift,0.5+shift,1.-0.5*DX/3.+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[1]):
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ",(eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)

print "   "
print "Shifted right on nodes"
shift=DX/6.
shapes=shape_functions(np.array([0.5*DX/3.+shift,0.5+shift,1.-0.5*DX/3.+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[1]):
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ",(eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)

print "   "
print "Shifted symetrically"
shift=0.1
shapes=shape_functions(np.array([0.5*DX/3.-shift,0.5,1.-0.5*DX/3.+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[1]):
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ", (eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)

print "   "
print "Shifted symmetrically to nodes"
shift=DX/6.
shapes=shape_functions(np.array([0.5*DX/3.-shift,0.5,1-0.5*DX/3.+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[1]):
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ", (eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)

# 4PPC
print "**************************************************************"
print "******************  4PPC discretization **********************"
print "**************************************************************"
print "   "

shapes=shape_functions(np.array([0.5*DX/4.,3*DX/8.,5*DX/8.,7*DX/8.]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[1]):
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ", (eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)

print "   "
print "Shifted"
shift=0.1
shapes=shape_functions(np.array([0.5*DX/4.+shift,3*DX/8.+shift,5*DX/8.+shift,7*DX/8.+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[1]):
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ", (eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)

print "   "
print "Shifted left on nodes"
shift=-DX/8.
shapes=shape_functions(np.array([0.5*DX/4.+shift,3*DX/8.+shift,5*DX/8.+shift,7*DX/8.+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[1]):
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ", (eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)

print "   "
print "Shifted right on nodes"
shift=DX/8.
shapes=shape_functions(np.array([0.5*DX/4.+shift,3*DX/8.+shift,5*DX/8.+shift,7*DX/8.+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[1]):
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ", (eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)


print "   "
shift=0.1
print "Shifted symetrically"
shapes=shape_functions(np.array([0.5*DX/4.-shift,3*DX/8.-shift,5*DX/8.+shift,7*DX/8.+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[1]):
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ", (eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)

print "   "
print "Shifted symetrically to nodes"
shift=DX/8.
shapes=shape_functions(np.array([0.5*DX/4.-shift,3*DX/8.-shift,5*DX/8.+shift,7*DX/8.+shift]))
eulerSolution=[]
rk2Solution=[]
for i in range(np.shape(shapes)[1]):
    eulerSolution.append(gridSearch(residualEuler(i,shapes,shapes)))
    rk2Solution.append(gridSearch(residualRK2(i,shapes,shapes)))
print "Euler solution, CFL= ", (eulerSolution)
print "RK2 solution, CFL= ", (rk2Solution)
