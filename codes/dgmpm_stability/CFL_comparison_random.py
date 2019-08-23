#!/usr/bin/python

import numpy as np
from scipy import optimize
from sympy import *
import matplotlib.pyplot as plt
import pdb
import random
import os

def export2DTeXFile(fileName,bins,fields,*kwargs):
    TeXFile=open(fileName,"w")
    
    n_fields = np.shape(fields)[0]
    n_labels = np.shape(kwargs)[1]
    # Define Paul Tol's colors (purple to red)
    color=['Blue','Red','Green','Red','black','black','black']
    marker=['+','x','star','+','none','none','none']
    size=['very thick','very thick','very thick','very thick','thin','thin',]
    line=['solid','solid','dashed','dashed']
    TeXFile.write(r'\begin{tikzpicture}')
    TeXFile.write('\n')
    TeXFile.write(r'\begin{axis}[xlabel=CFL,x label style={at={(axis description cs:0.5,-0.25)}},ylabel=Density,legend pos = north west,ybar interval=0.7,xticklabel interval boundaries,x tick label style = {rotate=90,anchor=east}]')
    TeXFile.write('\n')
    TeXFile.write('%%%%%%%%%%% EULER SOLUTION')
    TeXFile.write('\n')
    #pdb.set_trace()
    for i in range(n_fields):
        TeXFile.write(r'\addplot['+str(color[i])+',fill='+str(color[i])+'] coordinates {')
        for j in range(len(fields[i,:])):
            TeXFile.write('('+str(bins[j])+','+str(fields[i,j])+') ')
        TeXFile.write('(1.,0.) ')
        TeXFile.write('};\n')
        TeXFile.write(r'\addlegendentry{'+str(kwargs[0][i])+'}; \n')
        if i==0:
            TeXFile.write('%%%%%%%%%%% RK2 SOLUTION')
            TeXFile.write('\n')

    TeXFile.write(r'\end{axis}')
    TeXFile.write('\n')
    TeXFile.write('\end{tikzpicture}')
    TeXFile.write('\n')
    TeXFile.write('%%% Local Variables:')
    TeXFile.write('\n')
    TeXFile.write('%%% mode: latex')
    TeXFile.write('\n')
    TeXFile.write('%%% TeX-master: "../manuscript"')
    TeXFile.write('\n')
    TeXFile.write('%%% End:')
    TeXFile.write('\n')
    TeXFile.close()


def residual(point,shapes,shapes_prev,t_order):
    CFL = symbols('CFL')
    Res=0.
    Nmp=len(shapes)
    Nmpp=len(shapes_prev)
    S1=shapes ; Sum1=np.sum(S1)
    S2=1.-shapes; Sum2=np.sum(S2)

    Sp1=shapes_prev ; Sump1=np.sum(Sp1)
    Sp2=1.-shapes_prev ; Sump2=np.sum(Sp2)
    
        
    # Sum over material points in curent cell
    for k in range(Nmp):
        ## First order contributions
        D_mu = S1[k]*S1[point]/Sum1 + S2[k]*S2[point]/Sum2 + CFL*( S2[point]/Sum2 - S1[point]/Sum1 -Nmp*S2[k]*S2[point]/(Sum2**2) )
        ## Second order contributions
        if t_order==2:
            D_mu += 0.5*Nmp*(CFL**2)*((S2[k]/Sum2)*(S1[point]/Sum1-S2[point]/Sum2) + (S2[point]/Sum2)*(Nmp*S2[k]/Sum2-1.)/Sum2)
        Res = Res +np.abs(D_mu)
    # Sum over material points in previous cell
    for k in range(Nmpp):
        ## First order contributions
        D_mu = CFL*Nmp*Sp2[k]*S1[point]/(Sum1*Sump2)
        ## Second order contributions
        if t_order==2:
            D_mu +=0.5*Nmp*(CFL**2)*( S1[point]/(Sum1*Sump2)*(1-(Nmpp)*Sp2[k]/Sump2) -(Sp2[k]/Sump2)*(S1[point]/Sum1-S2[point]/Sum2) )
        Res=Res + np.abs(D_mu)    
    Residual = lambdify((CFL),Res-1.)
    return Residual


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

def Rand(start, end, num): 
    res = [] 
  
    for j in range(num): 
        res.append(random.randint(start, end)) 
  
    return np.asarray(res)

def RandPosition(numberOfPoints):
    res=[]
    for nPoints in(numberOfPoints):
        position=np.zeros((nPoints))
        for i in range(nPoints):
            position[i]=random.uniform(0., 1.)
        res.append(position)
    return res

# Symbolic function to evaluate shape functions
shape_functions=lambda x: np.matrix([(1-x)/DX,x/DX])

xn = np.array([0.,1.])
DX = 1.

## required for plotting residual
CFL=np.linspace(0.,1.,100.)

samples=1000
number_prev = Rand(1, 4, samples)
position_prev = RandPosition(number_prev)

number_curr = Rand(1, 4, samples)
position_curr = RandPosition(number_curr)


if not os.path.exists('eulerRandom.npy'):
    eulerSolution=[]
    rk2Solution=[]
    eulerSolution_id=[]
    rk2Solution_id=[]
    for i in range(samples):
        print "Computing critical CFL for sample ",i,": ",number_curr[i]," particles"
        shapes_prev=shape_functions(position_prev[i])
        shapes_curr=shape_functions(position_curr[i])
        solution_euler=[]
        solution_rk2=[]
        solution_euler_id=[]
        solution_rk2_id=[]
        for k in range(number_curr[i]):
            # if number_curr[i]<number_prev[i] :
            #     print "Attention ca va merder !!!!!!"
            # else:
            #     print "Ca va le faire..."
        
            res=residual(k,position_curr[i],position_prev[i],1)
            solution_euler.append(gridSearch(res))
            res=residual(k,position_curr[i],position_curr[i],1)
            solution_euler_id.append(gridSearch(res))
            res=residual(k,position_curr[i],position_prev[i],2)
            solution_rk2.append(gridSearch(res))
            res=residual(k,position_curr[i],position_curr[i],2)
            solution_rk2_id.append(gridSearch(res))

        eulerSolution.append(min(solution_euler))
        rk2Solution.append(min(solution_rk2))
        eulerSolution_id.append(min(solution_euler_id))
        rk2Solution_id.append(min(solution_rk2_id))
    np.save('eulerRandom.npy',eulerSolution)
    np.save('rk2Random.npy',rk2Solution)
    np.save('eulerRandom_id.npy',eulerSolution_id)
    np.save('rk2Random_id.npy',rk2Solution_id)
else :
    eulerSolution=np.load('eulerRandom.npy')
    rk2Solution=np.load('rk2Random.npy')
    eulerSolution_id=np.load('eulerRandom_id.npy')
    rk2Solution_id=np.load('rk2Random_id.npy')

import statistics
print "Mean CFL for euler periodic: ", statistics.mean(eulerSolution_id)
print "Mean CFL for euler non-periodic: ", statistics.mean(eulerSolution)
print "Mean CFL for rk2 periodic: ", statistics.mean(rk2Solution_id)
print "Mean CFL for rk2 non-periodic: ", statistics.mean(rk2Solution)
print " "
print "Median CFL for euler periodic: ", statistics.median(eulerSolution_id)
print "Median CFL for euler non-periodic: ", statistics.median(eulerSolution)
print "Median CFL for rk2 periodic: ", statistics.median(rk2Solution_id)
print "Median CFL for rk2 non-periodic: ", statistics.median(rk2Solution)
pdb.set_trace()
barsEuler=np.histogram(eulerSolution,bins=np.linspace(0.,1.,11))
barsRk2=np.histogram(rk2Solution,bins=np.linspace(0.,1.,11))

export2DTeXFile('cflStatistics.tex',barsEuler[1],np.array([barsEuler[0]/float(samples),barsRk2[0]/float(samples)]),['Euler','RK2'])

barsEuler2=np.histogram(eulerSolution_id,bins=np.linspace(0.,1.,11))
barsRk22=np.histogram(rk2Solution_id,bins=np.linspace(0.,1.,11))
pdb.set_trace()
export2DTeXFile('cflStatistics_id.tex',barsEuler2[1],np.array([barsEuler2[0]/float(samples),barsRk22[0]/float(samples)]),['Euler','RK2'])
