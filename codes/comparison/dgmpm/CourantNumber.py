from sympy import *
import numpy as np
import scipy.optimize as optimize
import pdb

def computeCourantNumber(mesh,parent,Map,t_order):
    print "Computing Courant number ..."
    sol=[]
    # Symbolic function to evaluate shape functions
    # sum over mesh cells to compute the lowest Courant number by solving the minimum amplification factor equation
    ###############################"Residual = residualEuler(0,shapes,shapes)
    mat_points=np.shape(Map)[1]
    ## sum over material points
    for alpha in range(mat_points):
        ## parent cell and nodes indices
        cell=parent[alpha]
        n1=2*cell+1 ; n2=2*cell+2
        ## material points in cell
        matpointN = np.where(parent==cell)[0]
        ## material points in previous cell
        matpointP = np.where(parent==cell-1)[0]
        ## particle indice in current cell
        point=np.where(matpointN==alpha)[0][0]
        if len(matpointP)==0:
            ## if no particle in previous cell, take the next one
            matpointP = np.where(parent==cell+1)[0]
        ## shape functions in cell i
        S1 = Map[n1,matpointN] ; S2 = Map[n2,matpointN]
        ## shape functions in cell i-1
        S1M = Map[n1,matpointN] ; S2M = Map[n2,matpointN]
        if t_order==1:
            residual=residualEuler(point,[S1,S2],[S1M,S2M])
        elif t_order==2:
            residual=residualRK2(point,[S1,S2],[S1M,S2M])
        solution=optimize.newton(residual,1.)
        sol.append(solution)
    solution_CFL=np.min(sol)
    print "Courant number set to :",solution_CFL
    return solution_CFL

def residualRK2(point,S,Sp):
    CFL = symbols('CFL')
    Res=0.
    ## Curent cell shape functions
    S1=S[0];S2=S[1]
    Sum1=np.sum(S1) ; Sum2=np.sum(S2)
    Nmp=len(S1) 
    ## Previous cell shape functions
    Sp1=Sp[0];Sp2=Sp[1]
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
    ## Curent cell shape functions
    S1=S[0];S2=S[1]
    Sum1=np.sum(S1) ; Sum2=np.sum(S2)
    Nmp=len(S1) 
    ## Previous cell shape functions
    Sp1=Sp[0];Sp2=Sp[1]
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
