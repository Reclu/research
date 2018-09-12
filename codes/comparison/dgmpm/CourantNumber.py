from sympy import *
import numpy as np
import scipy.optimize as optimize

def computeCourantNumber(mesh,parent,Map,t_order):
    print "Computing Courant number ..."
    sol=[]
    # sum over mesh cells to compute the lowest Courant number by solving the minimum amplification factor equation

    for i in range(np.shape(mesh.connect)[0]-5):
        # Number of material points in cell i
        matpointN = np.where(parent==i)[0]
        N = len(matpointN)
        # Number of material points in cell i-1
        matpointM = np.where(parent==i-1)[0]
        M = len(matpointM)
        # nodes indices
        n1 = 2*i+1 ; n2=2*i+2
        S1 = Map[n1,matpointN] ; S2 = Map[n2,matpointN]
        sum1=np.sum(S1) ; sum2=np.sum(S2)
        if N!=0:
            if M!=0 :
                S1M = Map[n1-2,matpointM] ; S2M= Map[n1-1,matpointM] 
                sum1M=np.sum(S1M) ; sum2M=np.sum(S2M)
            CFL = symbols('CFL')
            for alpha in range(N):
                Res=0.
                for mu in range(N): # sum over particles contained in the same cell
                    ## First order contributions
                    Dmu=S1[alpha]*S1[mu]/sum1 + S2[alpha]*S2[mu]/sum2
                    Dmu+=CFL*(S2[alpha]/sum2*(1.-N*S2[mu]/sum2) - S1[alpha]/sum1)
                    ## Second order contributions
                    if t_order==2:
                        Dmu+=0.5*N*(CFL**2)*((S2[mu]/sum2)*(S1[alpha]/sum1-S2[alpha]/sum2) + (S2[alpha]/sum2)*(N*S2[mu]/sum2-1.)/sum2)
                    Res+=np.abs(Dmu)
                for mu in range(M): # sum over particles contained in the same cell
                    ## First order contributions
                    Dmu=CFL*S1[alpha]*N*S2M[mu]/(sum1*sum2M)
                    ## Second order contributions
                    if t_order==2:
                        Dmu +=0.5*N*(CFL**2)*( S1[alpha]/(sum1*sum2M)*(1.-M*S2M[mu]/sum2M) -(S2M[mu]/sum2M)*(S1[alpha]/sum1-S2[alpha]/sum2) )    
                    Res+=np.abs(Dmu)
                Residual = lambdify((CFL),Res-1.)
                solution=optimize.newton(Residual,1.)
                sol.append(solution)
    solution_CFL=np.min(sol)
    print "Courant number set to :",solution_CFL
    return solution_CFL
