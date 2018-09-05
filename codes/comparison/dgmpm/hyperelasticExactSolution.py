#!/usr/bin/python

import numpy as np
import scipy.optimize

def compute_stationnary_solution(qL,qR,E,rho):
    JL = qL[0] ; JR = qR[0]
    vL = qL[1]/rho ; vR = qR[1]/rho
    rvL = qL[1] ; rvR = qR[1]
    if qL[0]-qR[0]>1.e-12:
        # 1-rarefaction-2-shock
        ## Residual form : J - (f(J) - g(J,JR) + A(JL) )/h(J) 
        R=lambda x: computeResidual(rho,E,rvL,rvR,JR,JL,x)
        # python root finder within an interval
        J_star = scipy.optimize.brentq(R,JR,JL)
        rv_star = computeQ2_star_1R2S(E,rho,rvL,JL,J_star)
    elif qR[0]-qL[0]>1.e-12:
        # 1-shock-2-rarefaction
        ## Residual form : J - (f(J) - g(J) + A )/h(J) 
        R=lambda x: computeResidual(rho,E,rvL,rvR,JL,JR,x)
        # python root finder within an interval
        J_star = scipy.optimize.brentq(R,JL,JR)
        rv_star = computeQ2_star_1S2R(E,rho,rvR,JR,J_star)
    else :
        # Simple wave
        fL = compute_F(JL) 
        hL = compute_H(JL) 
        rv_star = 0.5*(rvL+rvR) 
        A = 0.5*(rvR-rvL)*np.sqrt(8./(E*rho))
        R = lambda x : x*compute_H(x) - compute_F(x)/np.sqrt(3.) - JL*hL + fL/np.sqrt(3.)-A
        dR = lambda x : compute_H(x) + x*compute_Hprime(x) - compute_Fprime(x)/np.sqrt(3.)
        J_star = scipy.optimize.newton(R,JL,dR)
    return [J_star,rv_star]


def computeSolution(t,X,qL,qR,q_star,E,rho,impact_coor=[]):
    if qL[0]-qR[0]>1.e-12:
        # 1-rarefaction-2-shock
        #print q_star[0],q_star[1]/rho
        s = -(qR[1]/rho - q_star[1]/rho)/(qR[0]-q_star[0])
        if X/t - s >1.e-12 :
            # Right solution
            J = qR[0]
            v = qR[1]/rho
        else :
            # Stationnary solution
            J = q_star[0]
            v = q_star[1]/rho
    elif qR[0]-qL[0]>1.e-12:
        # 1-shock-2-rarefaction
        cL =np.sqrt(E*(3*qL[0]**2 -1)/(2.*rho))
        cR =np.sqrt(E*(3*qR[0]**2 -1)/(2.*rho))
        cstar =np.sqrt(E*(3*q_star[0]**2 -1)/(2.*rho))
        if X/t - cR > 0. :#1.e-12 :
            # Right solution
            J = qR[0]
            v = qR[1]/rho
        elif X/t - cstar < 0. :#1.e-12 :
            # Stationnary solution
            J = q_star[0]
            v = q_star[1]/rho
        else :
            # Inside rarefaction fan
            J = np.sqrt( 2.*rho*(X/t)*(X/t)/E + 1.)/np.sqrt(3.)
            v = qR[1] - np.sqrt(E*rho/24.)*( np.sqrt(3.)*(J*compute_H(J) - qR[0]*compute_H(qR[0]) ) - compute_F(J) + compute_F(qR[0]))
            v/=rho
    else :
        cL =np.sqrt(E*(3*qL[0]**2 -1)/(2.*rho))
        cR =np.sqrt(E*(3*qR[0]**2 -1)/(2.*rho))
        
        if (X > impact_coor - cL*t) :
            if (X < impact_coor + cR*t):
                # Stationnary solution
                J = q_star[0]
                v = q_star[1]/rho
            else :
                # Right solution
                J = qR[0]
                v = qR[1]/rho
        else :
            # Left solution
            J = qL[0]
            v = qL[1]/rho
    pi = E*J*(J**2-1.)/2.
    return [pi,v]

def compute_F(J):
    F = np.log( np.sqrt(3.)*J + np.sqrt(3.*J*J - 1.) )
    return F

def compute_Fprime(J):
    F = compute_F(J)
    Fp= (np.sqrt(3.) + 3.*J/np.sqrt(3.*J*J -1.) )/np.exp(F)
    return Fp

def compute_G(J,Jside):
    G = np.sqrt( (Jside-J)*(Jside**3-Jside-(J**3-J)) )
    return G

def compute_Gprime(J,Jside):
    G = compute_G(J,Jside)
    Gp = ((4*J-3*Jside)*J**2 + 2.*(Jside-J) - Jside**3)/(2*G)
    return Gp

def compute_H(J):
    H = np.sqrt( 3.*J**2 -1.)
    return H

def compute_Hprime(J):
    H = compute_H(J)
    Hp = 3.*J/H
    return Hp

def compute_A(rho,E,rvL,rvR,J):
    A=2.*np.sqrt(2./(E*rho))*(rvR-rvL) - compute_F(J)/np.sqrt(3.) + J*compute_H(J)
    return A
                    
def computeResidual(rho,E,rvL,rvR,J,J_A,Ji):
    Res = Ji-(compute_F(Ji)/np.sqrt(3.) - 2*compute_G(Ji,J) + compute_A(rho,E,rvL,rvR,J_A) )/compute_H(Ji)
    return Res

def compute_dRes(rho,E,rvL,rvR,J,J_A,Ji):
    f = compute_F(Ji) ; df = compute_Fprime(Ji)
    g = compute_G(Ji,J) ; dg = compute_Gprime(Ji,J)
    h = compute_H(Ji) ; dh = compute_Hprime(Ji)
    A = compute_A(rho,E,vL,vR,J_A)
    dR = 1 - (df/np.sqrt(3) -2*dg)/h + (f/np.sqrt(3.) -2.*g +A)*dh/(h**2)
    return dR

def computeQ2_star_1R2S(E,rho,rvL,JL,Ji):
    q2 = rvL + np.sqrt(E*rho/24.)*(np.sqrt(3.)*(Ji*compute_H(Ji) - JL*compute_H(JL) ) \
                    - compute_F(Ji) + compute_F(JL))
    return q2

def computeQ2_star_1S2R(E,rho,rvR,JR,Ji):
    q2 = rvR - np.sqrt(E*rho/24.)*(np.sqrt(3.)*(Ji*compute_H(Ji) - JR*compute_H(JR) ) \
                    - compute_F(Ji) + compute_F(JR))
    return q2

def celerity(E,rho0,J):
    return np.sqrt(E*(3.*J**2 -1.)/(2.*rho0))

def applyBoundaryCondition(E,rho,Jd,qR):
    JR = qR[0] ; rvR=qR[1]
    cstar = celerity(E,rho,Jd)
    cR = celerity(E,rho,JR)
    hyperbolicity= np.sqrt(1/3.)
    if cstar<cR :
        print "Boundary condition based on a 2-rarefaction"
        # 2-rarefaction
        rhs = JR*compute_H(JR) - Jd*compute_H(Jd) + (compute_F(Jd)-compute_F(JR))/np.sqrt(3.)
        R = lambda x: x - Jd - 0.25*rhs**2/(x**3 - x -Jd**3 + Jd)
        dR = lambda x: 1 + (0.25*(3.*x**2 - 1.)*rhs**2)/((x**3 - x -Jd**3 + Jd)**2)
        JL = scipy.optimize.newton(R,hyperbolicity,dR)
    else :
        print "Boundary condition based on a 2-shock"
        # 2-shock
        rhs = Jd*compute_H(Jd) - compute_F(Jd)/np.sqrt(3.) + 2*compute_G(Jd,JR)
        R=lambda x: x - (compute_F(x)/np.sqrt(3.) + rhs)/compute_H(x)
        dR=lambda x: 1. - (compute_Fprime(x)/np.sqrt(3.))/compute_H(x) + (compute_F(x)/np.sqrt(3.) + rhs)*compute_Hprime(x)/compute_H(x)**2
        JL = scipy.optimize.newton(R,JR,dR)
    return [JL,rvR]
