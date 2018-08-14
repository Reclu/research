#!/usr/bin/python

import numpy as np

def bilinear(x,u_n,u,EPn,Pn,E,Sigy,H):
    #initialization
    h = x[1:len(x)]-x[:(len(x)-1)]
    eps_n = (u_n[1:len(u_n)]-u_n[:(len(u_n)-1)])/h
    eps = (u[1:len(u)]-u[:(len(u)-1)])/h
    S = np.zeros(len(eps))
    EP = np.zeros(len(eps))
    P  = np.zeros(len(eps))
    TM  = np.zeros(len(eps))
    #Loop on integration points
    for i,DEFO in enumerate(eps):
        #(i) Elastic prediction
        Selas = E*(DEFO-EPn[i])
        #(ii) Compute the criterion 
        f = np.abs(Selas) - (Sigy+H*Pn[i])
        if (f<=0):
            #elastic step
            S[i] = Selas
            EP[i] = EPn[i]
            P[i] = Pn[i]
            TM[i] = E
        elif (f>0):
            #elastoplastic step: solve a nonlinear scalar equation
            dP = f/(E+H)
            P[i] = Pn[i]+dP
            EP[i] = EPn[i]+(P[i]-Pn[i])*np.sign(Selas)
            S[i] = E*(DEFO-EP[i])
            TM[i] = (E*H)/(E+H)
    return S,P,EP,TM
