#!/usr/bin/python

import numpy as np

def Lin_shp(idx,xi,Xe):
    if (idx == 1):
        return np.array([(1.0-xi)/2.0,(1.0+xi)/2.0])
    elif (idx == 2):
        return (np.array([-np.ones(len(xi)),np.ones(len(xi))])/(Xe[1]-Xe[0]))

def Lin_eval(xi,Xe):
    return ((Xe[0]*(-xi+1.0)+Xe[1]*(xi+1.0))/2.0)

def Cst_eval(xi,Xe):
    return Xe*np.ones(len(xi))

def Jac(Xe):
    return ((Xe[1]-Xe[0])/2.0)

def Linear_form(x,op,lin_fields,cst_fields):
    Nnodes = len(x)
    Nelem = Nnodes - 1
    T10 = np.array([np.arange(0,Nnodes-1,1),np.arange(1,Nnodes,1)]).T
    h = x[1:len(x)]-x[:len(x)-1]
    
    xi_integ = np.array([0.774596669241483,0.0,-0.774596669241483])
    w_integ = np.array([0.555555555555556,0.888888888888889,0.555555555555556])
    A = np.zeros(Nnodes)
    for i in range(Nelem):
        mapp = T10[i,:]
        Xe = x[mapp]
        loc_jac = Jac(Xe)
        
        M = Lin_shp(op+1,xi_integ,Xe)
        D = np.array(w_integ)
        if len(lin_fields)!=0:
            D*=Lin_eval(xi_integ,lin_fields[mapp])
        if len(cst_fields)!=0:
            D*=Cst_eval(xi_integ,cst_fields[i])
        A[mapp[0]:mapp[1]+1] += (loc_jac*np.dot(M,D.T))
    return A
