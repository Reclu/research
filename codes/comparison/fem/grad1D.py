#!/usr/bin/python
import numpy as np

def Grad1D(x,u):
    h = x[1:len(x)]-x[:(len(x)-1)]
    wg = 2.0
    grad_u_g = (u[1:len(u)]-u[:(len(u)-1)])/h
    #grad_u_x = np.concatenate((np.array([grad_u_g[0]]),(wg*grad_u_g[1:len(u)-1]+wg*grad_u_g[:(len(u)-2)])/(2.0*wg),np.array([grad_u_g[-1]])),axis=0)
    return grad_u_g
