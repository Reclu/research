#!/usr/bin/python

import numpy as np
import math as m

def buildApproximation1D(xp,xn,connect):
    Nn=len(xn)
    Np=np.shape(xp)[0]
    Dofs=np.zeros((Nn))
    Parent=np.zeros(Np)
    matpoints=np.zeros(np.shape(connect)[0])
    Map=np.zeros((Nn,Np))
    Grad=np.zeros((Nn,Np))
    Le=(xn[1]-xn[0])
    for Pt in range(Np):
        # detect parent element of current material point
        parent=m.ceil((xp[Pt,0]-xn[0])/Le)-1 
        d=np.array([connect[parent,:]]).astype(np.int64)
        # Active nodes' indices storage
        Dofs[d[0]]+=1
        Parent[Pt]=parent+1
        xi=(xp[Pt,0]-xn[d[0][0]])/Le;
        N=np.array([1-xi,xi])
        dN=np.array([-1/Le,1/Le])
        Map[d,Pt]+=N.T;
        Grad[d,Pt]+=dN.T;
    # Number of material point in each element to compute the particles' volumes
    for elem in range(np.shape(connect)[0]):
        matpoints[elem]=np.shape(Parent[Parent==elem+1])[0]
    d=[]
    for j in range(np.shape(Dofs)[0]):
        if Dofs[j]!=0:
            d.append(j)
    return Map,Grad,d

def buildGIMPApproximation1D(xp,xn,lp):
    Nn=len(xn)
    Np=np.shape(xp)[0]
    Dofs=np.zeros((Nn)) 
    Map=np.zeros((Nn,Np))
    Grad=np.zeros((Nn,Np))
    le=xn[1]-xn[0]
    d=[]
    for Pt in range(Np):
        xi=xp[Pt,0]
        #####################################
        # Particular cases for boundary nodes
        x=xn[0]
        Map[0,Pt]=(1-(xi-x)/le)*(((xi-x)>-lp and(xi-x)<lp)or((xi-x)>=lp and(xi-x)<=le-lp))+\
            ((lp+le-xi+x)**2/(4*lp*le))*((xi-x)>(le-lp)and(xi-x)<=lp+le)
        Grad[0,Pt]=-1/le*( ((xi-x)>-lp and(xi-x)<lp) or ((xi-x)>=lp and(xi-x)<=le-lp) ) -\
            (le+lp+x-xi)/(2*le*lp)*((xi-x)>(le-lp)and(xi-x)<=lp+le)
        x=xn[Nn-1]
        Map[Nn-1,Pt]=(1+(xi-x)/le)*( ((xi-x)>=(-le+lp)and(xi-x)<=-lp) or ((xi-x)>-lp and(xi-x)<lp) ) +\
            (le+lp+xi-x)**2/(4*le*lp)*((xi-x)>=-(le+lp)and(xi-x)<(-le+lp))
        Grad[Nn-1,Pt]=(x-xi)/(le*lp)*( ((xi-x)>=(-le+lp)and(xi-x)<=-lp) or ((xi-x)>-lp and(xi-x)<lp) ) +\
            (le+lp+(xi-x))/(2*le*lp)*((xi-x)>=(-le-lp)and(xi-x)<(-le+lp))
        #####################################
        for N in range(1,Nn-1):
            x=xn[N]
            Map[N,Pt]=(le+lp+xi-x)**2/(4*le*lp)*((xi-x)>=-(le+lp) and (xi-x)<(-le+lp))+\
                (1+(xi-x)/le)*( (xi-x)>=(-le+lp) and (xi-x)<=-lp )+\
                (1-((xi-x)**2+lp**2)/(2*lp*le))*( (xi-x)>-lp and (xi-x)<lp )+\
                (1-(xi-x)/le)*( (xi-x)>=lp and (xi-x)<=le-lp )+\
                ((lp+le-xi+x)**2/(4*lp*le))*( (xi-x)>(le-lp) and (xi-x)<=lp+le )        

            Grad[N,Pt]=(le+lp+(xi-x))/(2*le*lp)*((xi-x)>=(-le-lp)and(xi-x)<(-le+lp))+\
                1/le*((xi-x)>=(-le+lp)and(xi-x)<=-lp)+\
                (x-xi)/(le*lp)*((xi-x)>-lp and(xi-x)<lp)-\
                1/le*((xi-x)>=lp and(xi-x)<=le-lp)-\
                (le+lp+x-xi)/(2*le*lp)*((xi-x)>(le-lp)and(xi-x)<=lp+le) 
    for N in range(Nn):
        if not all(Map[N,:]==0): # node is added if its line in Mapping is not zero (i.e :it's connected to a MP)
            d.append(N)
    return Map,Grad,d

def buildApproximation2D(xp,xn,connect,nnx):
    nex=nnx-1
    Nn=np.shape(xn)[0]
    Np=np.shape(xp)[0]
    Map=np.zeros((2*Nn,2*Np))
    Grad=np.zeros((2*Nn,3*Np))
    Parent=np.zeros(Np)
    Dofs=np.zeros((Nn))
    matpoints=np.zeros(np.shape(connect)[0])
    Volume=np.zeros(Np)
    le=xn[1,0]-xn[0,0]
    he=le
    for Pt in range(Np):
        # detect parent element of current material point
        parent=m.ceil((xp[Pt,0]-xn[0,0])/le)-1 + (m.ceil((xp[Pt,1]-xn[0,1])/he))*nex-1
        d=np.array([connect[parent-1]]).astype(np.int64)
        # Active nodes' indices storage
        Dofs[d[0,:]-1]+=d[0,:]
        # Parent element storage (to compute each particle's volume)
        Parent[Pt]=parent
        ndof=np.array([2*d[0,0]-2,2*d[0,0]-1,2*d[0,1]-2,\
                       2*d[0,1]-1,2*d[0,2]-2,2*d[0,2]-1,\
                       2*d[0,3]-2,2*d[0,3]-1])
        ndof=ndof.astype(np.int64)
        pdof=np.array([2*Pt,2*Pt+1])
        pdof=pdof.astype(np.int64)
        pgdof=np.array([3*Pt,3*Pt+1,3*Pt+2]);
        pgdof=pgdof.astype(np.int64)
        xi=2*(xp[Pt,0]-xn[d[0,0]-1,0])/le - 1;
        et=2*(xp[Pt,1]-xn[d[0,0]-1,1])/he - 1;
        # Shape functions and derivatives
        N=np.array([(1-xi)*(1-et)/4,(1+xi)*(1-et)/4,(1+xi)*(1+et)/4,(1-xi)*(1+et)/4])
        dN=np.mat([[-(1-et)/4,-(1-xi)/4],\
                   [ (1-et)/4,-(1+xi)/4],\
                   [ (1+et)/4, (1+xi)/4],\
                   [-(1+et)/4, (1-xi)/4]])
        # Elementary mapping and gradient matrix
        S=np.array([[N[0],0.,N[1],0.,N[2],0.,N[3],0.],\
                  [0.,N[0],0.,N[1],0.,N[2],0.,N[3]]])              
        Coord=np.mat([ xn[d[0,:]-1,0].T , xn[d[0,:]-1,1].T])
        Jac=Coord*dN;
        Jinv=np.linalg.inv(Jac)
        dN*=Jinv
        B1=np.array([dN[0,0],0.,dN[1,0],0.,dN[2,0],0.,dN[3,0],0.])
        B2=np.array([0.,dN[0,1],0.,dN[1,1],0.,dN[2,1],0.,dN[3,1]])
        B3=np.array([dN[0,1],dN[0,0],dN[1,1],dN[1,0],dN[2,1],dN[2,0],dN[3,1],dN[3,0]])
        B=np.array([ B1.T,B2.T,B3.T ])
        # Global mapping and gradient matrix
        Map[ndof,pdof[0]]+=S[0];Map[ndof,pdof[1]]+=S[1]
        Grad[ndof,pgdof[0]]+=B[0];Grad[ndof,pgdof[1]]+=B[1];Grad[ndof,pgdof[2]]+=B[2]
    # Number of material point in each element to compute the particles' volumes
    for elem in range(np.shape(connect)[0]):
        matpoints[elem]=np.shape(Parent[Parent==elem+1])[0]
    # Volume of each material point
    for Pt in range(Np):
        Volume[Pt]=he*le/matpoints[Parent[Pt]-1]
    d=[]
    for j in range(np.shape(Dofs)[0]):
        if Dofs[j]!=0:
            d.append(2*j)
            d.append(2*j+1)
    return Map,Grad,Volume,d

def buildGIMPApproximation2D(xp,xn,connect,nnx):
    nex=nnx-1
    Nn=np.shape(xn)[0]
    Np=np.shape(xp)[0]
    Map=np.zeros((2*Nn,2*Np))
    Grad=np.zeros((2*Nn,3*Np))
    Parent=np.zeros(Np)
    Dofs=np.zeros((Nn))
    matpoints=np.zeros(np.shape(connect)[0])
    Volume=np.zeros(Np)
    le=xn[1,0]-xn[0,0]
    he=le
    for Pt in range(Np):
        # detect parent element of current material point
        parent=m.floor(xp[Pt,0]/le)+1 + (m.floor(xp[Pt,1]/he))*nex
        d=np.array([connect[parent-1]]).astype(np.int64)
        # Active nodes' indices storage
        Dofs[d[0,:]]+=d[0,:]
        K=Dofs[Dofs!=0].astype(np.int64)
        # Parent element storage (to compute each particle's volume)
        Parent[Pt]=parent
        ndof=np.array([2*d[0,0]-2,2*d[0,0]-1,2*d[0,1]-2,\
                       2*d[0,1]-1,2*d[0,2]-2,2*d[0,2]-1,\
                       2*d[0,3]-2,2*d[0,3]-1])
        ndof=ndof.astype(np.int64)
        pdof=np.array([2*Pt,2*Pt+1])
        pdof=pdof.astype(np.int64)
        pgdof=np.array([3*Pt,3*Pt+1,3*Pt+2]);
        pgdof=pgdof.astype(np.int64)
        xi=2*(xp[Pt,0]-xn[d[0,0]-1,0])/le - 1;
        et=2*(xp[Pt,1]-xn[d[0,0]-1,1])/he - 1;
        # Shape functions and derivatives
        N=np.array([(1-xi)*(1-et)/4,(1+xi)*(1-et)/4,(1+xi)*(1+et)/4,(1-xi)*(1+et)/4])
        dN=np.mat([[-(1-et)/4,-(1-xi)/4],\
                   [ (1-et)/4,-(1+xi)/4],\
                   [ (1+et)/4, (1+xi)/4],\
                   [-(1+et)/4, (1-xi)/4]])
        # Elementary mapping and gradient matrix
        S=np.array([[N[0],0.,N[1],0.,N[2],0.,N[3],0.],\
                  [0.,N[0],0.,N[0],0.,N[0],0.,N[0]]])              
        Coord=np.mat([ xn[d[0,:]-1,0].T , xn[d[0,:]-1,1].T])
        Jac=Coord*dN;
        Jinv=np.linalg.inv(Jac)
        dN*=Jinv
        B1=np.array([dN[0,0],0.,dN[1,0],0.,dN[2,0],0.,dN[3,0],0.])
        B2=np.array([0.,dN[0,1],0.,dN[1,1],0.,dN[2,1],0.,dN[3,1]])
        B3=np.array([dN[0,1],dN[0,0],dN[1,1],dN[1,0],dN[2,1],dN[2,0],dN[3,1],dN[3,0]])
        B=np.array([ B1.T,B2.T,B3.T ])
        # Global mapping and gradient matrix
        Map[ndof,pdof[0]]+=S[0];Map[ndof,pdof[1]]+=S[1]
        Grad[ndof,pgdof[0]]+=B[0];Grad[ndof,pgdof[1]]+=B[1];Grad[ndof,pgdof[2]]+=B[2]
    # Number of material point in each element to compute the particles' volumes
    for elem in range(np.shape(connect)[0]):
        matpoints[elem]=np.shape(Parent[Parent==elem+1])[0]
    # Volume of each material point
    for Pt in range(Np):
        Volume[Pt]=he*le/matpoints[Parent[Pt]-1]
    return Map,Grad,K,Volume
