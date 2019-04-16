# !/usr/bin/python

import numpy as np


class DGmesh:
    #Constructor
    def __init__(self,Mp,l,ppc,c,rho):
        nex = Mp/ppc
        nnx=2*nex+2
        lmp=l/(Mp-1)
        self.dx = ppc*l/nex
        self.xn = np.linspace(-lmp/2,l+lmp/2,nex+1)
        self.connect = np.array([np.arange(1,nnx-1,2),np.arange(2,nnx,2)]).T
        self.c = c
        self.rho = rho
        
        
    def buildApproximation(self,xp):
        xn =self.xn
        connect=self.connect
        Nn=2*len(connect) + 2
        Np=np.shape(xp)[0]
        Dofs=np.zeros((Nn))
        Parent=np.zeros(Np)
        matpoints=np.zeros(np.shape(connect)[0])
        self.Map=np.zeros((Nn,Np))
        self.Grad=np.zeros((Nn,Np))
        Le=(xn[1]-xn[0])
        for Pt in range(Np):
            # detect parent element of current material point
            parent=np.ceil((xp[Pt,0]-xn[0])/Le)-1 
            d=np.array([connect[int(parent),:]]).astype(np.int64)
            # Active nodes' indices storage
            Dofs[d[0,:]]+=d[0,:]
            xi=(xp[Pt,0]-xn[(d[0,0]-1)/2])/Le;
            N=np.array([1-xi,xi])
            dN=np.array([-1/Le,1/Le])
            self.Map[d,Pt]+=N.T
            self.Grad[d,Pt]+=dN.T
            # parent element storage
            Parent[Pt]=parent
        self.d=[]
        for j in range(np.shape(Dofs)[0]):
            if Dofs[j]!=0:
                self.d.append(j)
        return self.Map,self.Grad,self.d,Parent.astype(int)


    def computeFlux(self,U,dofs,md,dim,limiter):
        """
        Evaluation of interface fluxes by solving Riemann problems
        Computation of volumic fluxes over each cell
        Assembling of interfaces and volumes contributions
        """
        if dim==2:
            fluxes = self.computeInterfaceFlux(U,md)
            fint = self.computeInternalForces(U,dofs)
            f = self.computeTotalForces(fint,fluxes,U)
        if dim==1:
            if limiter!=-1:
                U = self.slopeLimiter(U,limiter)
            # Compute interface fluxes :
            Nnodes = np.shape(U)[0]
            Nelem = (Nnodes-2)/2
            #amdq = np.zeros(Nelem+1)
            apdq = np.zeros(Nelem+1)
            T10 = np.array([np.arange(1,Nnodes-1,2),np.arange(2,Nnodes,2)]).T
            #Loop on interfaces
            for i in range(Nelem+1):
                UL = U[2*i]
                apdq[i] =  self.c*UL
            ### Compute volumic term
            S = self.Kx
            Nnodes = np.shape(U)[0]
            fint = np.zeros(Nnodes)
            flux = self.c*U
            fint[dofs] = np.dot(S,flux[dofs])
            f = np.copy(fint)
            for i in range(Nelem):
                mapp = T10[i,:]
                mapp_interfaces = mapp - i -1
                f[mapp[0]]+=self.rho*apdq[mapp_interfaces[0]] # left node
                f[mapp[1]]-=self.rho*apdq[mapp_interfaces[1]] # right node
        return f

    def computeDelta(self,dU):
        delta1 = (dU[0]+self.c*dU[1])/(2*self.c)
        delta2 = -(dU[0]-self.c*dU[1])/(2*self.c)
        return [delta1,delta2]

    def computeElasticWaves(self,delta):
        waves = np.zeros((2,2))
        waves[:,0] = delta[0]*np.array([self.c,1])
        waves[:,1] = delta[1]*np.array([-self.c,1])
        return waves

    def computeInterfaceFlux(self,U,md):
        Nnodes = np.shape(U)[0]
        Nelem = (Nnodes-2)/2
        fluxes = np.zeros((Nelem+1,U.shape[1]))
        #Loop on interfaces
        for i in range(Nelem+1):
            UL = U[2*i,:]
            UR = U[2*i+1,:]
            dU = UR-UL
            alpha = self.computeDelta(dU)
            waves = self.computeElasticWaves(alpha)
            A=np.array([[0.,-self.c**2],[-1.,0.]])
            fluxes[i,:] = np.dot(A,UL) -self.c*waves[:,0]
        return fluxes

    def computeTotalForces(self,fint,fluxes,U):
        Nnodes = np.shape(fint)[0]
        Nelem = (Nnodes-2)/2
        # Loop on elements
        f=np.copy(fint) # Weak
        for i in range(Nelem):
            mapp = self.connect[i,:]
            mapp_inter = mapp - i -1
            f[mapp[0],:]+=self.rho*fluxes[mapp_inter[0],:] # Left node
            f[mapp[1],:]-=self.rho*fluxes[mapp_inter[1],:] # Right node
        return f
    
    def computeFluxVector(self,U):
        if dim==2:
            Nnodes = np.shape(U)[0]
            flux = np.zeros((Nnodes,2))
            for i in range(Nnodes):
                flux[i,:] = np.array([-self.c**2*U[i,1],-U[i,0]])
        if dim==1:
            flux = self.c*U
        return flux

    def computeInternalForces(self,U,dofs):
        S = self.Kx
        Nnodes = np.shape(U)[0]
        fint = np.zeros((Nnodes,2))
        flux = self.computeFluxVector(U)
        fint[dofs,0] = np.dot(S,flux[dofs,0])
        fint[dofs,1] = np.dot(S,flux[dofs,1])
        return fint

    def setMapping(self,Kx):
        self.Kx = Kx
        

    ## METHODS FOR LIMITERS

    def minmod(self,*args):
        s=0.
        for i in range(len(args)):
            s+=np.sign(args[i])
        s/=len(args)
        if (np.abs(s)==1) :
            return s*np.min(np.abs(args))
        else :
            return 0.

    def maxmod(self,*args):
        s=0.
        for i in range(len(args)):
            s+=np.sign(args[i])
        s/=len(args)
        if (np.abs(s)==1) :
            return s*np.max(np.abs(args))
        else :
            return 0.
        
    def slopeLimiter(self,U,limiter):
        dx = self.dx
        Nnodes = len(U)
        Nelem = Nnodes/2 -1
        for i in range(Nelem):
            # Build mean value vector
            Umean=0.5*(U[2*i+1]+U[2*i+2])
            if i==0 : UmeanL= U[0]
            else : UmeanL=0.5*(U[2*i-1]+U[2*i])
            if i==Nelem-1 : UmeanR= U[-1]
            else : UmeanR=0.5*(U[2*i+3]+U[2*i+4])
            dU = (U[2*i+2]-U[2*i+1])/dx
            dUL = (Umean-UmeanL)/dx
            dUR = (UmeanR-Umean)/dx
            if limiter==0 : s = self.minmodLimiter(dUL,dUR)
            elif limiter==1 : s = self.superbeeLimiter(dU,dUL,dUR)
            elif limiter==2 : s = self.musclLimiter(dU,dUL,dUR)
            elif limiter==3 : s = self.mcLimiter(dU,dUL,dUR)
            ul = Umean - 0.5*dx*s
            ur = Umean + 0.5*dx*s
            if not (np.abs(ul-U[2*i+1])<1.e-8 and np.abs(ur-U[2*i+2])<1.e-8):
                U[2*i+1]=ul
                U[2*i+2]=ur
        return U

    
    def mcLimiter(self,dU,dUL,dUR):
        return self.minmod(dU,2.*dUR,2.*dUL)

    def minmodLimiter(self,dUL,dUR):
        return self.minmod(dUR,dUL)
    
    def superbeeLimiter(self,dU,dUL,dUR):
        sl = self.minmod(dUR,2.*dUL,dU)
        sr = self.minmod(2.*dUR,dUL,dU)
        si = self.minmod(dUR,dUL,2.*dU)
        return self.maxmod(sl,si,sr)
    
    def musclLimiter(self,dU,dUL,dUR):
        return self.minmod(dU,dUR,dUL)


