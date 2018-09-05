# !/usr/bin/python

import numpy as np
import scipy.optimize

class DGmesh:
    #Constructor
    def __init__(self,Mp,l,ppc,rho0,E,exact):
        nex = Mp/ppc
        nnx=2*nex+2
        lmp=l/(Mp-1)
        self.dx = ppc*l/nex
        self.xn = np.linspace(-lmp/2,l+lmp/2,nex+1)
        self.connect = np.array([np.arange(1,nnx-1,2),np.arange(2,nnx,2)]).T
        self.rho0 = rho0
        self.E = E
        self.exact = exact
        
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
        Le=(xn[2]-xn[1])
        for Pt in range(Np):
            # detect parent element of current material point
            parent=np.ceil((xp[Pt,0]-xn[0])/Le)-1 
            #print Pt,xp[Pt,0],xn[0],parent
            d=np.array([connect[int(parent),:]]).astype(np.int64)
            # Active nodes' indices storage
            Dofs[d[0,:]]+=d[0,:]
            #xi=(xp[Pt,0]-xn[(d[0,0]-1)/2])/Le;
            xi=(xp[Pt,0]-xn[(d[0,0])/2])/Le;
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


    def computeFlux(self,U,W,dofs,u0,limit,limiter):
        """
        Evaluation of interface fluxes by solving Riemann problems
        Computation of volumic fluxes over each cell
        Assembling of interfaces and volumes contributions
        """
        if (limit) : 
            U = self.slopeLimiter(U,limiter)
            W = self.slopeLimiter(W,limiter)
        fluxes = self.computeInterfaceFlux(U,W)
        fint = self.computeInternalForces(U,dofs,u0)
        f = self.computeTotalForces(fint,fluxes)
        return f

    def computeDelta(self,dU,cL,cR):
        delta1 = (self.rho0*cR*dU[1]+dU[0])/(self.rho0*(cL+cR))
        delta2 = (self.rho0*cL*dU[1]-dU[0])/(self.rho0*(cL+cR))
        return [delta1,delta2]

    def computeElasticWaves(self,delta,cL,cR):
        waves = np.zeros((2,2))
        waves[:,0] = delta[0]*np.array([self.rho0*cL,1.])
        waves[:,1] = delta[1]*np.array([-self.rho0*cR,1.])
        return waves

    def approximateRiemannSolver(self,WL,WR,JL,JR):
        HL = self.E*(-1.+ 3.*JL*JL)/2.
        HR = self.E*(-1.+ 3.*JR*JR)/2.
        cL = np.sqrt(HL/self.rho0)
        cR = np.sqrt(HR/self.rho0)
        dW = WR-WL
        alpha = self.computeDelta(dW,cL,cR)
        waves = self.computeElasticWaves(alpha,cL,cR)
        Wstar = WL + waves[:,0]
        return Wstar
    
    def computeInterfaceFlux(self,U,W):
        Nnodes = np.shape(U)[0]
        Nelem = (Nnodes-2)/2
        fluxes = np.zeros((Nnodes,U.shape[1]))
        #Loop on interfaces
        for i in range(Nelem+1):
            # Take original quantity vector and not specific one
            WL = W[2*i,:]
            WR = W[2*i+1,:]
            JL = U[2*i,0]*self.rho0
            JR = U[2*i+1,0]*self.rho0
            Pi,v = self.approximateRiemannSolver(WL,WR,JL,JR)
            fluxes[2*i,:] = np.array([-v,-Pi])
            fluxes[2*i+1,:] = np.array([-v,-Pi])
        return fluxes

    def computeTotalForces(self,fint,fluxes):
        Nnodes = np.shape(fint)[0]
        Nelem = (Nnodes-2)/2
        # Loop on elements
        f=np.copy(fint) 
        for i in range(Nelem):
            mapp = self.connect[i,:]
            f[mapp[0],:]+=fluxes[mapp[0],:] # Left node
            f[mapp[1],:]-=fluxes[mapp[1],:] # Right node
        return f
    
    def computeFluxVector(self,U,u0):
        Nnodes = np.shape(U)[0]
        flux = np.zeros((Nnodes,2))
        for i in range(Nnodes):
            vs=U[i,1]/self.rho0
            J=U[i,0]*self.rho0
            flux[i,:] = np.array([-vs,-self.E*J*(J*J -1.)/(2.*self.rho0)])
        return flux

    def computeInternalForces(self,U,dofs,u0):
        S = self.Kx
        Nnodes = np.shape(U)[0]
        fint = np.zeros((Nnodes,2))
        flux = self.computeFluxVector(U,u0)
        fint[dofs,0] = np.dot(S,flux[dofs,0])
        fint[dofs,1] = np.dot(S,flux[dofs,1])
        return fint

    def setMapping(self,Kx):
        self.Kx = Kx
        
   ###################################################################################
   ###################################################################################
   ## METHODS FOR EXACT RIEMANN SOLVER

    def compute_F(self,J):
        F = np.log( np.sqrt(3.)*J + np.sqrt(3.*J*J - 1.) )
        return F

    def compute_Fprime(self,J):
        F = self.compute_F(J)
        Fp= (np.sqrt(3.) + 3.*J*(3.*J*J -1.)**(-0.5) )/np.exp(F)
        return Fp

    def compute_G(self,J,Jside):
        G = np.sqrt( (Jside-J)*(Jside**3-Jside-(J**3-J)) )
        return G

    def compute_Gprime(self,J,Jside):
        G = self.compute_G(J,Jside)
        Gp = ((4.*J-3*Jside)*J**2 + 2.*(Jside-J) - Jside**3)/(2.*G)
        return Gp

    def compute_H(self,J):
        H = np.sqrt( 3.*J**2 -1.)
        return H

    def compute_Hprime(self,J):
        H = self.compute_H(J)
        Hp = 3.*J/H
        return Hp

    def compute_A(self,rho,E,vL,vR,J):
        A=2.*np.sqrt(2./(E*rho))*(vR-vL)  - self.compute_F(J)/np.sqrt(3.) + J*self.compute_H(J)
        return A

    def computeResidual(self,rho,E,vL,vR,J,J_A,Ji):
        Res = Ji-(self.compute_F(Ji)/np.sqrt(3.) - 2.*self.compute_G(Ji,J) + self.compute_A(rho,E,vL,vR,J_A) )/self.compute_H(Ji)
        return Res

    def compute_dRes(self,rho,E,vL,vR,J,J_A,Ji):
        f = self.compute_F(Ji) ; df = self.compute_Fprime(Ji)
        g = self.compute_G(Ji,J) ; dg = self.compute_Gprime(Ji,J)
        h = self.compute_H(Ji) ; dh = self.compute_Hprime(Ji)
        A = self.compute_A(rho,E,vL,vR,J_A)
        dR = 1. - (df/np.sqrt(3.) -2.*dg)/h + (f/np.sqrt(3.) -2.*g +A)*dh/(h**2)
        return dR


    def computeQ2_star_1R2S(self,E,rho,v,J,Ji):
        q2 = v + np.sqrt(E*rho/24.)*(np.sqrt(3.)*(Ji*self.compute_H(Ji) - J*self.compute_H(J) ) - self.compute_F(Ji) + self.compute_F(J))
        return q2

    def computeQ2_star_1S2R(self,E,rho,v,J,Ji):
        q2 = v - np.sqrt(E*rho/8.)*(Ji*self.compute_H(Ji) - J*self.compute_H(J)  - (self.compute_F(Ji) - self.compute_F(J))/np.sqrt(3.))
        return q2

    def computeQ2_star_1S2R_2(self,E,rho,v,J,Ji):
        q2 = v + np.sqrt(E*rho/2.)*self.compute_G(Ji,J)
        return q2


   ###################################################################################
   ###################################################################################
   ## METHODS FOR LIMITERS

    def minmod(self,*args):
        s=np.zeros(2)
        for i in range(len(args)):
            s+=np.sign(args[i])
        s/=len(args)
        if (np.abs(s)==1).all() :
            return s*np.min(np.abs(args))
        else :
            return 0.

    def maxmod(self,*args):
        s=np.zeros(2)
        for i in range(len(args)):
            s+=np.sign(args[i])
        s/=len(args)
        if (np.abs(s)==1).all() :
            return s*np.max(np.abs(args))
        else :
            return 0.


    def slopeLimiter(self,U,limiter):
        dx = self.dx
        Nnodes = np.shape(U)[0]
        Nelem = Nnodes/2 -1
        for i in range(Nelem):
            # Build mean value vector
            Umean=0.5*(U[2*i+1,:]+U[2*i+2,:])
            if i==0 : UmeanL= U[0,:]
            else : UmeanL=0.5*(U[2*i-1,:]+U[2*i,:])
            if i==Nelem-1 : UmeanR= U[-1,:]
            else : UmeanR=0.5*(U[2*i+3,:]+U[2*i+4,:])
            dU = (U[2*i+2,:]-U[2*i+1,:])/dx
            dUL = (Umean-UmeanL)/dx
            dUR = (UmeanR-Umean)/dx
            if limiter==0 : s = self.minmodLimiter(dU,dUL,dUR)
            elif limiter==1 : s = self.superbeeLimiter(dU,dUL,dUR)
            elif limiter==2 : s = self.musclLimiter(dU,dUL,dUR)
            ul = Umean - 0.5*dx*s
            ur = Umean + 0.5*dx*s
            if not (np.abs(ul-U[2*i+1,:]).any()<1.e-8 and np.abs(ur-U[2*i+2,:]).any()<1.e-8):
                U[2*i+1,:]=ul
                U[2*i+2,:]=ur
        return U

    
    def minmodLimiter(self,dU,dUL,dUR):
        return self.minmod(dU,2.*dUR,2.*dUL)
    
    def superbeeLimiter(self,dU,dUL,dUR):
        sl = self.minmod(dUR,2.*dUL,dU)
        sr = self.minmod(2.*dUR,dUL,dU)
        si = self.minmod(dUR,dUL,2.*dU)
        return self.maxmod(sl,si,sr)
    
    def musclLimiter(self,dU,dUL,dUR):
        return self.minmod(dU,dUR,dUL)
