# !/usr/bin/python

import numpy as np
import pdb

class DGmesh:
    #Constructor
    def __init__(self,Mp,l,ppc,c,cp,rho,sigy,H,hardening):
        nex = Mp/ppc
        nnx=2*nex+2
        lmp=l/(Mp-1)
        self.dx = ppc*l/nex
        self.xn = np.linspace(-lmp/2,l+lmp/2,nex+1)
        self.connect = np.array([np.arange(1,nnx-1,2),np.arange(2,nnx,2)]).T
        self.c = c
        self.cp = cp
        self.rho = rho
        self.sigy=sigy
        self.H = H
        self.hard = hardening
        
    def set_sigy(self,sigy):
        self.sigy=sigy

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


    def computeFlux(self,U,W,dofs,md,epeq):
        """
        Evaluation of interface fluxes by solving Riemann problems
        Computation of volumic fluxes over each cell
        Assembling of interfaces and volumes contributions
        """
        fluxes = self.computeInterfaceFlux(W,epeq)
        fint = self.computeInternalForces(W,dofs)
        f = self.computeTotalForces(fint,fluxes)
        return f

    def computeDelta(self,dU,cL,cR):
        delta1 = (dU[0]+self.rho*cR*dU[1])/(self.rho*cR+self.rho*cL)
        delta2 = -(dU[0]-self.rho*cL*dU[1])/(self.rho*cR+self.rho*cL)
        return [delta1,delta2]

    def computeElasticWaves(self,delta):
        waves = np.zeros((2,2))
        waves[:,0] = delta[0]*np.array([self.rho*self.c,1.])
        waves[:,1] = delta[1]*np.array([-self.rho*self.c,1.])
        return waves

    def computeInterfaceFlux(self,W,EPeq):
        debug = False
        Nnodes = np.shape(W)[0]
        Nelem = (Nnodes-2)/2
        fluxes = np.zeros((Nnodes,W.shape[1]))
        flux = np.zeros((Nnodes,W.shape[1]))
        #Loop on interfaces
        for i in range(Nelem+1):
            WL = W[2*i,:] ; SL = WL[0]
            WR = W[2*i+1,:] ; SR = WR[0]
            dW = WR-WL
            #Elastic prediction
            delta = self.computeDelta(dW,self.c,self.c)
            Strial= SL + delta[0]*self.c*self.rho
            #Tests on the criterion
            if self.hard == 'isotropic':
                yieldL=self.H*EPeq[2*i]+self.sigy
                yieldR=self.H*EPeq[2*i+1]+self.sigy
                StrialL= Strial
                StrialR= Strial
            elif self.hard == 'kinematic':
                yieldL=self.sigy
                yieldR=self.sigy
                StrialL= Strial - self.H*EPeq[2*i]
                StrialR= Strial - self.H*EPeq[2*i+1]
                print "Carefull, kinematic hardening selected although it is not correctly coded here !!"
                pdb.set_trace()
            fL = np.abs(StrialL)- yieldL
            fR = np.abs(StrialR)- yieldR
            
            if debug: print np.abs(Strial),self.sigy,fL,fR,
            if (fL>0.0) and (fR<0.0):
                if debug : print i," plastic-elastic"
                #Leftward plastic wave
                delta1=(np.sign(StrialL)*(yieldL)-SL)/(self.rho*self.c)
                R = WR- WL + delta1*np.array([self.rho*self.c,1.0])
                deltaP = self.computeDelta(R,self.cp,self.c)
                Wstar = WL + delta1*np.array([self.rho*self.c,1.])+ deltaP[0]*np.array([self.rho*self.cp,1.])
            elif (fL<0.0) and (fR>0.0):
                if debug : print i," elastic-plastic"
                #Rightward plastic wave
                delta2=(np.sign(StrialR)*(yieldR)-SR)/(self.rho*self.c)
                R = WR-delta2*np.array([-self.rho*self.c,1.0]) -WL
                deltaP = self.computeDelta(R,self.c,self.cp)
                Wstar = WR - delta2*np.array([-self.rho*self.c,1.])- deltaP[1]*np.array([-self.rho*self.cp,1.])
            elif (fL>0.0) and (fR>0.0):
                if debug : print i," plastic-plastic"
                #Four waves
                delta1=(np.sign(StrialL)*(yieldL)-SL)/(self.rho*self.c)
                delta2=(np.sign(StrialR)*(yieldR)-SR)/(self.rho*self.c)
                R = WR - delta2*np.array([-self.rho*self.c,1.0]) - ( WL + delta1*np.array([self.rho*self.c,1.0]))
                deltaP = self.computeDelta(R,self.cp,self.cp)
                Wstar = WR - delta2*np.array([-self.rho*self.c,1.])- deltaP[1]*np.array([-self.rho*self.cp,1.])
            else:
                if debug : print i," elastic-elastic"
                # Elastic fluctuations
                waves = self.computeElasticWaves(delta)
                Wstar = WL + waves[:,0]
            fluxes[2*i,:] = np.array([-Wstar[1],-Wstar[0]])
            fluxes[2*i+1,:] = np.array([-Wstar[1],-Wstar[0]])
        return fluxes

    def computeTotalForces(self,fint,fluxes):
        Nnodes = np.shape(fint)[0]
        Nelem = (Nnodes-2)/2
        # Loop on elements
        
        f=np.copy(fint) # Weak
        for i in range(Nelem):
            mapp = self.connect[i,:]
            f[int(mapp[0]),:]+=fluxes[int(mapp[0]),:] # Left node
            f[int(mapp[1]),:]-=fluxes[int(mapp[1]),:] # Right node
        return f
    
    def computeFluxVector(self,W):
        Nnodes = np.shape(W)[0]
        flux = np.zeros((Nnodes,2))
        for i in range(Nnodes)[1:-1]:
            flux[i,:] = np.array([-W[i,1]/self.rho,-W[i,0]/self.rho])
        return flux

    def computeInternalForces(self,W,dofs):
        S = self.Kx
        Nnodes = np.shape(W)[0]
        fint = np.zeros((Nnodes,2))
        flux = self.computeFluxVector(W)
        fint[dofs,0] = np.dot(S,flux[dofs,0])
        fint[dofs,1] = np.dot(S,flux[dofs,1])
        return fint

    def setMapping(self,Kx):
        self.Kx = Kx
        
    

