#!/usr/bin/python
import methodes as m
import numpy as np

### Callable object to compute the diretized flux vector associated to a field U given as argument
class Flux:
    """
    Callable stockant un compteur que nous pouvons incrementer a chaque
    appel en lui donnant comme argument un entier.
    """

    def __init__(self,dim,E,rho,c,x):
        """
        Initialisation du compteur.
        """
        self.dim=dim
        self.E = E
        self.c = c
        self.rho=rho
        self.x=x
        
    def __call__(self, U):
        """
        Methode executee si l'objet est appele comme une fonction.
        """
        if self.dim==2:
            amdq,apdq = m.computeInterfaceFlux(U,self.rho,self.c)
            fint = m.computeInternalForces(self.x,U,self.E,self.rho)
            f = m.computeTotalForces(fint,amdq,apdq)
            return f
        elif self.dim==1:
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
            xi_integ = np.array([0.774596669241483,0.0,-0.774596669241483])
            w_integ = np.array([0.555555555555556,0.888888888888889,0.555555555555556])
            f = np.zeros(Nnodes)
            for i in range(Nelem):
                mapp = T10[i,:]
                mapp_interfaces = mapp - i -1
                Xe = self.x[[i,i+1]]
                loc_jac = m.Jac(Xe)
                M = m.Lin_shp(2,xi_integ,Xe)
                D = np.array(w_integ)  # Gauss weights for Flux first component for all nodes
                flux = self.c*U[mapp] # Evaluate flux at gauss points
                D*=m.Lin_eval(xi_integ,flux) # Multiply flux values by weights   
                f[mapp] = (loc_jac*np.dot(M,D))
                ### Assemble RHS 
                f[mapp[0]]+=apdq[mapp_interfaces[0]] # left node
                f[mapp[1]]-=apdq[mapp_interfaces[1]] # right node
            return f
