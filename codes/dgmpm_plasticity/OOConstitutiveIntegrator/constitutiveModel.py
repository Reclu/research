# !/usr/bin/python
import scipy.optimize
import numpy as np
import pdb

class planeWaveState:
    #Constructor
    def __init__(self,Grad,Flux,Internal):
        #assert ((type(Grad) is float) or (type(Grad) is float (np.float64))), "Epsilon dimension should be one for plane waves!"
        assert len(Flux)==2, "Sigma dimension should be two for plane waves!"
        assert len(Internal)==2, "Internal variables should two one for plane waves!"
        self.grad = np.zeros((3,3))  
        self.flux = np.zeros((3,3))  
        self.internal = [np.zeros((3,3)),0.] 
        self.setState(Grad,Flux,Internal)
        
    def setState(self,grad,flux,internal):
        #pdb.set_trace()
        #assert ((type(grad) is float) or (type(grad) is float (np.float64))), "Epsilon dimension should be one for plane waves!"
        assert len(flux)==2, "Sigma dimension should be two for plane waves!"
        assert len(internal)==2, "Internal variables should be two for plane waves!"
        ## Strain components
        self.grad[0,0]=np.copy(grad)
        ## Stress components
        self.flux[0,0]=np.copy(flux[0])
        self.flux[1,1]=np.copy(flux[1])
        self.flux[2,2]=np.copy(flux[1])
        ## Plastic strain components
        self.internal[0][0,0]=np.copy(internal[0])
        self.internal[0][1,1]=-np.copy(internal[0])/2.
        self.internal[0][2,2]=-np.copy(internal[0])/2.
        ## Equivalent plastic strain
        self.internal[1]=internal[1]
        
class J2Plasticity:
    #Constructor for isotropic hardening
    def __init__(self,lamb,mu,sigy,H,n):
        self.lamb=lamb
        self.mu=mu
        self.sigy=sigy
        self.Hardening=H
        self.IH_power=n

    def deviatoricPart(self,sig):
        traceI=np.sum(np.diagonal(sig))*np.eye(3)
        return sig-traceI/3.

    def doubleContract(self,tens1,tens2):
        prod=np.dot(tens1,tens2.T)
        return np.sum(np.diagonal(prod))

    def trialStress(self,state0,state1):
        eps_elas=state1.grad-state0.internal[0]
        return 2.*self.mu*(eps_elas) + self.lamb*np.sum(np.diagonal(eps_elas))*np.eye(3)
    
    def yieldSurface(self,state0,state1):
        s_trial=self.deviatoricPart(self.trialStress(state0,state1))
        sEq=np.sqrt(3./2.)*np.sqrt(self.doubleContract(s_trial,s_trial))
        return s_trial,sEq-self.sigy-self.Hardening*state0.internal[1]**self.IH_power
    
    def radialReturn(self,f,s_trial,state0,state1):
        s_trial_eq=np.sqrt((3./2.)*self.doubleContract(s_trial,s_trial))
        #Res=lambda x: x - (s_trial_eq -self.sigy-self.Hardening*(state0.internal[1]+x)**self.IH_power)/(3.*self.mu)
        Res=lambda x: x - f/(3.*self.mu + self.IH_power*self.Hardening*(state0.internal[1]+x)**(self.IH_power-1.))
        solution=scipy.optimize.fsolve(Res,0.,full_output=1)
        DeltaP=solution[0][0]
        assert DeltaP>=0.,"Plastic multiplier lower than 0!!"
        p=state0.internal[1] + DeltaP
        flow = np.sqrt(3./2.)*s_trial/s_trial_eq
        epsp = state0.internal[0] + np.sqrt(3./2.)*DeltaP*flow
        epsl = state1.grad-epsp
        sig = 2.*self.mu*epsl + self.lamb*np.sum(np.diagonal(epsl))*np.eye(3)
        return p,epsp,sig,solution[2]


## State (needs hardening->nIntVar() )
## Hardening (gives internal variables, needs Kinematics to determine the number of components)
## Yield surface (needs hardening->getCenter() and hardening->getRadius() )
## 
## Tangent modulus ; Kinematics (python matrix with "flat") ; material (parameters)
## elasticity ; plasticity
