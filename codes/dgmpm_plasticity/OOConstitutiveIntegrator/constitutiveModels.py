# !/usr/bin/python
import scipy.optimize
import numpy as np
import pdb
from hardeningModels import *
from state import *
from kinematics import *
from material import *

class elasticity:

    def __init__(self,material,kinematic):
        self.material=material
        self.kin = kinematic
        self.fourthOrderIdentity()
        self.hydrostaticProjector()
        self.deviatoricProjector()
        
    def fourthOrderIdentity(self):
        self.I=np.eye(9)
        #self.I=np.asmatrix(I)

    def hydrostaticProjector(self):
        PH=np.zeros((9,9))
        #PH=np.array(PH)
        PH[0,:]=[1,0,0,0,1,0,0,0,1] ; PH[4,:]=[1,0,0,0,1,0,0,0,1] ; PH[8,:]=[1,0,0,0,1,0,0,0,1]
        self.PH=PH*(1./3.)
        
    def deviatoricProjector(self):
        self.PD = self.I - self.PH
        
    def computeCauchyStress(self,eps,sig):
        epsl = np.asmatrix(eps).A1 ; sigma = np.asmatrix(sig).A1
        ## epsl = [e11 e12 e13 e21 e22 e23 e31 e32 e33]
        result = 2.*self.material.mu*np.dot(self.PD,epsl) + 3.*self.material.kappa*np.dot(self.PH,epsl)
        sigma[self.kin.fluxMap]=result[self.kin.fluxMap]
        
class plasticity(elasticity):
    def __init__(self,material,kinematic,hardeningModel):
        self.material=material
        self.kin = kinematic
        self.fourthOrderIdentity()
        self.hydrostaticProjector()
        self.deviatoricProjector()
        self.hardeningModel = hardeningModel
        
class J2Plasticity(plasticity):

    def constitutiveUpdate(self,state0,state1):
        ## Strain driven problem
        self.computeCauchyStress(state1.grad-state0.get('PLASTIC_STRAIN'),state0.flux)
        s_trial,f = self.yieldSurface(state0)
        if (f<=0.):
            state1.internal = state0.internal
            state1.flux = state0.flux
            return 0.
        else :
            ## Radial return mapping
            X=np.asmatrix(self.hardeningModel.getCenter(state0)).A1
            s_trial_eq=np.sqrt((3./2.)*np.dot(s_trial-X,s_trial-X))
            Res=lambda x: x - f/(3.*self.material.mu + self.hardeningModel.slope(state0.get("CUMULATED_PLASTIC_STRAIN")+x))
            solution=scipy.optimize.fsolve(Res,0.,full_output=1)
            DeltaP=solution[0][0]
            assert DeltaP>=0.,"Plastic multiplier lower than 0!!"
            if not solution[2]:
                print "Problem in radial return mapping!"
                pdb.set_trace()
            flow = np.sqrt(3./2.)*s_trial/s_trial_eq
            epsp = np.asmatrix(state0.get("PLASTIC_STRAIN")).A1
            epsp[self.kin.fluxMap] += np.sqrt(3./2.)*DeltaP*flow[self.kin.fluxMap]
            state1.internal[0]=state0.internal[0]
            epsl = np.asmatrix(state1.grad).A1-epsp
            self.computeCauchyStress(state1.grad-state1.get('PLASTIC_STRAIN'),state1.flux)
            state1.setState(state1.grad,state1.flux,[state1.get("PLASTIC_STRAIN"),state0.get('CUMULATED_PLASTIC_STRAIN')+DeltaP],self.hardeningModel)
            return DeltaP
        
    
    def yieldSurface(self,state):
        sig=np.asmatrix(state.flux).A1 ; 
        s_trial=np.dot(self.PD,sig)
        X=np.asmatrix(self.hardeningModel.getCenter(state)).A1
        sEq=np.sqrt(3./2.)*np.sqrt(np.dot(s_trial-X,s_trial-X))
        return s_trial,sEq-self.hardeningModel.getRadius(state.get("CUMULATED_PLASTIC_STRAIN"))
    


## State (needs hardening->nIntVar() )
## Hardening (gives internal variables, needs Kinematics to determine the number of components)
## Yield surface (needs hardening->getCenter() and hardening->getRadius() )
## 
## Tangent modulus ; Kinematics (python matrix with "flat") ; material (parameters)
## elasticity ; plasticity
