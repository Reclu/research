# !/usr/bin/python
import numpy as np
import pdb
from state import *
from material import *

class baseHardening:
    def __init__(self,material,kinematic_hardening_modulus=None,isotropic_hardening_modulus=None,isotropic_power_exponent=None):
        assert material.isPlastic, "Plastic material must be given"
        self.mat=material
        self.K = kinematic_hardening_modulus
        self.H = isotropic_hardening_modulus
        self.IH_power = isotropic_power_exponent
        self.IntVarNames = ['']
        self.IntVarSizes = []
        
    def getIntVarNames(self):
        return self.IntVarNames

    def getIntSizes(self):
        return self.IntVarSizes
    
    def updateInternalVariable(self,dp):
        ## Better to di this here?
        return 0.

    def getCenter(self,state):
        return state.get("PLASTIC_STRAIN")*0.#(2./3.)*self.K*state.internal[0]

    def getRadius(self,p):
        return self.mat.sigy


class isotropicHardening(baseHardening):
    def __init__(self,material,isotropic_hardening_modulus,isotropic_power_exponent=1.):
        assert material.isPlastic, "Plastic material must be given"
        self.mat=material
        #self.K = kinematic_hardening_modulus
        self.H = isotropic_hardening_modulus
        self.IH_power = isotropic_power_exponent
        self.IntVarNames = ['PLASTIC_STRAIN','CUMULATED_PLASTIC_STRAIN']
        self.IntVarSizes = [9,1]

    def getRadius(self,p):
        return self.mat.sigy + self.H*(p)**self.IH_power

    def slope(self,p):
        return self.IH_power*self.H*(p)**(self.IH_power-1.)

class kinematicHardening(baseHardening):
    ## Linear model
    def __init__(self,material,kinematic_hardening_modulus):
        assert material.isPlastic, "Plastic material must be given"
        self.mat=material
        self.K = kinematic_hardening_modulus
        self.IntVarNames = ['PLASTIC_STRAIN','CUMULATED_PLASTIC_STRAIN']
        self.IntVarSizes = [9,1]

    def getCenter(self,state):
        return (2./3.)*self.K*state.get("PLASTIC_STRAIN")

    def slope(self,p):
        return 3./2.*self.K
