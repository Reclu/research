# !/usr/bin/python
import numpy as np
import pdb

class elasticMaterial:
    def __init__(self,Young,Poisson,Density):
        self.isPlastic=False
        self.E = Young
        self.nu = Poisson
        self.rho = Density
        self.kappa = Young/(3.*(1.-2.*nu))
        self.mu = Young/(2.*(1.+nu))
        self.lamb = self.kappa - 2.*self.mu/3.
        
class plasticMaterial:
    def __init__(self,Young,Poisson,Density,yieldStress):
        self.isPlastic=True
        self.E = Young
        self.nu = Poisson
        self.rho = Density
        self.kappa = Young/(3.*(1.-2.*Poisson))
        self.mu = Young/(2.*(1.+Poisson))
        self.lamb = self.kappa - 2.*self.mu/3.
        self.sigy = yieldStress
