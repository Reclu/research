# !/usr/bin/python
import numpy as np
import pdb

"""
Based on the python flat storage of numpy.matrix, which ascribe to the component of a matrix the following indices:
|0 1 2|
|3 4 5|
|6 7 8|,
this class provides the index of active components for a given kinematic (i.e. bar, plane strain etc.).
Then, all the tensor quantities are written as 3x3 but only the active components are changed.
"""

class kinematic:
    ## General 3D problem
    def __init__(self):
        self.vectorMap = [0,1,2]
        self.fluxMap = [0,1,2,3,4,5,6,7,8]
        self.gradMap = [0,1,2,3,4,5,6,7,8]

class barKinematic:
    def __init__(self):
        self.vectorMap = [0]
        self.fluxMap = [0]
        self.gradMap = [0]

class planeWaveKinematic:
    def __init__(self):
        self.vectorMap = [0]
        self.fluxMap = [0,4,8]
        self.gradMap = [0]

class planeStrainKinematic:
    def __init__(self):
        self.vectorMap = [0,1]
        self.fluxMap = [0,1,3,4,8]
        self.gradMap = [0,1,3,4]

class planeStressKinematic:
    def __init__(self):
        self.vectorMap = [0,1]
        self.fluxMap = [0,1,3,4]
        self.gradMap = [0,1,3,4,8]

