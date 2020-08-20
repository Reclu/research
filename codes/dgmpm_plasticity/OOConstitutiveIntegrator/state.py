# !/usr/bin/python
import numpy as np
import pdb
from hardeningModels import *

class state:
    def __init__(self,Grad,Flux,Internal,hardeningModel):
        assert Grad.shape==(3,3), "Epsilon dimension should (3,3)!"
        assert Flux.shape==(3,3), "Sigma dimension should (3,3)!"

        # self.grad = np.copy(np.asmatrix(Grad))
        # self.flux = np.copy(np.asmatrix(Flux))
        self.grad = np.copy(Grad)
        self.flux = np.copy(Flux)

        
        self.internal=[]
        self.IntSizeList=hardeningModel.getIntSizes()
        self.IntNameList=hardeningModel.getIntVarNames()
        assert len(Internal)==len(self.IntSizeList), "The number of internal variables given is not correct!"
        for i in range(len(Internal)):
            if  (type(Internal[i]) is np.float64) or (type(Internal[i]) is float) :
                assert self.IntSizeList[i]==1, "Scalar value given for an internal variable that is not scalar!"
                self.internal.append(Internal[i])
            else :
                assert Internal[i].shape==(3,3), "Internal variable "+str( hardeningModel.getIntVarNames()[i])+" dimension should (3,3)!"
                self.internal.append(np.asmatrix(Internal[i]))

    def get(self,key):
        if key=="STRAIN":
            return self.grad
        elif key=="STRESS":
            return self.flux

        elif key in self.IntNameList :
            return self.internal[self.IntNameList.index(key)]
        else:
            assert 0., "Unkown key to get state variables!!"

    def setState(self,Grad,Flux,Internal,hardeningModel):
        assert Grad.shape==(3,3), "Epsilon dimension should (3,3)!"
        assert Flux.shape==(3,3), "Sigma dimension should (3,3)!"

        # self.grad = np.copy(np.asmatrix(Grad))
        # self.flux = np.copy(np.asmatrix(Flux))
        self.grad = np.copy(Grad)
        self.flux = np.copy(Flux)

        
        self.internal=[]
        self.IntSizeList=hardeningModel.getIntSizes()
        self.IntNameList=hardeningModel.getIntVarNames()
        assert len(Internal)==len(self.IntSizeList), "The number of internal variables given is not correct!"
        for i in range(len(Internal)):
            if  (type(Internal[i]) is np.float64) or (type(Internal[i]) is float) :
                assert self.IntSizeList[i]==1, "Scalar value given for an internal variable that is not scalar!"
                self.internal.append(Internal[i])
            else :
                assert Internal[i].shape==(3,3), "Internal variable "+str( hardeningModel.getIntVarNames()[i])+" dimension should (3,3)!"
                self.internal.append(np.asmatrix(Internal[i]))
