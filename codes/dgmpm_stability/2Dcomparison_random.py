#!/usr/bin/python

import numpy as np
from scipy import optimize
from sympy import *
import matplotlib.pyplot as plt
import random
import pdb
import os

def export2DTeXFile(fileName,xField,fields,*kwargs):
    TeXFile=open(fileName,"w")
    n_fields = np.shape(fields)[0]
    n_labels = np.shape(kwargs)[0]
    # Define Paul Tol's colors (purple to red)
    color=['Blue','Red','Green','Red','black','black','black']
    marker=['+','x','star','+','none','none','none']
    size=['very thick','very thick','very thick','very thick','thin','thin',]
    line=['solid','solid','dashed','dashed']
    TeXFile.write(r'\begin{tikzpicture}[scale=0.5]')
    TeXFile.write('\n')
    TeXFile.write(r'\begin{axis}[xlabel=$s_1/s_2$,ymajorgrids=true,xmajorgrids=true,xmin=1,xmax=41,xtick={1,10,20,30,40}]')
    TeXFile.write('\n')
    TeXFile.write('%%%%%%%%%%% NATURAL CONFIGURATION')
    TeXFile.write('\n')
    #pdb.set_trace()
    for i in range(np.shape(fields)[0]):
        TeXFile.write(r'\addplot['+str(color[i])+',mark='+str(marker[i])+',very thick,mark size=5pt] coordinates {')
        for j in range(len(fields[i,:])):
            TeXFile.write('('+str(xField[j])+','+str(fields[i,j])+') ')
        
        TeXFile.write('};\n')
        if i==0:
            TeXFile.write('%%%%%%%%%%% MODIFIED CONFIGURATION')
            TeXFile.write('\n')

    TeXFile.write(r'\end{axis}')
    TeXFile.write('\n')
    TeXFile.write('\end{tikzpicture}')
    TeXFile.write('\n')
    TeXFile.write('%%% Local Variables:')
    TeXFile.write('\n')
    TeXFile.write('%%% mode: latex')
    TeXFile.write('\n')
    TeXFile.write('%%% TeX-master: "../../mainManuscript"')
    TeXFile.write('\n')
    TeXFile.write('%%% End:')
    TeXFile.write('\n')
    TeXFile.close()


# Symbolic function to evaluate shape functions
shape_functions=lambda x,y: np.array([(1.-x)*(1.-y)/4.,(1.+x)*(1.-y)/4.,(1.+x)*(1.+y)/4.,(1.-x)*(1.+y)/4.])
grad_xi=lambda y:np.array([-(1.-y)/4.,(1.-y)/4.,(1.+y)/4.,-(1.+y)/4.])
grad_eta=lambda x:np.array([-(1.-x)/4.,-(1.+x)/4.,(1.+x)/4.,(1.-x)/4.])

# shapes=| N1(Xp1) N1(Xp2) ... N1(XNp) |
#        | N2(Xp1) N2(Xp2) ... N2(XNp) |
#        | N3(Xp1) N3(Xp2) ... N3(XNp) |
#        | N4(Xp1) N4(Xp2) ... N4(XNp) |

# grad_z=| N1_z(Xp1) N1_z(Xp2) ... N1_z(XNp) |
#        | N2_z(Xp1) N2_z(Xp2) ... N2_z(XNp) |
#        | N3_z(Xp1) N3_z(Xp2) ... N3_z(XNp) |
#        | N4_z(Xp1) N4_z(Xp2) ... N4_z(XNp) |

# where Ni(Xj) is the shape function of node i evaluated at the jth particles position


def symbolResidual(point,dx,cx,cy,XC,XB,XL,XBL=0):
    transverse=True
    if XBL==0: transverse=False
    shapesC=shape_functions(XC[0],XC[1])
    dSxi_C=grad_xi(XC[1])
    dSeta_C=grad_eta(XC[0])
    shapesB=shape_functions(XB[0],XB[1])
    dSxi_B=grad_xi(XB[1])
    dSeta_B=grad_eta(XB[0])
    shapesL=shape_functions(XL[0],XL[1])
    dSxi_L=grad_xi(XL[1])
    dSeta_L=grad_eta(XL[0])
    ## Number of material points in cells
    NmpC=len(XC[0])
    NmpL=len(XL[0])
    NmpB=len(XB[0])
    if XBL!=0:
        shapesBL=shape_functions(XBL[0],XBL[1])
        dSxi_BL=grad_xi(XBL[1])
        dSeta_BL=grad_eta(XBL[0])
        NmpBL=len(XBL[0])
    else:
        NmpBL=0

    dt = symbols('dt')
    ## sum_i^K = np.sum(shapesK[i,:]) with cell K and node i
    ## shape functions evaluated at edges centers to weight fluxes contributions
    ## o -- 3 -- o
    ## |         |
    ## 4         2
    ## |         |
    ## o -- 1 -- o
    shapeOnEdge=shape_functions(np.array([0.,1.,0.,-1.]),np.array([-1.,0.,1.,0.]))
    ## Define the normal to edges
    Nx=np.array([0.,1.,0.,-1.])
    Ny=np.array([-1.,0.,1.,0.])
    Nnodes=4
    Nedges=4
    Res=0.
    for P in range(NmpC):
        ## Contributions of material points sharing the same cell
        D_PI=0.
        for i in range(Nnodes):
            # 0th-order contributions
            wheightings=shapesC[i,point]/np.sum(shapesC[i,:])
            D_PI+=wheightings*shapesC[i,P]
            # 1st-order contributions
            for j in range(Nnodes):
                D_PI+=2.*dt*wheightings*(shapesC[j,P]/np.sum(shapesC[j,:]))*(cx*np.dot(dSxi_C[i,:],shapesC[j,:])/dx + cy*np.dot(dSeta_C[i,:],shapesC[j,:])/dx)
            # Contributions of edges 2 and 3
            #pdb.set_trace()
            D_PI-=0.5*(dt/dx)*wheightings*shapeOnEdge[i,1]*NmpC*cx*(shapesC[1,P]/np.sum(shapesC[1,:])+shapesC[2,P]/np.sum(shapesC[2,:]))
            D_PI-=0.5*(dt/dx)*wheightings*shapeOnEdge[i,2]*NmpC*cy*(shapesC[2,P]/np.sum(shapesC[2,:])+shapesC[3,P]/np.sum(shapesC[3,:]))
            # Transverse contributions
            if transverse:
                D_PI+= 0.25*(dt/dx)**2*wheightings*shapeOnEdge[i,1]*NmpC*cx*cy*(shapesC[0,P]/np.sum(shapesC[0,:])+shapesC[1,P]/np.sum(shapesC[1,:]))
                D_PI+= 0.25*(dt/dx)**2*wheightings*shapeOnEdge[i,2]*NmpC*cx*cy*(shapesC[0,P]/np.sum(shapesC[0,:])+shapesC[3,P]/np.sum(shapesC[3,:]))
        Res+=np.abs(D_PI)
    ## Contributions of material points of left cell
    for P in range(NmpL):
        D_PI=0.
        for i in range(Nnodes):
            wheightings=shapesC[i,point]/np.sum(shapesC[i,:])
            ## edge 4 contribution
            D_PI+= 0.5*(dt/dx)*wheightings*shapeOnEdge[i,3]*NmpC*cx*(shapesL[1,P]/np.sum(shapesL[1,:])+shapesL[2,P]/np.sum(shapesL[2,:]))
            if transverse:
                D_PI-=0.25*(dt/dx)**2*wheightings*shapeOnEdge[i,3]*NmpC*cx*cy*(shapesL[0,P]/np.sum(shapesL[0,:])+shapesL[1,P]/np.sum(shapesL[1,:]))
                ## edge 3 contribution
                D_PI-=0.25*(dt/dx)**2*wheightings*shapeOnEdge[i,2]*NmpC*cy*cx*(shapesL[1,P]/np.sum(shapesL[1,:])+shapesL[2,P]/np.sum(shapesL[2,:]))
        Res+=np.abs(D_PI)
    ## Contributions of material points of bottom cell            
    for P in range(NmpB):
        D_PI=0.
        for i in range(Nnodes):
            wheightings=shapesC[i,point]/np.sum(shapesC[i,:])
            ## edge 1 contribution
            D_PI+= 0.5*(dt/dx)*wheightings*shapeOnEdge[i,0]*NmpC*cy*(shapesB[2,P]/np.sum(shapesB[2,:])+shapesB[3,P]/np.sum(shapesB[3,:]))
            if transverse:
                D_PI-=0.25*(dt/dx)**2*wheightings*shapeOnEdge[i,0]*NmpC*cy*cx*(shapesB[0,P]/np.sum(shapesB[0,:])+shapesB[3,P]/np.sum(shapesB[3,:]))
                ## edge 2 contribution
                D_PI-=0.25*(dt/dx)**2*wheightings*shapeOnEdge[i,1]*NmpC*cx*cy*(shapesB[2,P]/np.sum(shapesB[2,:])+shapesB[3,P]/np.sum(shapesB[3,:]))
        Res+=np.abs(D_PI)
    ## Contributions of material points of bottom-left cell
    for P in range(NmpBL):
        D_PI=0.
        for i in range(Nnodes):
            wheightings=shapesC[i,point]/np.sum(shapesC[i,:])
            ## edge 1 contribution
            D_PI+=0.25*(dt/dx)**2*wheightings*shapeOnEdge[i,0]*NmpC*cy*cx*(shapesBL[1,P]/np.sum(shapesBL[1,:])+shapesBL[2,P]/np.sum(shapesBL[2,:]))
            ## edge 4 contribution
            D_PI+=0.25*(dt/dx)**2*wheightings*shapeOnEdge[i,3]*NmpC*cx*cy*(shapesBL[2,P]/np.sum(shapesBL[2,:])+shapesBL[3,P]/np.sum(shapesBL[3,:]))
        Res+=np.abs(D_PI)
    Residual = lambdify((dt),Res-1.)
    return Residual

def gridSearch(function,dx,cx,tol=1.e-2):
    samples=10000
    # Find the bigest root of the residual by grid search algorithm
    CFL=np.linspace(0.,1.,samples)
    for i in range(samples):
        value=CFL[-1-i]
        a0=function(value*dx/cx)
        if a0<tol:
            return value
        else:
            continue
    return 0.

def Rand(start, end, num): 
    res = [] 
  
    for j in range(num): 
        res.append(random.randint(start, end)) 
  
    return np.asarray(res)

def RandPosition(numberOfPoints):
    res=[]
    for nPoints in(numberOfPoints):
        position=np.zeros((nPoints,2))
        for i in range(nPoints):
            position[i,0]=random.uniform(-1., 1.)
            position[i,1]=random.uniform(-1., 1.)
        res.append(position)
    return res


# samples=20
# cx=np.linspace(2.,80.,samples)
# cy=cx[0]
cx=2.
cy=2.
dx=2.

samples=1000
number_left = Rand(1, 4, samples)
position_left = RandPosition(number_left)

number_bott = Rand(1, 4, samples)
position_bott = RandPosition(number_bott)

number_curr = Rand(1, 4, samples)
position_curr = RandPosition(number_curr)

number_botle = Rand(1, 4, samples)
position_botle = RandPosition(number_botle)

if not os.path.exists('dcuRandom.npy'):
    dcuSolution=[]
    dcuSolution_id=[]
    ctuSolution=[]
    ctuSolution_id=[]
    for i in range(samples):
        print "Computing critical CFL for sample ",i,": ",number_curr[i]," particles"
        solution_dcu=[]
        solution_dcu_id=[]
        solution_ctu=[]
        solution_ctu_id=[]
        for k in range(number_curr[i]):
            # if number_curr[i]<number_prev[i] :
            #     print "Attention ca va merder !!!!!!"
            # else:
            #     print "Ca va le faire..."
            XL = position_left[i][:,0] ; YL = position_left[i][:,1]
            XB = position_bott[i][:,0] ; YB = position_bott[i][:,1]
            XBL = position_botle[i][:,0] ; YBL = position_botle[i][:,1]
            XC = position_curr[i][:,0] ; YC = position_curr[i][:,1]
            
            res=symbolResidual(k,dx,cx,cy,(XC,YC),(XB,YB),(XL,YL))
            solution_dcu.append(gridSearch(res,dx,cx))
            res=symbolResidual(k,dx,cx,cy,(XC,YC),(XC,YC),(XC,YC))
            solution_dcu_id.append(gridSearch(res,dx,cx))

            res=symbolResidual(k,dx,cx,cy,(XC,YC),(XB,YB),(XL,YL),(XBL,YBL))
            solution_ctu.append(gridSearch(res,dx,cx))
            res=symbolResidual(k,dx,cx,cy,(XC,YC),(XC,YC),(XC,YC),(XC,YC))
            solution_ctu_id.append(gridSearch(res,dx,cx))
            
        dcuSolution.append(min(solution_dcu))
        dcuSolution_id.append(min(solution_dcu_id))
        ctuSolution.append(min(solution_ctu))
        ctuSolution_id.append(min(solution_ctu_id))
        
    np.save('dcuRandom.npy',dcuSolution)
    np.save('dcuRandom_id.npy',dcuSolution_id)
    np.save('ctuRandom.npy',ctuSolution)
    np.save('ctuRandom_id.npy',ctuSolution_id)
else :
    dcuSolution=np.load('dcuRandom.npy')
    dcuSolution_id=np.load('dcuRandom_id.npy')
    ctuSolution=np.load('ctuRandom.npy')
    ctuSolution_id=np.load('ctuRandom_id.npy')

import statistics
plt.figure()
plt.hist(dcuSolution,bins='auto',color='blue')
plt.grid()

plt.figure()
plt.hist(dcuSolution_id,bins='auto',color='red')
plt.grid()

plt.show()
pdb.set_trace()

plt.figure()
plt.hist(ctuSolution,bins='auto',color='blue')
plt.grid()

plt.figure()
plt.hist(ctuSolution_id,bins='auto',color='red')
plt.grid()

plt.show()
