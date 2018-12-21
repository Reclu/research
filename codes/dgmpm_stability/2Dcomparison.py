#!/usr/bin/python

import numpy as np
from scipy import optimize
from sympy import *
import matplotlib.pyplot as plt
import pdb

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


def symbolResidual(point,cx,cy,XC,XB,XL,XBL=0):
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
        shapesBL=shape_functions(XL[0],XL[1])
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
                D_PI+=dt*wheightings*(shapesC[j,P]/np.sum(shapesC[j,:]))*(cx*np.dot(dSxi_C[i,:],shapesC[j,:]) + cy*np.dot(dSeta_C[i,:],shapesC[j,:]))
            # Contributions of edges 2 and 3
            D_PI-=0.25*dt*wheightings*shapeOnEdge[i,1]*NmpC*cx*(shapesC[1,P]/np.sum(shapesC[1,:])+shapesC[2,P]/np.sum(shapesC[2,:]))
            D_PI-=0.25*dt*wheightings*shapeOnEdge[i,2]*NmpC*cx*(shapesC[2,P]/np.sum(shapesC[2,:])+shapesC[3,P]/np.sum(shapesC[3,:]))
            # Transverse contributions
            if transverse:
                D_PI+= (0.25*dt)**2*wheightings*shapeOnEdge[i,1]*NmpC*cx*cy*(shapesC[0,P]/np.sum(shapesC[0,:])+shapesC[1,P]/np.sum(shapesC[1,:]))
                D_PI+= (0.25*dt)**2*wheightings*shapeOnEdge[i,2]*NmpC*cx*cy*(shapesC[0,P]/np.sum(shapesC[0,:])+shapesC[3,P]/np.sum(shapesC[3,:]))
        Res+=np.abs(D_PI)
    ## Contributions of material points of left cell
    for P in range(NmpL):
        D_PI=0.
        for i in range(Nnodes):
            wheightings=shapesC[i,point]/np.sum(shapesC[i,:])
            ## edge 4 contribution
            D_PI+= 0.25*dt*wheightings*shapeOnEdge[i,3]*NmpC*cx*(shapesL[2,P]/np.sum(shapesL[2,:])+shapesL[3,P]/np.sum(shapesL[3,:]))
            if transverse:
                D_PI-=(0.25*dt)**2*wheightings*shapeOnEdge[i,3]*NmpC*cx*cy*(shapesL[0,P]/np.sum(shapesL[0,:])+shapesL[1,P]/np.sum(shapesL[1,:]))
                ## edge 3 contribution
                D_PI-=(0.25*dt)**2*wheightings*shapeOnEdge[i,2]*NmpC*cy*cx*(shapesL[1,P]/np.sum(shapesL[1,:])+shapesL[2,P]/np.sum(shapesL[2,:]))
        Res+=np.abs(D_PI)
    ## Contributions of material points of bottom cell            
    for P in range(NmpB):
        D_PI=0.
        for i in range(Nnodes):
            wheightings=shapesC[i,point]/np.sum(shapesC[i,:])
            ## edge 1 contribution
            D_PI+= 0.25*dt*wheightings*shapeOnEdge[i,0]*NmpC*cy*(shapesB[2,P]/np.sum(shapesB[2,:])+shapesB[3,P]/np.sum(shapesB[3,:]))
            if transverse:
                D_PI-=(0.25*dt)**2*wheightings*shapeOnEdge[i,0]*NmpC*cy*cx*(shapesB[0,P]/np.sum(shapesB[0,:])+shapesB[3,P]/np.sum(shapesB[3,:]))
                ## edge 2 contribution
                D_PI-=(0.25*dt)**2*wheightings*shapeOnEdge[i,1]*NmpC*cx*cy*(shapesB[2,P]/np.sum(shapesB[2,:])+shapesB[3,P]/np.sum(shapesB[3,:]))
        Res+=np.abs(D_PI)
    ## Contributions of material points of bottom-left cell
    for P in range(NmpBL):
        D_PI=0.
        for i in range(Nnodes):
            wheightings=shapesC[i,point]/np.sum(shapesC[i,:])
            ## edge 1 contribution
            D_PI+=(0.25*dt)**2*wheightings*shapeOnEdge[i,0]*NmpC*cy*cx*(shapesBL[1,P]/np.sum(shapesBL[1,:])+shapesBL[2,P]/np.sum(shapesBL[2,:]))
            ## edge 4 contribution
            D_PI+=(0.25*dt)**2*wheightings*shapeOnEdge[i,3]*NmpC*cx*cy*(shapesBL[2,P]/np.sum(shapesBL[2,:])+shapesBL[3,P]/np.sum(shapesBL[3,:]))
        Res+=np.abs(D_PI)
    Residual = lambdify((dt),Res-1.)
    return Residual

def gridSearch(function,tol=1.e-7):
    samples=100000
    # Find the bigest root of the residual by grid search algorithm
    CFL=np.linspace(0.,1.,samples)
    for i in CFL:
        if i==samples-1:return i
        a0=function(i)
        if a0<tol:
            continue
        else:
            return i



samples=20
cx=np.linspace(2.,80.,samples)
cy=cx[0]
############### 1PPC
print "**************************************************************"
print "******************  1PPC discretization **********************"
print "**************************************************************"
print "=== Symmetric horizontal ==="
# Xp=np.array([-0.])
# Yp=np.array([0.])
# DCU=np.zeros(samples)
# CTU=np.zeros(samples)
# for i in range(samples):
#     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp))
#     solution=gridSearch(residual)
#     DCU[i]=cx[i]*solution/2.#.append(cx[i]*solution/2.)
#     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp),(Xp,Yp))
#     solution=optimize.newton(residual,1.)
#     #solution=gridSearch(residual)
#     CTU[i]=cx[i]*solution/2.#.append(cx[i]*solution/2.)
# plt.plot(cx/cy,DCU,'b')
# plt.plot(cx/cy,CTU,'r')
# plt.grid()
# plt.show()

############### 2PPC
print "**************************************************************"
print "******************  2PPC discretization **********************"
print "**************************************************************"
print "=== Symmetric horizontal ==="
# Xp=np.array([-0.5,0.5])
# Yp=np.array([0.,0.])
# Xp2=np.array([-0.25,0.25])
# Yp2=np.array([0.,0.])

# DCU_natural=np.zeros((2,samples))
# CTU_natural=np.zeros((2,samples))
# DCU_shifted=np.zeros((2,samples))
# CTU_shifted=np.zeros((2,samples))
# for i in range(samples):
#     for p in range(2):
#         residual=symbolResidual(p,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp))
#         solution=gridSearch(residual)
#         DCU_natural[p,i]=cx[i]*solution/2.#.append(cx[i]*solution/2.)
#         residual=symbolResidual(p,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp),(Xp,Yp))
#         solution=gridSearch(residual)
#         CTU_natural[p,i]=cx[i]*solution/2.#.append(cx[i]*solution/2.)

#         # Shift symmetrically
#         residual=symbolResidual(p,cx[i],cy,(Xp2,Yp2),(Xp2,Yp2),(Xp2,Yp2))
#         solution=gridSearch(residual)
#         DCU_shifted[p,i]=cx[i]*solution/2.#.append(cx[i]*solution/2.)
#         residual=symbolResidual(p,cx[i],cy,(Xp2,Yp2),(Xp2,Yp2),(Xp2,Yp2),(Xp2,Yp2))
#         solution=gridSearch(residual)
#         CTU_shifted[p,i]=cx[i]*solution/2.#.append(cx[i]*solution/2.)

# plt.plot(cx/cy,DCU_natural[0,:],'b')
# plt.plot(cx/cy,DCU_natural[1,:],'b--')
# plt.plot(cx/cy,CTU_natural[0,:],'r')
# plt.plot(cx/cy,CTU_natural[1,:],'r--')
# plt.grid()
# plt.show()

# plt.plot(cx/cy,DCU_shifted[0,:],'b')
# plt.plot(cx/cy,DCU_shifted[1,:],'b--')
# plt.plot(cx/cy,CTU_shifted[0,:],'r')
# plt.plot(cx/cy,CTU_shifted[1,:],'r--')
# plt.grid()
# plt.show()    

# export2DTeXFile('2dCFL_2ppcH_DCU.tex',cx/cy,np.array([DCU_natural,DCU_shifted]))
# export2DTeXFile('2dCFL_2ppcH_CTU.tex',cx/cy,np.array([CTU_natural,CTU_shifted]))

# print "=== Symmetric Vertical ==="
# Yp=np.array([-0.5,0.5])
# Xp=np.array([0.,0.])

# DCU_natural=[]
# CTU_natural=[]
# DCU_shifted=[]
# CTU_shifted=[]
# for i in range(samples):
#     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp))
#     solution=optimize.newton(residual,1.)
#     DCU_natural.append(cx[i]*solution/2.)
#     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp),(Xp,Yp))
#     solution=optimize.newton(residual,1.)
#     CTU_natural.append(cx[i]*solution/2.)

# Yp=np.array([-0.25,0.25])
# Xp=np.array([0.,0.])
# # Shift symmetrically
# for i in range(samples):
#     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp))
#     solution=optimize.newton(residual,1.)
#     DCU_shifted.append(cx[i]*solution/2.)
#     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp),(Xp,Yp))
#     solution=optimize.newton(residual,1.)
#     CTU_shifted.append(cx[i]*solution/2.)

# export2DTeXFile('2dCFL_2ppcV_DCU.tex',cx/cy,np.array([DCU_natural,DCU_shifted]))
# export2DTeXFile('2dCFL_2ppcV_CTU.tex',cx/cy,np.array([CTU_natural,CTU_shifted]))


############### 4PPC
# print "**************************************************************"
# print "******************  4PPC discretization **********************"
# print "**************************************************************"
# print "=== Symmetric ==="
# Xp=np.array([-0.5,0.5,0.5,-0.5])
# Yp=np.array([-0.5,-0.5,0.5,0.5])

# Xp2=np.array([-0.25,0.25,0.25,-0.25])
# Yp2=np.array([-0.25,-0.25,0.25,0.25])

# Xp3=np.array([-1.,1.,1.,-1.])
# Yp3=np.array([-1.,-1.,1.,1.])

# DCU_natural=np.zeros((1,samples))
# CTU_natural=np.zeros((1,samples))
# DCU_shifted=np.zeros((1,samples))
# CTU_shifted=np.zeros((1,samples))
# DCU_shifted2=np.zeros((1,samples))
# CTU_shifted2=np.zeros((1,samples))
# fig = plt.figure()
# ax1=plt.subplot2grid((1,2),(0,0))
# ax2=plt.subplot2grid((1,2),(0,1))
# col=['b','r','g','k','m','y','c']
# for p in range(1):
#     for i in range(samples):
#         residual=symbolResidual(p,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp))
#         solution=gridSearch(residual)
#         DCU_natural[p,i]=cx[i]*solution/2.#.append(cx[i]*solution/2.)
#         residual=symbolResidual(p,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp),(Xp,Yp))
#         solution=gridSearch(residual)
#         CTU_natural[p,i]=cx[i]*solution/2.#.append(cx[i]*solution/2.)

#         # Shift symmetrically
#         residual=symbolResidual(p,cx[i],cy,(Xp2,Yp2),(Xp2,Yp2),(Xp2,Yp2))
#         solution=gridSearch(residual)
#         DCU_shifted[p,i]=cx[i]*solution/2.#.append(cx[i]*solution/2.)
#         residual=symbolResidual(p,cx[i],cy,(Xp2,Yp2),(Xp2,Yp2),(Xp2,Yp2),(Xp2,Yp2))
#         solution=gridSearch(residual)
#         CTU_shifted[p,i]=cx[i]*solution/2.#.append(cx[i]*solution/2.)

#         # Shift symmetrically on nodes
#         residual=symbolResidual(p,cx[i],cy,(Xp3,Yp3),(Xp3,Yp3),(Xp3,Yp3))
#         solution=optimize.newton(residual,1.)
#         DCU_shifted2[p,i]=cx[i]*solution/2.#.append(cx[i]*solution/2.)
#         residual=symbolResidual(p,cx[i],cy,(Xp3,Yp3),(Xp3,Yp3),(Xp3,Yp3),(Xp3,Yp3))
#         solution=optimize.newton(residual,1.)
#         CTU_shifted2[p,i]=cx[i]*solution/2.#.append(cx[i]*solution/2.)
#     ax1.plot(cx/cy,DCU_natural[p,:],col[p],label='nat'+str(p))
#     ax1.plot(cx/cy,DCU_shifted[p,:],col[p],label='shift'+str(p))
#     ax1.plot(cx/cy,DCU_shifted2[p,:],col[p],label='shift2'+str(p))
#     ax2.plot(cx/cy,CTU_natural[p,:],col[p],label='nat'+str(p))
#     ax2.plot(cx/cy,CTU_shifted[p,:],col[p],linestyle='--',label='shift'+str(p))
#     ax2.plot(cx/cy,CTU_shifted2[p,:],col[p],linestyle='-.',label='shift2'+str(p))

# ax1.grid();ax2.grid()
# ax1.legend();ax2.legend()
# plt.show()

# # plt.plot(cx/cy,DCU_natural[0,:],'b')
# # plt.plot(cx/cy,DCU_natural[1,:],'b--')
# # plt.plot(cx/cy,CTU_natural[0,:],'r')
# # plt.plot(cx/cy,CTU_natural[1,:],'r--')
# # plt.grid()
# # plt.show()

# # plt.plot(cx/cy,DCU_shifted[0,:],'b')
# # plt.plot(cx/cy,DCU_shifted[1,:],'b--')
# # plt.plot(cx/cy,CTU_shifted[0,:],'r')
# # plt.plot(cx/cy,CTU_shifted[1,:],'r--')
# # plt.grid()
# # plt.show()    

# # DCU_natural=[]
# # CTU_natural=[]
# # DCU_shifted=[]
# # CTU_shifted=[]
# # for i in range(samples):
# #     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp))
# #     solution=gridSearch(residual)
# #     DCU_natural.append(cx[i]*solution/2.)
# #     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp),(Xp,Yp))
# #     solution=gridSearch(residual)
# #     CTU_natural.append(cx[i]*solution/2.)

# # Shift symmetrically
# # for i in range(samples):
# #     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp))
# #     solution=gridSearch(residual)
# #     DCU_shifted.append(cx[i]*solution/2.)
# #     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp),(Xp,Yp))
# #     solution=gridSearch(residual)
# #     CTU_shifted.append(cx[i]*solution/2.)

# export2DTeXFile('2dCFL_4ppcS_DCU.tex',cx/cy,np.array([DCU_natural,DCU_shifted,DCU_shifted2]))
# export2DTeXFile('2dCFL_4ppcS_CTU.tex',cx/cy,np.array([CTU_natural,CTU_shifted,CTU_shifted2]))

# pdb.set_trace()
# Xp=np.array([-0.5,0.5,0.5,-0.5])
# Yp=np.array([-0.5,-0.5,0.5,0.5])+0.25

# DCU_natural=[]
# CTU_natural=[]
# DCU_shifted=[]
# CTU_shifted=[]
# DCU_shifted2=[]
# CTU_shifted2=[]
# for i in range(samples):
#     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp))
#     solution=gridSearch(residual)
#     DCU_natural.append(cx[i]*solution/2.)
#     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp),(Xp,Yp))
#     solution=gridSearch(residual)
#     CTU_natural.append(cx[i]*solution/2.)

# Xp=np.array([-0.5,0.5,0.5,-0.5])
# Yp=np.array([-0.5,-0.5,0.5,0.5])-0.25
# # Shift symmetrically
# for i in range(samples):
#     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp))
#     solution=gridSearch(residual)
#     DCU_shifted.append(cx[i]*solution/2.)
#     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp),(Xp,Yp))
#     solution=gridSearch(residual)
#     CTU_shifted.append(cx[i]*solution/2.)

# export2DTeXFile('2dCFL_4ppcV_DCU.tex',cx/cy,np.array([DCU_natural,DCU_shifted]))
# export2DTeXFile('2dCFL_4ppcV_CTU.tex',cx/cy,np.array([CTU_natural,CTU_shifted]))


# Yp=np.array([-0.5,0.5,0.5,-0.5])
# Xp=np.array([-0.5,-0.5,0.5,0.5])+0.25

# DCU_natural=[]
# CTU_natural=[]
# DCU_shifted=[]
# CTU_shifted=[]
# for i in range(samples):
#     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp))
#     solution=gridSearch(residual)
#     DCU_natural.append(cx[i]*solution/2.)
#     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp),(Xp,Yp))
#     solution=gridSearch(residual)
#     CTU_natural.append(cx[i]*solution/2.)

# Yp=np.array([-0.5,0.5,0.5,-0.5])
# Xp=np.array([-0.5,-0.5,0.5,0.5])-0.25
# # Shift symmetrically
# for i in range(samples):
#     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp))
#     solution=gridSearch(residual)
#     DCU_shifted.append(cx[i]*solution/2.)
#     residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp),(Xp,Yp))
#     solution=gridSearch(residual)
#     CTU_shifted.append(cx[i]*solution/2.)



# export2DTeXFile('2dCFL_4ppcH_DCU.tex',cx/cy,np.array([DCU_natural,DCU_shifted]))
# export2DTeXFile('2dCFL_4ppcH_CTU.tex',cx/cy,np.array([CTU_natural,CTU_shifted]))

############### 9PPC
print "**************************************************************"
print "******************  9PPC discretization **********************"
print "**************************************************************"
print "=== Symmetric ==="
Xp=np.array([-2./3.,2./3.,2./3.,-2./3.,-2./3.,0.,2./3.,0.,0.])
Yp=np.array([-2./3.,-2./3.,2./3.,2./3.,0.,0.,0.,-2./3.,2./3.])

Xp2=np.array([-1./3.,1./3.,1./3.,-1./3.,-1./3.,0.,1./3.,0.,0.])
Yp2=np.array([-1./3.,-1./3.,1./3.,1./3.,0.,0.,0.,-1./3.,1./3.])


Xp3=np.array([-1.,1.,1.,-1.,-1.,0.,1.,0.,0.])
Yp3=np.array([-1.,-1.,1.,1.,0.,0.,0.,-1.,1.])

DCU_natural=[]
CTU_natural=[]
DCU_shifted=[]
CTU_shifted=[]
DCU_shifted2=[]
CTU_shifted2=[]

for i in range(samples):
    residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp))
    solution=gridSearch(residual)
    DCU_natural.append(cx[i]*solution/2.)
    residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp),(Xp,Yp))
    solution=gridSearch(residual)
    CTU_natural.append(cx[i]*solution/2.)
    
    # Shift symmetrically
    residual=symbolResidual(0,cx[i],cy,(Xp2,Yp2),(Xp2,Yp2),(Xp2,Yp2))
    solution=gridSearch(residual)
    DCU_shifted.append(cx[i]*solution/2.)
    residual=symbolResidual(0,cx[i],cy,(Xp2,Yp2),(Xp2,Yp2),(Xp2,Yp2),(Xp2,Yp2))
    solution=gridSearch(residual)
    CTU_shifted.append(cx[i]*solution/2.)

    # Shift symmetrically 2
    residual=symbolResidual(0,cx[i],cy,(Xp3,Yp3),(Xp3,Yp3),(Xp3,Yp3))
    solution=gridSearch(residual)
    DCU_shifted2.append(cx[i]*solution/2.)
    residual=symbolResidual(0,cx[i],cy,(Xp3,Yp3),(Xp3,Yp3),(Xp3,Yp3),(Xp3,Yp3))
    solution=gridSearch(residual)
    CTU_shifted2.append(cx[i]*solution/2.)

export2DTeXFile('2dCFL_9ppcS_DCU.tex',cx/cy,np.array([DCU_natural,DCU_shifted,DCU_shifted2]))
export2DTeXFile('2dCFL_9ppcS_CTU.tex',cx/cy,np.array([CTU_natural,CTU_shifted,CTU_shifted2]))

Xp=np.array([-2./3.,2./3.,2./3.,-2./3.,-2./3.,0.,2./3.,0.,0.])
Yp=np.array([-2./3.,-2./3.,2./3.,2./3.,0.,0.,0.,-2./3.,2./3.])+0.1


DCU_natural=[]
CTU_natural=[]
DCU_shifted=[]
CTU_shifted=[]
for i in range(samples):
    residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp))
    solution=gridSearch(residual)
    DCU_natural.append(cx[i]*solution/2.)
    residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp),(Xp,Yp))
    solution=gridSearch(residual)
    CTU_natural.append(cx[i]*solution/2.)

Xp=np.array([-2./3.,2./3.,2./3.,-2./3.,-2./3.,0.,2./3.,0.,0.])
Yp=np.array([-2./3.,-2./3.,2./3.,2./3.,0.,0.,0.,-2./3.,2./3.])-0.1
# Shift symmetrically
for i in range(samples):
    residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp))
    solution=gridSearch(residual)
    DCU_shifted.append(cx[i]*solution/2.)
    residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp),(Xp,Yp))
    solution=gridSearch(residual)
    CTU_shifted.append(cx[i]*solution/2.)

export2DTeXFile('2dCFL_9ppcV_DCU.tex',cx/cy,np.array([DCU_natural,DCU_shifted]))
export2DTeXFile('2dCFL_9ppcV_CTU.tex',cx/cy,np.array([CTU_natural,CTU_shifted]))


Xp=np.array([-2./3.,2./3.,2./3.,-2./3.,-2./3.,0.,2./3.,0.,0.])+0.1
Yp=np.array([-2./3.,-2./3.,2./3.,2./3.,0.,0.,0.,-2./3.,2./3.])

DCU_natural=[]
CTU_natural=[]
DCU_shifted=[]
CTU_shifted=[]
for i in range(samples):
    residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp))
    solution=gridSearch(residual)
    DCU_natural.append(cx[i]*solution/2.)
    residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp),(Xp,Yp))
    solution=gridSearch(residual)
    CTU_natural.append(cx[i]*solution/2.)

Xp=np.array([-2./3.,2./3.,2./3.,-2./3.,-2./3.,0.,2./3.,0.,0.])-0.1
Yp=np.array([-2./3.,-2./3.,2./3.,2./3.,0.,0.,0.,-2./3.,2./3.])
# Shift symmetrically
for i in range(samples):
    residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp))
    solution=gridSearch(residual)
    DCU_shifted.append(cx[i]*solution/2.)
    residual=symbolResidual(0,cx[i],cy,(Xp,Yp),(Xp,Yp),(Xp,Yp),(Xp,Yp))
    solution=gridSearch(residual)
    CTU_shifted.append(cx[i]*solution/2.)

export2DTeXFile('2dCFL_9ppcH_DCU.tex',cx/cy,np.array([DCU_natural,DCU_shifted]))
export2DTeXFile('2dCFL_9ppcH_CTU.tex',cx/cy,np.array([CTU_natural,CTU_shifted]))
