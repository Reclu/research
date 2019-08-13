# !/usr/bin/python 

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from DGMesh import*
from exportParaview import*
import types
import pickle
import os
import pdb
from sympy import *
from scipy import optimize

###############################################################
def plot_domain(xp,xn):
    plt.plot(xp[:,0],xp[:,1],'ro',label='Material points',markersize=1.5)
    plt.plot(xn[:,0],xn[:,1],'b+',label='Mesh nodes')
    plt.grid()
    plt.axis('equal')
    plt.legend()
    plt.show()

def readMaterialPointSet(File):
    #Read a gmsh mesh and buid the tables coor, connec and boundary
    dataFile = open(File,"r")
    coor = np.zeros(2)
    matboundary = np.empty((1,2),dtype=int)
    Stop = False
    while True:
        line = dataFile.readline().strip('\n')
        if line == '$Nodes':
            NbNodes = np.int(dataFile.readline().strip('\n'))
            line = dataFile.readline().strip('\n')
            k=0
            while True:
                listLine = line.split(' ')
                coor = np.resize(coor,(k+1,2))
                coor[k,0] = listLine[1]
                coor[k,1] = listLine[2]
                k+=1
                line = dataFile.readline().strip('\n')
                if line == '$EndNodes':
                    break
        if line == '$Elements':
            NbElements =  np.int(dataFile.readline().strip('\n'))
            line = dataFile.readline().strip('\n')
            i=0
            while True:
                listLine = line.split(' ')
                if (np.int(listLine[1]) == 1):
                    for j in range(2):
                        write=True
                        line = np.array([np.int(listLine[5+j])-1,np.int(listLine[3])])
                        for p in range(matboundary.shape[0]):
                            if (line[0]==matboundary[p,0]) \
                                    and (line[1]==matboundary[p,1]):
                                write=False
                        if (write):
                            matboundary = np.resize(matboundary,(i+1,2))
                            matboundary[i,:]=line
                            i+=1
                    line = dataFile.readline().strip('\n')
                else :
                    line = dataFile.readline().strip('\n')
                if line == '$EndElements':
                    Stop = True
                    break
        if Stop:
            break
    dataFile.close()
    return coor,matboundary

def readBoundaryConditions(File):
    #Read a data file and associate a physical line to it's load parameters
    dataFile = open(File,"r")
    bcs = []
    Stop = False
    while True:
        line = dataFile.readline().strip('\n')
        k=0
        if line == '$BOUNDARY_CONDITIONS':
            line = dataFile.readline().strip('\n')
            while True:
                listLine = line.split(' ')
                if listLine[0][0]=="#":
                    line = dataFile.readline().strip('\n')
                    continue
                if listLine[1]=="'FREE'" :
                    bcs.append([np.int(listLine[0]),str(listLine[1]),0])
                else :
                    bcs.append([np.int(listLine[0]),str(listLine[1]),float(listLine[2])])
                line = dataFile.readline().strip('\n')
                if line == '$END':
                    Stop=True
                    break
        if Stop:
            break
    dataFile.close()
    return np.asarray(bcs)

def findNextAndPrev(mesh,cell,edge_num):
    # Adjacent cells to neighbour edges
    prev_LeftAdjacent = mesh.edgeToDofs[mesh.cellToEdges[cell,edge_num-1],4]
    prev_RightAdjacent = mesh.edgeToDofs[mesh.cellToEdges[cell,edge_num-1],5]
    next_LeftAdjacent = mesh.edgeToDofs[mesh.cellToEdges[cell,edge_num-3],4]
    next_RightAdjacent = mesh.edgeToDofs[mesh.cellToEdges[cell,edge_num-3],5]
    if prev_LeftAdjacent==cell:
        if prev_RightAdjacent==-1:
            prev_edge=-1
        else :
            prev_edge = mesh.cellToEdges[prev_RightAdjacent,edge_num]
    else :
        if prev_LeftAdjacent==-1:
            prev_edge=-1
        else :
            prev_edge = mesh.cellToEdges[prev_LeftAdjacent,edge_num]
    if next_LeftAdjacent==cell:
        if next_RightAdjacent==-1:
            next_edge=-1
        else :
            next_edge = mesh.cellToEdges[next_RightAdjacent,edge_num]
    else :
        if next_LeftAdjacent==-1:
            next_edge=-1
        else :
            next_edge = mesh.cellToEdges[next_LeftAdjacent,edge_num]
    return prev_edge,next_edge

def applyBoundaryCondition(u,mesh,boundary_conditions,parent,matboundary,bc_edges,rho,Map):
    bc_interface=[]
    inter = np.zeros((len(bc_edges)+4,6))
    count = 0
    done=np.empty(1)
    ed=0
    cellToEdges=mesh.cellToEdges
    edgeToDofs =mesh.edgeToDofs
    dofToNodes=mesh.dofToNodes
    corners=mesh.cornerNodes
    dof = np.shape(dofToNodes)[0]
    for q in range(len(bc_edges)):
        edge = bc_edges[q][0]
        cell = bc_edges[q][1]
        if np.in1d(edge,done):
            continue
        done=np.resize(done,(ed+1))
        done[ed]=edge
        ed+=1
        edge_num = np.where(cellToEdges[cell,:]==edge)[0][0]
        # Adjacent edges (direct sens in the cell)
        prev_edge,next_edge = findNextAndPrev(mesh,cell,edge_num)
        # Find the correct condition to apply
        number = len(np.where(bc_edges[:,0]==edge)[0])
        lines = bc_edges[np.where(bc_edges[:,0]==edge)[0],2]
        if number>1 : 
            # the edge belongs to a corner cell so it belongs to two physical lines
            if prev_edge==-1:
                # next_edge belongs to the same physical line 
                phys_line = bc_edges[np.where(bc_edges[:,0]==next_edge)[0][0],2]
                # preparing corner condition 
                phys_line_prev = lines[np.where(lines!=phys_line)[0][0]]
            else :
                # prev_edge belongs to the same physical line
                phys_line = bc_edges[np.where(bc_edges[:,0]==prev_edge)[0][0],2]
                # preparing corner condition 
                phys_line_next = lines[np.where(lines!=phys_line)[0][0]]
        else :
            phys_line = bc_edges[q][2]
        # If many material points linked to this edge, we have to build an interpolation
        # first : find the boundary points that belong to the same phys_line as the edge
        points = np.where(parent==cell)[0]
        same = matboundary[np.where(matboundary[:,1]==phys_line)[0],0]
        bc_points = np.intersect1d(points,same)
        
        # Remark : if the current edge belongs to a cell having 3 boundary edge ... It will be more difficult
        # Classical boundary condition on current edge
        # (1) identify ghost nodes
        if edgeToDofs[edge,4]==-1 : # left cell is a ghost cell
            g=[edgeToDofs[edge,0],edgeToDofs[edge,3]] # ghost nodes on this edge
            r=[edgeToDofs[edge,1],edgeToDofs[edge,2]] # real nodes on this edge
        else : # right cell is a ghost cell
            g=[edgeToDofs[edge,1],edgeToDofs[edge,2]]
            r=[edgeToDofs[edge,0],edgeToDofs[edge,3]]
        # (2) Apply boundary condition to ghost nodes according to phys_line
        condition=np.int(np.where(boundary_conditions[:,0].astype(int)==phys_line)[0][0])
        value = float(boundary_conditions[condition,2])/rho
        u[g[0]]=value
        u[g[1]]=value
        
        # (3) if corner cell, find corner node and apply condition on it
        # i = [nodeL,nodeR,edgeL,edgeR,nx_gedge,ny_gedge]
        # Rq : in global frame
        if prev_edge==-1 and next_edge!=-1:
            # "(2)" if corner is on the left of current edge when looked toward inside domain
            c = np.where(corners[:,3]==edge)[0][0] +dof#np.shape(mesh.dofToNodes)[0]
            #print edge,c
            if edge_num==0 :   # left-bottom corner (2)
                g = edgeToDofs[edge,3] 
                i1 = [c,g,-1,edge,0.,1.]  
                i2 = [edgeToDofs[edge,0],edgeToDofs[next_edge,3],edge,next_edge,0.,1.]
            elif edge_num==1 : # right-bottom corner (2)
                g = edgeToDofs[edge,3] 
                i1 = [c,g,-1,edge,-1.,0.]
                i2 = [edgeToDofs[edge,2],edgeToDofs[next_edge,1],edge,next_edge,-1.,0.]
            elif edge_num==2 : # right-top corner (2)
                g = edgeToDofs[edge,1] 
                i1 = [edgeToDofs[next_edge,1],edgeToDofs[edge,2],next_edge,edge,0.,-1.]
                i2 = [g,c,edge,-1,0.,-1.]
            elif edge_num==3 : # left-top corner (2)
                g = edgeToDofs[edge,3] 
                i1 = [edgeToDofs[next_edge,3],edgeToDofs[edge,0],next_edge,edge,1.,0.]
                i2 = [g,c,edge,-1,1.,0.]
            else :
                print "Warning : Unknown corner boundary conditions !!"
            condition=np.int(np.where(boundary_conditions[:,0].astype(int)==phys_line_prev)[0][0])
            value = float(boundary_conditions[condition,2])/rho
            u[c]= value
        elif next_edge==-1 and prev_edge!=-1:
            # "(1)" if corner is on the right of current edge when looked toward inside domain
            c = np.where(corners[:,2]==edge)[0][0]+np.shape(mesh.dofToNodes)[0]
            if edge_num==0 : # right-bottom corner (1)
                g = edgeToDofs[edge,0] 
                i1 = [edgeToDofs[prev_edge,0],edgeToDofs[edge,3],prev_edge,edge,0.,1.]
                i2 = [g,c,edge,-1,0.,1.]
            elif edge_num==1 : # right-top corner (1)
                g = edgeToDofs[edge,0] 
                i1 = [edgeToDofs[prev_edge,2],edgeToDofs[edge,1],prev_edge,edge,-1.,0.]
                i2 = [g,c,edge,-1,-1.,0.]
            elif edge_num==2 : # left-top corner (1)
                g = edgeToDofs[edge,2] 
                i1 = [c,g,-1,edge,0.,-1.]
                i2 = [edgeToDofs[edge,1],edgeToDofs[prev_edge,2],edge,prev_edge,0.,-1.]
            elif edge_num==3 : # left-bottom corner (1)
                g = edgeToDofs[edge,0] 
                i1 = [c,g,-1,edge,1.,0.]
                i2 = [edgeToDofs[edge,3],edgeToDofs[prev_edge,0],edge,prev_edge,1.,0.]
            else : 
                print "Warning : Unknown corner boundary conditions (2) !!"
            condition=np.int(np.where(boundary_conditions[:,0].astype(int)==phys_line_next)[0][0])
            value = float(boundary_conditions[condition,2])/rho
            u[c]= value
        
        else : # both!=-1 (or both==-1 if it was coded !)
            #classical boundary edge for transverse Riemann solver
            # i = [nodeL,nodeR,edgeL,edgeR,nx_realedge,ny_realedge]
            # Rq : in global frame
            if edge_num==0 :
                i1 = [edgeToDofs[prev_edge,0],edgeToDofs[edge,3],prev_edge,edge,0.,1.]
                i2 = [edgeToDofs[edge,0],edgeToDofs[next_edge,3],edge,next_edge,0.,1.]
            elif edge_num==1 :
                i1 = [edgeToDofs[prev_edge,2],edgeToDofs[edge,1],prev_edge,edge,-1.,0.]
                i2 = [edgeToDofs[edge,2],edgeToDofs[next_edge,1],edge,next_edge,-1.,0.]
            elif edge_num==2 :
                i1 = [edgeToDofs[next_edge,1],edgeToDofs[edge,2],next_edge,edge,0.,-1.]
                i2 = [edgeToDofs[edge,1],edgeToDofs[prev_edge,2],edge,prev_edge,0.,-1.]
            else :
                i1 = [edgeToDofs[next_edge,3],edgeToDofs[edge,0],next_edge,edge,1.,0.]
                i2 = [edgeToDofs[edge,3],edgeToDofs[prev_edge,0],edge,prev_edge,1.,0.]
        # Check if interface has already be written
        if not np.in1d(i1[0:2],inter[:,0:2]).all():
            inter[count,:]=i1
            count+=1
        if not np.in1d(i2[0:2],inter[:,0:2]).all():
            inter[count,:]=i2
            count+=1
        
    return u,inter      

def computeMaterialPointMass(rho,volume,parent):
    mass = np.zeros(parent.shape[0])
    for i in range(parent.shape[0]):
        cell = parent[i]
        number = len(np.where(parent==cell)[0])
        mass[i]=rho*volume[int(cell)]/number
    return mass

def UpdateState(dt,dofs,inter,Ml,u):
    f=mesh.computeFlux(u,dofs,inter,dt)
    u[dofs]+=dt*(f[dofs]/Ml[dofs])
    return u

def UpdateStateRK2(dt,dofs,inter,Ml,u):
    k1=np.zeros(np.shape(u)[0])
    k2=np.zeros(np.shape(u)[0])
    # first step : compute flux and intermediate state w
    f=mesh.computeFlux(u,dofs,inter,dt)
    k1[dofs]=dt*(f[dofs]/Ml[dofs])
    w = u+k1*0.5
    w,inter = applyBoundaryCondition(w,mesh,boundary_conditions,cells,MatBoundary,bc_edges,rho,Map)
    # second step : compute flux and update U
    f=mesh.computeFlux(w,dofs,inter,dt)
    k2[dofs]=dt*(f[dofs]/Ml[dofs])
    u+=k2
    return u


def UpdateStateRK2TVD(dt,dofs,inter,Ml,u):
    k1=np.copy(u)
    # first step : compute flux and intermediate state w
    f=mesh.computeFlux(u,dofs,inter,dt)
    k1[dofs]+=dt*(f[dofs]/Ml[dofs])
    # second step : compute flux and update U
    f=mesh.computeFlux(k1,dofs,inter,dt)
    k2=np.copy(k1)
    k2[dofs]+=dt*(f[dofs]/Ml[dofs])
    u=0.5*(u+k2)
    return u

def buildDiscreteOperator(mesh,Map,gradx,grady,parent,invParent,neighbors,transverse,cx,cy,dx):
    dt = symbols('dt')
    H_matrix=zeros(Map.shape[1])
    shape_functions=lambda x,y: np.array([(1.-x)*(1.-y)/4.,(1.+x)*(1.-y)/4.,(1.+x)*(1.+y)/4.,(1.-x)*(1.+y)/4.])
    
    ## shape functions evaluated at edges centers to weight fluxes contributions
    ## o -- 3 -- o
    ## |         |
    ## 4         2
    ## |         |
    ## o -- 1 -- o
    shapeOnEdge=shape_functions(np.array([0.,1.,0.,-1.]),np.array([-1.,0.,1.,0.]))
    for p in range(np.shape(Map)[1]):
        elem = int(parent[p])
        sharing = invParent[elem]
        left = invParent[neighbors[elem,0]]
        bottom = invParent[neighbors[elem,1]]
        bottomLeft = invParent[neighbors[elem,2]]
        ## List of dofs connected through the current cell
        
        dofs=np.array([ mesh.edgeToDofs[mesh.cellToEdges[elem,0],2] , mesh.edgeToDofs[mesh.cellToEdges[elem,1],0] , mesh.edgeToDofs[mesh.cellToEdges[elem,2],0] , mesh.edgeToDofs[mesh.cellToEdges[elem,3],2] ])
        ## Contributions of material points sharing the same cell
        for k in sharing:
            for i,node in enumerate(dofs):
                weightings=Map[node,p]/np.sum(Map[node,:])
                H_matrix[p,k]+=weightings*Map[node,k]
                for j in dofs:
                    H_matrix[p,k]+=dt*weightings*(Map[j,k]/np.sum(Map[j,:]))*(cx*np.dot(gradx[node,:],Map[j,:]) + cy*np.dot(grady[node,:],Map[j,:]))

                # Contributions of edges 2 and 3
                H_matrix[p,k]-=0.5*(dt/dx)*weightings*shapeOnEdge[i,1]*len(sharing)*cx*(Map[dofs[1],k]/np.sum(Map[dofs[1],:])+Map[dofs[2],k]/np.sum(Map[dofs[2],:]))
                H_matrix[p,k]-=0.5*(dt/dx)*weightings*shapeOnEdge[i,2]*len(sharing)*cy*(Map[dofs[2],k]/np.sum(Map[dofs[2],:])+Map[dofs[3],k]/np.sum(Map[dofs[3],:]))
                # Transverse contributions
                if transverse:
                    # across edge 2
                    H_matrix[p,k]+= 0.25*(dt/dx)**2*weightings*shapeOnEdge[i,1]*len(sharing)*cx*cy*(Map[dofs[0],k]/np.sum(Map[dofs[0],:])+Map[dofs[1],k]/np.sum(Map[dofs[1],:]))
                    # across edge 3
                    H_matrix[p,k]+= 0.25*(dt/dx)**2*weightings*shapeOnEdge[i,2]*len(sharing)*cx*cy*(Map[dofs[0],k]/np.sum(Map[dofs[0],:])+Map[dofs[3],k]/np.sum(Map[dofs[3],:]))
        ## Contributions of material points of left cell 
        leftElem=neighbors[elem,0]
        for k in left:
            if leftElem==-1:
                break
            ## List of dofs connected through the left cell
            leftDofs=np.array([ mesh.edgeToDofs[mesh.cellToEdges[leftElem,0],2] , mesh.edgeToDofs[mesh.cellToEdges[leftElem,1],0] , mesh.edgeToDofs[mesh.cellToEdges[leftElem,2],0] , mesh.edgeToDofs[mesh.cellToEdges[leftElem,3],2] ])
        
            for i,node in enumerate(dofs):
                weightings=Map[node,p]/np.sum(Map[node,:])
                ## edge 4 contribution
                H_matrix[p,k]+= 0.5*(dt/dx)*weightings*shapeOnEdge[i,3]*len(sharing)*cx*(Map[leftDofs[1],k]/np.sum(Map[leftDofs[1],:])+Map[leftDofs[2],k]/np.sum(Map[leftDofs[2],:]))
                if transverse:
                    ## edge 4 contribution
                    H_matrix[p,k]-=0.25*(dt/dx)**2*weightings*shapeOnEdge[i,3]*len(sharing)*cx*cy*(Map[leftDofs[0],k]/np.sum(Map[leftDofs[0],:])+Map[leftDofs[1],k]/np.sum(Map[leftDofs[1],:]))
                    ## edge 3 contribution
                    H_matrix[p,k]-=0.25*(dt/dx)**2*weightings*shapeOnEdge[i,2]*len(sharing)*cy*cx*(Map[leftDofs[1],k]/np.sum(Map[leftDofs[1],:])+Map[leftDofs[2],k]/np.sum(Map[leftDofs[2],:]))
        ## Contributions of material points of bottom cell
        bottomElem=neighbors[elem,1]
        for k in bottom:
            if bottomElem==-1:
                break
            ## List of dofs connected through the bottom cell
            bottomDofs=np.array([ mesh.edgeToDofs[mesh.cellToEdges[bottomElem,0],2] , mesh.edgeToDofs[mesh.cellToEdges[bottomElem,1],0] , mesh.edgeToDofs[mesh.cellToEdges[bottomElem,2],0] , mesh.edgeToDofs[mesh.cellToEdges[bottomElem,3],2] ])
            for i,node in enumerate(dofs):
                weightings=Map[node,p]/np.sum(Map[node,:])
                ## edge 1 contribution
                H_matrix[p,k]+= 0.5*(dt/dx)*weightings*shapeOnEdge[i,0]*len(sharing)*cy*(Map[bottomDofs[2],k]/np.sum(Map[bottomDofs[2],:])+Map[bottomDofs[3],k]/np.sum(Map[bottomDofs[3],:]))
                if transverse:
                    ## edge 1 contribution
                    H_matrix[p,k]-=0.25*(dt/dx)**2*weightings*shapeOnEdge[i,0]*len(sharing)*cy*cx*(Map[bottomDofs[0],k]/np.sum(Map[bottomDofs[0],:])+Map[bottomDofs[3],k]/np.sum(Map[bottomDofs[3],:]))
                    ## edge 2 contribution
                    H_matrix[p,k]-=0.25*(dt/dx)**2*weightings*shapeOnEdge[i,1]*len(sharing)*cx*cy*(Map[bottomDofs[2],k]/np.sum(Map[bottomDofs[2],:])+Map[bottomDofs[3],k]/np.sum(Map[bottomDofs[3],:]))
        if transverse:
            ## Contributions of material points of bottom-left cell
            bottomLeftElem=neighbors[elem,2]
            for k in bottomLeft:
                if bottomLeftElem==-1: break
                ## List of dofs connected through the bottom-left cell
                bottomLeftDofs=np.array([ mesh.edgeToDofs[mesh.cellToEdges[bottomLeftElem,0],2] , mesh.edgeToDofs[mesh.cellToEdges[bottomLeftElem,1],0] , mesh.edgeToDofs[mesh.cellToEdges[bottomLeftElem,2],0] , mesh.edgeToDofs[mesh.cellToEdges[bottomLeftElem,3],2] ])
                for i,node in enumerate(dofs):
                    weightings=Map[node,p]/np.sum(Map[node,:])
                    ## edge 1 contribution
                    H_matrix[p,k]+=0.25*(dt/dx)**2*weightings*shapeOnEdge[i,0]*len(sharing)*cy*cx*(Map[bottomLeftDofs[1],k]/np.sum(Map[bottomLeftDofs[1],:])+Map[bottomLeftDofs[2],k]/np.sum(Map[bottomLeftDofs[2],:]))
                    ## edge 4 contribution
                    H_matrix[p,k]+=0.25*(dt/dx)**2*weightings*shapeOnEdge[i,3]*len(sharing)*cx*cy*(Map[bottomLeftDofs[2],k]/np.sum(Map[bottomLeftDofs[2],:])+Map[bottomLeftDofs[3],k]/np.sum(Map[bottomLeftDofs[3],:]))
    

    Hoperator = lambdify((dt),H_matrix)
    return H_matrix,Hoperator

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

def computeCriticalCFL(Map,H,invParent,parent,neighbor,cx,cy,dx):
    dt=symbols('dt')
    sol=[]
    for cell in range(np.shape(neighbor)[0]):
        if (neighbor[cell,:]==-1).any(): continue
        for p in invParent[cell]:
            res = 0
            for k in range(np.shape(Map)[1]):
                res+=np.abs(H[p,k])
            # solve the residual
            residual=lambdify((dt),res-1.)
            solution=gridSearch(residual)
            #solution=optimize.root(residual,1.,method='hybr',options={'xtol':1.e-4}).x[0]
            if abs(residual(solution))>1.e-3: print "CAUTION: residual norm after solution is", abs(residual(solution))
            sol.append(solution*max(cx,cy)/dx)
        break
    Courant=min(sol)
    print "Critical Courant number set to ",Courant
    return Courant


def UpdateStateDiscreteOperator(U,H,BC,dt,dx,cx,cy,parent,invParent,transverse,neighbor,mesh,Map,gradx,grady):
    
    U_updated=np.zeros(np.shape(U))
    for p in range(np.shape(U)[0]):
        for k in range(np.shape(U)[0]):
            U_updated[p]+=H[p,k]*U[k]
    ## next, enforce the BC at the left and bottom boundaries
    
    Left=np.where(neighbor[:,0]==-1)[0]
    Bottom=np.where(neighbor[:,1]==-1)[0]
    BottomLeft=np.where(neighbor[:,2]==-1)[0]
    shape_functions=lambda x,y: np.array([(1.-x)*(1.-y)/4.,(1.+x)*(1.-y)/4.,(1.+x)*(1.+y)/4.,(1.-x)*(1.+y)/4.])
    shapeOnEdge=shape_functions(np.array([0.,1.,0.,-1.]),np.array([-1.,0.,1.,0.]))
    
    for p in range(np.shape(Map)[1]):
        elem = int(parent[p])
        sharing = invParent[elem]
        
        left = neighbor[elem,0]
        bottom = neighbor[elem,1]
        bottomLeft = neighbor[elem,2]
        if (left!=-1 and bottom!=-1 and bottomLeft!=-1):continue
        ## List of dofs connected through the current cell
        dofs=np.array([ mesh.edgeToDofs[mesh.cellToEdges[elem,0],2] , mesh.edgeToDofs[mesh.cellToEdges[elem,1],0] , mesh.edgeToDofs[mesh.cellToEdges[elem,2],0] , mesh.edgeToDofs[mesh.cellToEdges[elem,3],2] ])
        
        for k in sharing:
            Hpk=0.
            if left==-1:
                for i,node in enumerate(dofs):
                    weightings=Map[node,p]/np.sum(Map[node,:])
                    ## edge 4 contribution
                    Hpk+= 0.5*(dt/dx)*weightings*shapeOnEdge[i,3]*len(sharing)*cx*(Map[dofs[1],k]/np.sum(Map[dofs[1],:])+Map[dofs[2],k]/np.sum(Map[dofs[2],:]))
                    if transverse:
                        Hpk-=0.25*(dt/dx)**2*weightings*shapeOnEdge[i,3]*len(sharing)*cx*cy*(Map[dofs[0],k]/np.sum(Map[dofs[0],:])+Map[dofs[1],k]/np.sum(Map[dofs[1],:]))
                        ## edge 3 contribution
                        Hpk-=0.25*(dt/dx)**2*weightings*shapeOnEdge[i,2]*len(sharing)*cy*cx*(Map[dofs[1],k]/np.sum(Map[dofs[1],:])+Map[dofs[2],k]/np.sum(Map[dofs[2],:]))
                U_updated[p]+=Hpk*BC
            Hpk=0.
            if bottom==-1:
                for i,node in enumerate(dofs):
                    weightings=Map[node,p]/np.sum(Map[node,:])
                    ## edge 1 contribution
                    Hpk+= 0.5*(dt/dx)*weightings*shapeOnEdge[i,0]*len(sharing)*cy*(Map[dofs[2],k]/np.sum(Map[dofs[2],:])+Map[dofs[3],k]/np.sum(Map[dofs[3],:]))
                    if transverse:
                        Hpk-=0.25*(dt/dx)**2*weightings*shapeOnEdge[i,0]*len(sharing)*cy*cx*(Map[dofs[0],k]/np.sum(Map[dofs[0],:])+Map[dofs[3],k]/np.sum(Map[dofs[3],:]))
                        ## edge 2 contribution
                        Hpk-=0.25*(dt/dx)**2*weightings*shapeOnEdge[i,1]*len(sharing)*cx*cy*(Map[dofs[2],k]/np.sum(Map[dofs[2],:])+Map[dofs[3],k]/np.sum(Map[dofs[3],:]))
                U_updated[p]+=Hpk*0.#U[k]
            Hpk=0.
            if bottomLeft==-1 and transverse :
                for i,node in enumerate(dofs):
                    weightings=Map[node,p]/np.sum(Map[node,:])
                    ## edge 1 contribution
                    Hpk+=0.25*(dt/dx)**2*weightings*shapeOnEdge[i,0]*len(sharing)*cy*cx*(Map[dofs[1],k]/np.sum(Map[dofs[1],:])+Map[dofs[2],k]/np.sum(Map[dofs[2],:]))
                    ## edge 4 contribution
                    Hpk+=0.25*(dt/dx)**2*weightings*shapeOnEdge[i,3]*len(sharing)*cx*cy*(Map[dofs[2],k]/np.sum(Map[dofs[2],:])+Map[dofs[3],k]/np.sum(Map[dofs[3],:]))
                if (neighbor[elem,0]==-1 and neighbor[elem,1]!=-1):
                    U_updated[p]+=Hpk*BC
    
    return U_updated

def computeLpNorm(Unum,Uexact,dx,p):
    return (((dx**2*((np.abs(Unum-Uexact))**p)).sum())**(1.0/p))

def computeRelativeError(Unum,Uexact,dx,p):
    return (computeLpNorm(Unum,Uexact,dx,p)/computeLpNorm(np.zeros(len(Unum)),Uexact,dx,p))

###############################################################

# Mesh file
#############################################
MeshFile = 'geometry/mesh_4ppc.msh'
MeshFile = 'geometry/mesh_9ppc.msh'
#############################################

# Material points file
#############################################
MPoints_set = 'geometry/MPplate.msh'
#############################################

# Load file
#############################################
LoadFile='domain.dat'
#############################################

# MATERIAL PROPERTIES
rho = 1.

cp = 10. # Horizontal wave speed
cs = 5. # Vertical wave speed


mesh = DGmesh(MeshFile,MPoints_set,cp,cs,rho)
#print mesh.dofToNodes

# MATERIAL POINTS BUIDLING
print "Building material domain ...", 
xp,MatBoundary = readMaterialPointSet(MPoints_set)
print "Material domain built !"

#plot_domain(xp,mesh.coor)


# LOADING PARAMETERS
boundary_conditions=readBoundaryConditions(LoadFile)

# Minimum cell length to compute CFL number
dx = np.min(mesh.edgeLength[:,0])
L = np.max(xp[:,0])-np.min(xp[:,0])



## Counters
T=0.
n=0

## Options
update_position = False
transverse = True
mesh.setTransverse(transverse)

# APPROXIMATION BUILDING
print "Building approximation ...", 
Map,Gradx,Grady,dofs,cells = mesh.buildApproximation(xp)
print " Approximation built ! "
invParent=[]
for i in range(np.shape(mesh.connec)[0]):
    invParent.append(np.where(cells==i)[0])

# MATERIAL POINT'S MASS 
volume = mesh.computeVolumes()
mass = computeMaterialPointMass(rho,volume,cells)

neighbors=mesh.findNeighborCells()


if not os.path.exists("Hsyms_"+str(MeshFile[14:-4])+".npy"):
    print "Building the discrete operators..."
    Hsym,HOperator= buildDiscreteOperator(mesh,Map,Gradx,Grady,cells,invParent,neighbors,transverse,cp,cs,dx)
    File=open("Hsyms_"+str(MeshFile[14:-4])+".npy","w")
    pickle.dump(Hsym,File)
    File.close()
    print "Discrete operators built!"
else:
    print "Reading the discrete operators..."
    File = open("Hsyms_"+str(MeshFile[14:-4])+".npy","r")
    Hsym = pickle.load(File)
    File.close()
    dt=symbols('dt')
    HOperator=lambdify((dt),Hsym)
    print "Discrete operators read!"

CFL = computeCriticalCFL(Map,Hsym,invParent,cells,neighbors,cp,cs,dx)
# ALGORITHMIC PARAMETERS
#CFL = 0.7#/(1.+cs/cp)
dt=CFL*dx/max(cp,cs) 
tfinal=1.*L/max(cp,cs)
inc=round(tfinal/dt)
t_order= 1


bc_cells,bc_edges = mesh.buildBoundary(cells,MatBoundary)

# FIELDS CREATION
## Material points
Mp=xp.shape[0]
U = np.zeros(Mp)
Uh = np.zeros(Mp)

Md=mass*np.eye(Mp)

## Nodes
u = np.zeros(mesh.ndofs)

# FIELDS STORAGE ALLOCATION
Quantity=np.zeros((Mp,int(inc)+2))
Quantityh=np.zeros((Mp,int(inc)+2))
iteration=np.zeros(int(inc)+2)

Quantity[:,0]=rho*U
Quantityh[:,0]=rho*U
iteration[0]=T

mg=np.dot(np.dot(Map,Md),Map.T)
md=np.sum(mg,axis=1)
Kx=np.dot(np.dot(Gradx[dofs,:],Md),Map[dofs,:].T)
Ky=np.dot(np.dot(Grady[dofs,:],Md),Map[dofs,:].T)

# FLUX EVALUATOR CREATION
mesh.setMapping(Kx,Ky)

left=np.where(neighbors[:,0]==-1)[0]
for cell in left:
    matPoints=invParent[cell]
    U[matPoints]=10.
    Uh[matPoints]=10.
nameFile = 'results/dgmpm'+str(n)+'.vtu'
#exportParaview(nameFile,mesh,stress_11=S11[:,n],stress_22=S22[:,n],stress_12=S12[:,n],velocity_x=V1[:,n],velocity_y=V2[:,n])
exportDelaunay(nameFile,mesh.exportCoor ,mesh.exportConnec,q=U*rho,error=abs(U-Uh)*rho)
nameFile = 'results/dgmpmD'+str(n)+'.vtu'
exportDelaunay(nameFile,mesh.exportCoor ,mesh.exportConnec,q=Uh*rho,error=abs(U-Uh)*rho)

# COMPUTATION BEGINING
print '... computing ...'
while T<tfinal:
    
    # Mapping from material points to nodes
    u_glob = np.dot(Map,np.dot(Md,U))
    u[dofs]=u_glob[dofs]/md[dofs]
    
    
    # Boundary conditions
    u,inter = applyBoundaryCondition(u,mesh,boundary_conditions,cells,MatBoundary,bc_edges,rho,Map)
    
    # Conservative update
    if t_order==1 :
        u=UpdateState(dt,dofs,inter,md,u)
    if t_order==2 :
        u=UpdateStateRK2(dt,dofs,inter,md,u)
    if t_order==3 :
        u=UpdateStateRK2TVD(dt,dofs,inter,md,u)
    
    Uh=UpdateStateDiscreteOperator(Uh,HOperator(dt),0.,dt,dx,cp,cs,cells,invParent,transverse,neighbors,mesh,Map,Gradx,Grady)
    
    # Mapping back to the material points
    U=np.dot(Map[dofs,:].T,u[dofs])


    print 'Increment =', n, 't = ', T,' s.'
    n+=1
    T+=dt
    Quantity[:,n] = rho*U
    Quantityh[:,n] = rho*Uh
    
    
    u = np.zeros(mesh.ndofs)
    Error = abs(Quantity[:,n]-Quantityh[:,n])
    nameFile = 'results/dgmpm'+str(n)+'.vtu'
    #exportPointData(nameFile,xp,q=Quantity[:,n])
    exportDelaunay(nameFile,mesh.exportCoor ,mesh.exportConnec,q=Quantity[:,n],error=Error)
    nameFile = 'results/dgmpmD'+str(n)+'.vtu'
    exportDelaunay(nameFile,mesh.exportCoor ,mesh.exportConnec,q=Quantityh[:,n],error=Error)

## Compute the error between the two numerical solutions
error = computeRelativeError(Quantity[:,n],Quantityh[:,n],dx,2)
print "Error between the two numerical procedures: ",error

    


