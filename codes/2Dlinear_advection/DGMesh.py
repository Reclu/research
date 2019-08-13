# !/usr/bin/python

import numpy as np
import pdb

""" 
The class Mesh builds a mesh object from a mesh built with gmsh. Objects instancied from this class contain:
- The table of coordinates of nodes: coor(x , y)
- The table of connectivity of the mesh: connec(node1,node2,node3,node4)
- The table of edges of the mesh: edges(dof 1, dof 2, dof 3, dof 4, cell L, cell R)
- The table of inverse connectivity: invConnec(Nb elements connected to node, connected element number )
- The table of edges connected to one element: cellToEdges(edge numbers)
- The table of nodes associated to each dofs: dofsToNodes(Dof number,node number,element,edge number)
- The table of corner nodes: cornerNodes(corner number,prev ghost dof,next ghost dof,prev edge,next edge)
- The table of dofs connected to one edge edgeToDofs(dof1, dof2, dof3, dof4, left cell, right cell) (bottom left, bottom right, top right, top left)
"""

class DGmesh:
    #Constructor
    def __init__(self,File,pointFile,cp,cs,rho):
        self.cp = cp ; self.cs=cs ; self.rho = rho
        self.transverse = False
        
        # Read mesh from gmsh file
        self.coor,self.connec=self.readMesh(File)
        self.connec,self.barycenter=self.checkConnecOrientation(self.coor,self.connec)
        self.nodeToCell = self.connecInverse(self.coor,self.connec)
        self.edges,self.edgeToDofs,self.edgeLength,self.dofToNodes,self.cellToEdges = self.tableOfEdges(self.coor,self.connec,self.nodeToCell)
        self.cornerNodes=self.createCornerNodes()

        self.ndofs = np.shape(self.cornerNodes)[0]+np.shape(self.dofToNodes)[0]
        
        # Export mesh for paraview
        self.exportCoor ,self.exportConnec = self.readOriginalMesh(pointFile)

    def findNeighborCells(self):
        dx = np.min(self.edgeLength[:,0])
        nCells=self.connec.shape[0]
        neighbors=-np.ones((nCells,3),dtype=int)
        for i in range(nCells):
            xG=self.barycenter[i,0];yG=self.barycenter[i,1];
            for k in range(nCells):
                xGk=self.barycenter[k,0];yGk=self.barycenter[k,1];
                if (abs(xG-(xGk+dx))<1.e-4) and (abs(yG-yGk)<1.e-4):
                    # left neighbor
                    neighbors[i,0]=k
                if (abs(xG-(xGk+dx))<1.e-4) and (abs(yG-(yGk+dx))<1.e-4):
                    # left-bottom neighbor
                    neighbors[i,2]=k
                if (abs(xG-xGk)<1.e-4) and (abs(yG-(yGk+dx))<1.e-4):
                    # bottom neighbor
                    neighbors[i,1]=k
        return neighbors
    def checkCellOrientation(self,coordinates,nodes):
        """ 
        GMSH does not use the same storage of nodes in the connectivity table. Changes are then required in the connectivities (the cross product of bottom boundary and left boundary direction vectors must be positive)
        The method also compute the barycenter of elements that is used to localize the material points within the mesh
        """
        X = np.zeros(4);Y = np.zeros(4);posNode=np.zeros((2,4))
        for idp,p in enumerate(nodes):
            X[idp] = coordinates[p,0];Y[idp] = coordinates[p,1]
            posNode[0,idp]=X[idp];posNode[1,idp]=Y[idp]
        #Test orientation: cross product
        N1N2 = posNode[:,1]-posNode[:,0]
        N1N4 = posNode[:,3]-posNode[:,0]
        buf = np.cross(N1N2,N1N4)
        if buf>0.0:
            pass
        elif buf<0.0:
            #Correction
            nodes[1],nodes[3]=nodes[3],nodes[1]
        xb= (coordinates[nodes[1],0]+coordinates[nodes[0],0])/2.
        yb= (coordinates[nodes[2],1]+coordinates[nodes[1],1])/2.
        barycenter=np.array([xb,yb])
        return nodes,barycenter

    def readOriginalMesh(self,File):
        """Read the material points set as a gmsh mesh and buid the tables coor and connec"""
        ## Skip preamble
        dataFile = open(File,"r")
        line=''
        while line!='$Nodes':
            line = dataFile.readline().strip('\n')
        Nnodes= int(dataFile.readline().strip('\n'))
        ## Initialization of coordinates table
        coor = np.zeros((Nnodes,2))
        for i in range(Nnodes):
            line=dataFile.readline().strip('\n')
            listLine = line.split(' ')
            coor[i,0]=listLine[1] ; coor[i,1]=listLine[2]
            if line=='$EndNodes':
                print "Caution: number of nodes in msh file not consistent with the number of coordinates given"
                break
        ## Read connectivity table and build the table of edges
        connec = np.zeros((1,4),dtype=int)
        i=0;numberOfEdges=0
        # Initialization of barycenter table
        while line!='$EndElements':
            line = dataFile.readline().strip('\n')
            listLine = line.split(' ')
            if len(listLine)<2:
                continue
            if (int(listLine[1]) == 3):
                connec  = np.resize(connec,(i+1,4))
                connec[i,0]=int(listLine[5])-1 ; connec[i,1]=int(listLine[6])-1
                connec[i,2]=int(listLine[7])-1 ; connec[i,3]=int(listLine[8])-1
                connec[i,:],bary = self.checkCellOrientation(coor,connec[i,:])
                i+=1
            else:
                continue
        dataFile.close()
        return coor,connec


    def computeVolumes(self):
        nCells = self.connec.shape[0]
        volume = np.zeros(nCells)
        for i in range(nCells):
            # norm of the cross product of edge 1 and 4 vectors
            e1 = self.cellToEdges[i,0] #edge 1 index
            n1e1 = self.edges[e1,0] # edge 1 first node index
            n2e1 = self.edges[e1,1] # edge 1 second node index
            v1 = np.array([self.coor[n1e1,0]-self.coor[n2e1,0],self.coor[n1e1,1]-self.coor[n2e1,1]])
            e4 = self.cellToEdges[i,3] #edge 4 index
            n1e4 = self.edges[e4,0] # edge 4 first node index
            n2e4 = self.edges[e4,1] # edge 4 second node index
            v4 = np.array([self.coor[n1e4,0]-self.coor[n2e4,0],self.coor[n1e4,1]-self.coor[n2e4,1]])
            # cross product
            volume[i]=np.abs(v1[0]*v4[1]-v1[1]*v4[0])
            
        return volume

    def readMesh(self,File):
        #Read a gmsh mesh and buid the tables coor, connec and boundary
        dataFile = open(File,"r")
        coor = np.zeros(2)
        connec = np.zeros((1,4),dtype=int)
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
                k=0;i=0
                while True:
                    listLine = line.split(' ')
                    if (np.int(listLine[1]) == 3):
                        connec  = np.resize(connec,(i+1,4))
                        for j in range(4):
                            connec[i,j]=np.int(listLine[5+j])-1   #tables start at 0 in python
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
        return coor,connec

    def checkConnecOrientation(self,coor,connec):
        barycenter=np.zeros((connec.shape[0],2))
        for i in range(connec.shape[0]):
            X = np.zeros(4);Y = np.zeros(4);posNode=np.zeros((2,4))
            for idp,p in enumerate(connec[i,:]):
                X[idp] = coor[p,0];Y[idp] = coor[p,1]
                posNode[0,idp]=X[idp];posNode[1,idp]=Y[idp]
            #Test orientation: cross product
            N1N2 = posNode[:,1]-posNode[:,0]
            N1N4 = posNode[:,3]-posNode[:,0]
            buf = np.cross(N1N2,N1N4)
            if buf>0.0:
                pass
            elif buf<0.0:
                #Correction
                connec[i,1],connec[i,3]=connec[i,3],connec[i,1]
            xb= (coor[connec[i,1],0]+coor[connec[i,0],0])/2.
            yb= (coor[connec[i,2],1]+coor[connec[i,1],1])/2.
            barycenter[i,:]=xb,yb
        return connec,barycenter
    
    def connecInverse(self,coor,connec):
        #Compute the table of inverse connectivity of the mesh
        Nnodes = coor.shape[0]
        NCells = connec.shape[0]
        Nmax = 5
        isInCell = np.zeros((NCells,Nnodes))
        invConnec = []
        dofsToNodes = []
        #Loop on Cells
        for i in range(NCells):
            isInCell[i,connec[i,:]]=True
        for i in range(Nnodes):
            location = np.where(isInCell[:,i]==True)[0]
            invConnec.append(location)
        return np.asarray(invConnec)

    def neighbourCell(self,cell,edge_index,N0,N1,connec,invConnec):
        Neighbour1 = invConnec[N0]
        Neighbour2 = invConnec[N1]
        same_cell=np.intersect1d(Neighbour1,Neighbour2)
        if len(same_cell)==1 :
            other_cell=-1
        else :
            other_cell=same_cell[np.where(same_cell!=cell)[0][0]]
        if (edge_index==0) or (edge_index==3) :
            left = other_cell
            right = cell
        else :
            left = cell
            right = other_cell
        return left,right

    def dofManager(self,e,dofToNodes,N0,N1,left,right,dof):
        dofs = np.array([[N0,left,e],[N0,right,e],[N1,right,e],[N1,left,e]])
        write = np.array([True,True,True,True])
        edge = np.zeros(4)
        # Check if dofs already exist
        for i in range(dofToNodes.shape[0]):
            for k in range(4):
                if np.array_equal(dofs[k,0:2],dofToNodes[i,0:2]) and dofs[k,1]!=-1:
                    write[k] = False
                    edge[k]=i
            if (write==False).all():
                break
        ToBeCreated = np.where(write==True)[0]
        for i in range(len(ToBeCreated)):
            dof+=1
            dofToNodes = np.resize(dofToNodes,(dof+1,3))
            dofToNodes[dof,:]=dofs[ToBeCreated[i],:]
            edge[ToBeCreated[i]]=dof
        return dofToNodes,dof,edge
    
    def tableOfEdges(self,coor,connec,invConnec):
        #Computes the table of edges of the mesh
        Nnodes = coor.shape[0]
        NCells = connec.shape[0]
        edgeToDofs = np.zeros((1,6),dtype=int)
        edgeLength = np.zeros((1,3))
        edge = np.zeros((1,4),dtype=int)
        dofToNodes=np.empty((1,3),dtype=int)
        cellToEdges = np.zeros((NCells,4),dtype=int)
        k=-1
        dof=-1
        #Loop on Cells
        for i in range(NCells):
            nbCellNodes = connec.shape[1]
            N0=connec[i,0];N1=connec[i,1];N2=connec[i,2];N3=connec[i,3]  
            edges = np.array([[N1,N0],[N1,N2],[N2,N3],[N0,N3]])
            for j in range(4):
                write = True
                BC = False
                N0 = edges[j,0]
                N1 = edges[j,1] 
                #Does this edge already exist in the table?        
                for p in range(edge.shape[0]):
                    if (((edge[p,0]==N0) and (edge[p,1]==N1))
                       or ((edge[p,0]==N1) and (edge[p,1]==N0))):
                        write = False
                        edge_num=p
                if (write):
                    k+=1
                    edgeToDofs = np.resize(edgeToDofs,(k+1,6))
                    edge = np.resize(edge,(k+1,4))
                    left,right=self.neighbourCell(i,j,N0,N1,connec,invConnec)
                    edge[k,0] = N0 ; edge[k,1] = N1
                    edge[k,2] = left ; edge[k,3] = right
                    # Do dofs already exist ?
                    dofToNodes,dof,edge_dofs = self.dofManager(k,dofToNodes,N0,N1,left,right,dof)
                    # Compute length and normal vector
                    dX = coor[N1,0]-coor[N0,0];dY = coor[N1,1]-coor[N0,1]
                    ds =  np.sqrt(dX**2+dY**2)
                    cT = np.round(dY/ds,decimals=9);sT=-np.round(dX/ds,decimals=9)
                    nx=np.abs(cT) ; ny=np.abs(sT)
                    edgeToDofs[k,:]= np.array([edge_dofs[0],edge_dofs[1],edge_dofs[2],edge_dofs[3],left,right])
                    #edgeLength.append(np.array([ds,nx,ny]))
                    edgeLength = np.resize(edgeLength,(k+1,3))
                    edgeLength[k,:]=np.array([ds,nx,ny])
                    # Increment cellToEdges table 
                    cellToEdges[i,j]=k
                else :
                    cellToEdges[i,j]=edge_num  
        return edge,edgeToDofs,edgeLength,dofToNodes,cellToEdges

    def buildBoundary(self,parent,matboundary):
        cellToEdges = self.cellToEdges
        edgeToDofs = self.edgeToDofs
        bound=[]
        edge=[]
        # bound = [cell , physical_line , edge 1, edge 2, ... ]
        # edge = [edge, cell , material point]
        cells=[]
        for i,Pt in enumerate(matboundary[:,0]) :
            cell = int(parent[Pt])
            cells.append(cell)
            line=matboundary[i,1]
            # find neighbour cells indices 
            cellL = edgeToDofs[cellToEdges[cell,3],4]
            cellR = edgeToDofs[cellToEdges[cell,1],5]
            cellT = edgeToDofs[cellToEdges[cell,2],5]
            cellB = edgeToDofs[cellToEdges[cell,0],4]
            # check if they're active or not
            Lactive = np.in1d(cellL,parent)
            Ractive = np.in1d(cellR,parent)
            Tactive = np.in1d(cellT,parent)
            Bactive = np.in1d(cellB,parent)
            #print Lactive,Ractive,Tactive,Bactive
            # position of boundary edge in the current cell
            loc_bc_edges = np.where(np.array([Bactive[0],Ractive[0],Tactive[0],Lactive[0]])==False)
            # indices of those edges in the global mesh
            glob_bc_edges = cellToEdges[cell,loc_bc_edges[0]]
            for k in range(len(glob_bc_edges)):
                write=True
                for p in range(len(edge)):
                    if (edge[p]==np.array([glob_bc_edges[k],cell,line])).all():
                        write = False
                        break
                if write :
                    edge.append(np.array([glob_bc_edges[k],cell,line]).astype(int))
        return bound,np.asarray(edge)


    
    def createCornerNodes(self):
        # If a node has only three connected dofs and at least one dof belonging to a -1 element, we have to create a corner node
        ## Find nodes corresponding to conditions mentioned above
        c=-1
        # corner_node(corner number,prev ghost node,next ghost node,prev edge, next edge)
        corner_node=np.empty((1,2),dtype=int)
        occur = np.zeros((np.shape(self.coor)[0],1),dtype=int)
        for i in range (np.shape(self.coor)[0]):
            dofs = np.where(self.dofToNodes[:,0]==i)[0]
            occur[i]=len(dofs) ; bound = len(np.where(self.dofToNodes[dofs,1]==-1)) 
            if occur[i]==3 and bound!=0:
                c_edge=np.where(self.edges[:,0:2]==i)[0]
                c_cell=np.intersect1d(np.where(self.cellToEdges==c_edge[0])[0],np.where(self.cellToEdges==c_edge[1])[0])
                index1 = np.where(self.cellToEdges[c_cell,:][0]==c_edge[0])[0][0]
                index2 = np.where(self.cellToEdges[c_cell,:][0]==c_edge[1])[0][0]
                if index1>index2 :
                    index=np.sort(np.array([index1,index2]))
                    c_edge[0],c_edge[1]=c_edge[1],c_edge[0]
                else :
                    index = np.array([index1,index2])
                # find which corner it is and build table
                if (index==[0,3]).all() : 
                    # bottom-left corner
                    c+=1
                    corner_node=np.resize(corner_node,(c+1,4))
                    corner_node[c,:]=self.edgeToDofs[c_edge[1],0],self.edgeToDofs[c_edge[0],3],c_edge[1],c_edge[0]
                elif (index==[0,1]).all() : 
                    # bottom-right corner
                    c+=1
                    corner_node=np.resize(corner_node,(c+1,4))
                    corner_node[c,:]=self.edgeToDofs[c_edge[0],0],self.edgeToDofs[c_edge[1],1],c_edge[0],c_edge[1]
                elif (index==[1,2]).all() : 
                    # top-right corner
                    c+=1
                    corner_node=np.resize(corner_node,(c+1,4))
                    corner_node[c,:]=self.edgeToDofs[c_edge[0],2],self.edgeToDofs[c_edge[1],1],c_edge[0],c_edge[1]
                else :
                    # top-left corner
                    c+=1
                    corner_node=np.resize(corner_node,(c+1,4))
                    corner_node[c,:]=self.edgeToDofs[c_edge[0],2],self.edgeToDofs[c_edge[1],3],c_edge[0],c_edge[1]
        return corner_node


    def detectCell(self,X,barycenter):
        # point to barycenter vector
        x_b = barycenter- X
        dist = np.sqrt(x_b[:,0]*x_b[:,0]+x_b[:,1]*x_b[:,1])
        parent=np.where(dist==np.min(dist))[0]
        return parent
        
    def buildApproximation(self,xp):
        # First idea : compute each cell barycenter coordinates. For a given material point, the parent cell is the one whose distance barycenter-point is minimum
        xn =self.coor
        connect=self.connec
        edges = self.edgeToDofs
        cell2edges=self.cellToEdges
        Nn=np.shape(self.dofToNodes)[0]
        Np=np.shape(xp)[0]
        Map=np.zeros((Nn,Np))
        Gradx=np.zeros((Nn,Np))
        Grady=np.zeros((Nn,Np))
        Dofs=np.zeros(Nn)
        Parent=np.zeros(Np)
        for Pt in range(Np):
            # detect parent element of current material point
            # A node on a cell boundary owns to the right one // upper one 
            parent=self.detectCell(xp[Pt,:],self.barycenter)
            nodes = connect[parent,:]
            dx=self.edgeLength[cell2edges[parent,0][0],0]
            dy=self.edgeLength[cell2edges[parent,1][0],0]
            # Dofs corresponding to parent element
            d=np.array([ edges[cell2edges[parent,0][0],2] , edges[cell2edges[parent,1][0],0] , edges[cell2edges[parent,2][0],0] , edges[cell2edges[parent,3][0],2] ])
            # Active dofs' indices storage
            Dofs[d]+=1
            # Active cell's indices storage
            Parent[Pt]=parent
            # Local coordinates of the point
            xi=2.*(xp[Pt,0]-xn[self.dofToNodes[d[0],0],0])/dx - 1.
            et=2.*(xp[Pt,1]-xn[self.dofToNodes[d[0],0],1])/dy - 1.
            # Shape functions and derivatives
            Shapes=np.array([(1-xi)*(1-et)/4.,(1+xi)*(1-et)/4.,(1+xi)*(1+et)/4.,(1-xi)*(1+et)/4.])
            dN=np.mat([[-(1-et)/4.,-(1-xi)/4.],\
                       [ (1-et)/4.,-(1+xi)/4.],\
                       [ (1+et)/4., (1+xi)/4.],\
                       [-(1+et)/4., (1-xi)/4.]])
            # Elementary gradient matrix
            Coord=np.zeros((2,4))
            Coord[0,:]= xn[nodes,0]
            Coord[1,:]= xn[nodes,1]
            #Coord=np.mat([ xn[nodes,0].T , xn[nodes,1].T])
            Jac=np.dot(dN.T,Coord.T)
            Jinv=np.linalg.inv(Jac)
            dN=np.dot(dN,Jinv)
            # Global mapping and gradient matrix
            dN1=np.asarray(dN[:,0])
            dN2=np.asarray(dN[:,1])
            Map[d,Pt]+=Shapes.T
            Gradx[d,Pt]+=dN1[:,0]
            Grady[d,Pt]+=dN2[:,0]
        d=[]
        for j in range(np.shape(Dofs)[0]):
            if Dofs[j]!=0:
                d.append(j)
        return Map,Gradx,Grady,d,Parent

    def computeFlux(self, U,dofs,boundary_interfaces,dt):
        """
        Evaluation of interface fluxes by solving Riemann problems
        Computation of volumic fluxes over each cell
        Assembling of interfaces and volumes contributions
        """
        fluxes = self.computeInterfaceFlux(U,boundary_interfaces,dt)
        fint = self.computeInternalForces(U,dofs)
        f = self.computeTotalForces(fint,fluxes,U)
        return f

    def findIndexes(self,index):
        indexprev=index[0][0]-1
        indexnext=index[0][0]-3
        return indexprev,indexnext


    def computeInterfaceFlux(self,U,boundary_interfaces,dt):
        cp = self.cp ; cs = self.cs ; rho = self.rho
        ninterfaces = np.shape(self.edges)[0]
        fluxes = np.zeros(np.shape(U)[0])
        
        for i in range(ninterfaces):
            LeftCell = self.edgeToDofs[i,4]
            RightCell = self.edgeToDofs[i,5]
            ## nodes connected through the current edge
            n1 = self.edgeToDofs[i,0] ; n2 = self.edgeToDofs[i,1]
            n3 = self.edgeToDofs[i,2] ; n4 = self.edgeToDofs[i,3]
            ## Length and outward normal vector to the left cell of the current edge
            ledge,nx,ny = self.edgeLength[i,:]
            
            UL = rho*(U[n1]+U[n4])/2.             
            UR = rho*(U[n2]+U[n3])/2.             
            # Positive flux on the edge
            apdq =  (cs*ny + cp*nx)*(UR - UL)
            flux =  (cs*ny + cp*nx)*UR - apdq
            
            # Interpolation from edge to the nodes (+ or - according to the outward normal vector to the cell)
            fluxes[n1]+=flux
            fluxes[n2]-=flux
            fluxes[n3]-=flux
            fluxes[n4]+=flux
            
            # TRANSVERSE RIEMANN PROBLEMS INSIDE DOMAIN
            if (self.transverse) :
                Modapdq = np.copy(apdq)
                # Split of the right-going fluctuations
                if RightCell!=-1 :
                    index= np.where(self.cellToEdges[RightCell,:]==i)
                    indexprev,indexnext = self.findIndexes(index)
                    
                    if index[0][0] == 0 :
                        Edge = self.cellToEdges[RightCell,indexnext]
                    else :
                        Edge = self.cellToEdges[RightCell,indexprev]
                    
                    ## Length and normal vector of the current edge
                    ledge,nx,ny = self.edgeLength[Edge,:]
                    
                    # Compute transverse contribution
                    bpapdq = 0.5*(cs*ny + cp*nx)*Modapdq*(dt/ledge)
                    
                    ## nodes of edge
                    n1 = self.edgeToDofs[Edge,0] ; n2 = self.edgeToDofs[Edge,1]
                    n3 = self.edgeToDofs[Edge,2] ; n4 = self.edgeToDofs[Edge,3]

                    # Assembling from edge to the nodes (+ or - according to the outward normal vector to the cell)
                    
                    fluxes[n1]-=bpapdq
                    fluxes[n2]+=bpapdq
                    fluxes[n3]+=bpapdq 
                    fluxes[n4]-=bpapdq 
            
        if (self.transverse):
            for q in range(boundary_interfaces.shape[0]):
                nL = int(boundary_interfaces[q,0]) 
                nR = int(boundary_interfaces[q,1])
                nx = boundary_interfaces[q,4]
                ny = boundary_interfaces[q,5]
                LeftEdge = int(boundary_interfaces[q,2])
                RightEdge = int(boundary_interfaces[q,3])

                
                if abs(nx)<1.e-12: # vertical ghost edge
                    n1=1. ; n2=0.
                else : # horizontal ghost edge
                    n1=0. ; n2=1.
                # RIEMANN PROBLEM
                UL = rho*U[nL] ; UR = rho*U[nR]
                apdq =  (cp*n1 + cs*n2)*(UR-UL)

                # Split of the right-going fluctuations and only keep left and bottom boundary 
                if RightEdge!=-1 :
                    if nR==self.edgeToDofs[RightEdge,0] :
                        # 4th boundary 
                        bapdq = 0.5*(dt/ledge)*cs*apdq
                        # Assembling the contribution on active nodes only
                        nodes = [self.edgeToDofs[RightEdge,1],self.edgeToDofs[RightEdge,2]]
                    elif nR==self.edgeToDofs[RightEdge,3]:
                        # 1st boundary
                        bapdq = 0.5*(dt/ledge)*cp*apdq
                        # Assembling the contribution on active nodes only
                        nodes = [self.edgeToDofs[RightEdge,1],self.edgeToDofs[RightEdge,2]]
                    fluxes[nodes]+=bapdq
                       
        # 0.5 comes from integration (dx = ledge*0.5*dxi)
        fluxes*=0.5*ledge
        
        return fluxes

    def computeFluxVectors(self,U,cp,cs):
        flux1 = cp*U
        flux2 = cs*U
        return flux1,flux2

    def computeInternalForces(self,U,dofs):
        cp = self.cp ; cs = self.cs ; rho=self.rho
        Sx = self.Kx ; Sy = self.Ky
        
        Nnodes = np.shape(U)[0]
        Nelem = (Nnodes-2)/2
        fint = np.zeros(Nnodes)
        flux1,flux2 = self.computeFluxVectors(U,cp,cs)        
        fint[dofs] = np.dot(Sx,flux1[dofs]) + np.dot(Sy,flux2[dofs])
        return fint

    def computeTotalForces(self,fint,fluxes,U):
        f=fint-fluxes
        return f

    def setMapping(self,Kx,Ky):
        self.Kx = Kx
        self.Ky = Ky
        
    def setTransverse(self,boolean):
        self.transverse = boolean
        
    def setParameters(self,cp,cs,rho):
        self.cp=cp
        self.cs=cs
        self.rho=rho
