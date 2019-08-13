#!/usr/bin/python

import numpy as np
import pdb


def exportParaview(nameFile,myMesh,**fields):
    dataFile = open(nameFile,"w")
    dataFile.write('<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n')
    dataFile.write('<UnstructuredGrid>\n')
    dataFile.write('<Piece NumberOfPoints=\"'+str(myMesh.exportCoor.shape[0])+'\" NumberOfCells=\"'+str(myMesh.exportConnec.shape[0])+'\">\n')
    dataFile.write('<Points>\n')
    dataFile.write('<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n')
    for i in range(myMesh.exportCoor.shape[0]):
        dataFile.write(str(myMesh.exportCoor[i,0])+' '+str(myMesh.exportCoor[i,1])+' 0 ')
    dataFile.write('\n')
    dataFile.write('</DataArray>\n')
    dataFile.write('</Points>\n')
    dataFile.write('<Cells>\n')
    dataFile.write('<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n')
    for i in range(myMesh.exportConnec.shape[0]):
        for j in range(myMesh.exportConnec.shape[1]):
            dataFile.write(str(myMesh.exportConnec[i,j])+' ')
    dataFile.write('\n')
    dataFile.write('</DataArray>\n')
    dataFile.write('<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n')
    offset = 0
    for i in range(myMesh.exportConnec.shape[0]):
        offset+=len(myMesh.exportConnec[i,:])
        dataFile.write(str(offset)+' ')
    dataFile.write('\n')
    dataFile.write('</DataArray>\n')
    dataFile.write('<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n')
    for i in range(myMesh.exportConnec.shape[0]):
        dataFile.write('9 ')
    dataFile.write('\n')
    dataFile.write('</DataArray>\n')
    dataFile.write('</Cells>\n')
    #Insert the writing of fields here
    dataFile.write('<CellData>\n')
    for name,vector in fields.items():
        writeCellField(name,vector,dataFile)
    #Close the file
    dataFile.write('</CellData>\n')
    dataFile.write('</Piece>\n')
    dataFile.write('</UnstructuredGrid>\n')
    dataFile.write('</VTKFile>\n')
    dataFile.close()


def writeCellField(name,scalarField,dataFile):
    dataFile.write('<DataArray type=\"Float32\" Name=\"'+str(name)+'\" format=\"ascii\">\n')
    for val in scalarField:
        dataFile.write(str(val)+' ')
    dataFile.write('\n')
    dataFile.write('</DataArray>\n') 
    
def exportPointData(nameFile,coor,**fields):
    dataFile = open(nameFile,"w")
    dataFile.write('<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n')
    dataFile.write('<UnstructuredGrid>\n')
    dataFile.write('<Piece NumberOfPoints=\"'+str(coor.shape[0])+'\" NumberOfCells=\"'+str(0)+'\">\n')
    dataFile.write('<Points>\n')
    dataFile.write('<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n')
    for i in range(coor.shape[0]):
        dataFile.write(str(coor[i,0])+' '+str(coor[i,1])+' 0 ')
    dataFile.write('\n')
    dataFile.write('</DataArray>\n')
    dataFile.write('</Points>\n')
    dataFile.write('<PointData  Scalar="scalar">\n')
    for name,vector in fields.items():
        writePointField(name,vector,dataFile)
    dataFile.write('</PointData>')

    dataFile.write('<Cells>\n')
    dataFile.write('<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n')
    dataFile.write('</DataArray>\n')
    dataFile.write('<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n')
    dataFile.write('</DataArray>\n')
    dataFile.write('<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n')
    dataFile.write('</DataArray>\n')
    dataFile.write('</Cells>\n')
    dataFile.write('</Piece>\n')
    dataFile.write('</UnstructuredGrid>\n')
    dataFile.write('</VTKFile>\n')
    dataFile.close()


def writePointField(name,scalarField,dataFile):
    dataFile.write('<DataArray type=\"Float32\" Name=\"'+str(name)+'\" format=\"ascii\">\n')
    for val in scalarField:
        dataFile.write(str(val)+' ')
    dataFile.write('\n')
    dataFile.write('</DataArray>\n') 

def exportDelaunay(nameFile,coor,connec,**fields):
    dataFile = open(nameFile,"w")
    dataFile.write('<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n')
    dataFile.write('<UnstructuredGrid>\n')
    dataFile.write('<Piece NumberOfPoints=\"'+str(coor.shape[0])+'\" NumberOfCells=\"'+str(connec.shape[0])+'\">\n')
    dataFile.write('<Points>\n')
    dataFile.write('<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n')
    for i in range(coor.shape[0]):
        dataFile.write(str(coor[i,0])+' '+str(coor[i,1])+' 0 ')
    dataFile.write('\n')
    dataFile.write('</DataArray>\n')
    dataFile.write('</Points>\n')
    dataFile.write('<Cells>\n')
    dataFile.write('<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n')
    for i in range(connec.shape[0]):
        for j in range(connec.shape[1]):
            dataFile.write(str(connec[i,j])+' ')
    dataFile.write('\n')
    dataFile.write('</DataArray>\n')
    dataFile.write('<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n')
    offset = 0
    for i in range(connec.shape[0]):
        offset+=len(connec[i,:])
        dataFile.write(str(offset)+' ')
    dataFile.write('\n')
    dataFile.write('</DataArray>\n')
    dataFile.write('<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n')
    for i in range(connec.shape[0]):
        dataFile.write('9 ')
    dataFile.write('\n')
    dataFile.write('</DataArray>\n')
    dataFile.write('</Cells>\n')
    # #Insert the writing of cell fields here
    # dataFile.write('<CellData>\n')
    # for name,vector in fields.items():
    #     writeCellField(name,vector,dataFile)
    # dataFile.write('</CellData>\n')
    dataFile.write('<PointData>\n')
    for name,vector in fields.items():
        if len(vector.shape)==1:
            writePointField(name,vector,dataFile)
        else:
            writeArrayField(name,vector,dataFile)
    dataFile.write('</PointData>\n')
    #Close the file
    dataFile.write('</Piece>\n')
    dataFile.write('</UnstructuredGrid>\n')
    dataFile.write('</VTKFile>\n')
    dataFile.close()

