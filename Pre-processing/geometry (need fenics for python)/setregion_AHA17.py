import sys
sys.path.append("/mnt/home/caichen3/lab/")
import os
import vtk
import vtk_py as vtk_py
from dolfin import *
import numpy as np
import math as math

def partitionmeshforEP(ugrid, nsector=6, nz=2, meas=None, center=None, xaxis=[1,0]):

    if(center == None):
    	midx = 0.5*(ugrid.GetBounds()[0] + ugrid.GetBounds()[1])
    	midy = 0.5*(ugrid.GetBounds()[2] + ugrid.GetBounds()[3])
    	center = [midx, midy]
    
    cellcenter = vtk.vtkCellCenters()
    cellcenter.SetInputData(ugrid)
    cellcenter.Update()
    
    zpartition = np.linspace(ugrid.GetBounds()[4], ugrid.GetBounds()[5], nz+1)
    print "zpartition =", zpartition
    apartition = np.linspace(-math.pi, math.pi, nsector+1)
    regid = vtk.vtkIntArray()
    data = vtk.vtkFloatArray()

    for cellid in range(0, ugrid.GetNumberOfCells()):
     	x = cellcenter.GetOutput().GetPoints().GetPoint(cellid)[0]
     	y = cellcenter.GetOutput().GetPoints().GetPoint(cellid)[1]
     	z = cellcenter.GetOutput().GetPoints().GetPoint(cellid)[2]

        if -1.95416705 < z < 0:
                apartition = np.linspace(-math.pi, math.pi, nsector+1)
    	        zloc = np.argmax(zpartition>z)
    	# Determine position in theta direction
                #print "zloc = ", zloc
    	        norm = np.linalg.norm([(x - midx), (y - midy)])
    	        angle = np.arctan2((y - midy)/norm, (x - midx)/norm)
    	        sloc = np.argmax(apartition>angle)
                if sloc == 1:
                        regloc = (zloc-3)*nsector + 3
                elif sloc == 2:
                        regloc = (zloc-3)*nsector + 4
                elif sloc == 3:
                        regloc = (zloc-3)*nsector + 5
                elif sloc == 4:
                        regloc = (zloc-3)*nsector + 6
                elif sloc == 5:
                        regloc = (zloc-3)*nsector + 1
                else:
                        regloc = (zloc-3)*nsector + 2                

    	        #regloc = (zloc-3)*nsector + sloc 
    	        regid.InsertNextValue(regloc)
                        #print coordinate with regiona id
                        #print x, y, z, regloc  
                       
        elif -3.9083341 < z < -1.95416705:
                apartition = np.linspace(-math.pi, math.pi, nsector+1)
                zloc = np.argmax(zpartition>z)
        # Determine position in theta direction
                #print "zloc = ", zloc
                norm = np.linalg.norm([(x - midx), (y - midy)])
                angle = np.arctan2((y - midy)/norm, (x - midx)/norm)
                sloc = np.argmax(apartition>angle)
                
                if sloc == 1:
                        regloc = (zloc-1)*nsector + 3
                elif sloc == 2:
                        regloc = (zloc-1)*nsector + 4
                elif sloc == 3:
                        regloc = (zloc-1)*nsector + 5
                elif sloc == 4:
                        regloc = (zloc-1)*nsector + 6
                elif sloc == 5:
                        regloc = (zloc-1)*nsector + 1
                else:
                        regloc = (zloc-1)*nsector + 2               

                #regloc = (zloc-1)*nsector + sloc        
                regid.InsertNextValue(regloc)
           #print x, y, z, regloc        
        elif -5.2 < z < -3.9083341: 
                apartition = np.linspace(-math.pi+math.pi/4, math.pi-math.pi/4, 4+0)
                zloc = np.argmax(zpartition>z)
        # Determine position in theta direction
                #print "zloc = ", zloc
                norm = np.linalg.norm([(x - midx), (y - midy)])
                angle = np.arctan2((y - midy)/norm, (x - midx)/norm)
                sloc = np.argmax(apartition>angle)

                if sloc == 0:
                        regloc = (zloc+1)*nsector + 2
                elif sloc == 1:
                        regloc = (zloc+1)*nsector + 3
                elif sloc == 2:
                        regloc = (zloc+1)*nsector + 4
                else:
                        regloc = (zloc+1)*nsector + 1
                      
                #print "sloc = ", regloc
                #regloc = (zloc+1)*nsector + sloc
                regid.InsertNextValue(regloc)
        else:
                regloc = 17
                regid.InsertNextValue(regloc)   
 
    regid.SetName("Regionid")
    #data.SetName("EP measurements")
    
    ugrid.GetCellData().AddArray(regid)
    #ugrid.GetCellData().AddArray(data)
    
    vtk_py.writeXMLUGrid(ugrid, "test.vtu")

 
isepiflip = True	
isendoflip = True

mesh = Mesh()
comm_common = mesh.mpi_comm() 

casename = "tomtec399019_refine_rot"
LVangle = 60

meshfilename = casename + ".hdf5"
f = HDF5File(comm_common, meshfilename, 'r')
f.read(mesh, casename, False)
f.close() 

#cmd = "gmsh"+" -3 ellipsoidal_baseline_mesh1.geo -o " + meshfilename
#os.system(cmd)

#ugrid = vtk_py.readUGrid(meshfilename)
ugrid = vtk_py.convertXMLMeshToUGrid(mesh)

xmlgrid = vtk_py.convertUGridToXMLMesh(ugrid)

Qelem = dolfin.FiniteElement("Quadrature", xmlgrid.ufl_cell(), degree=4, quad_scheme="default")
segFS = dolfin.FunctionSpace(xmlgrid, Qelem)
segid = dolfin.Function(segFS)
partitionmeshforEP(ugrid, nsector=6, nz=3)

cnt = 0
for cell in dolfin.cells(xmlgrid):

        idx = int(ugrid.GetCellData().GetArray("Regionid").GetTuple(cnt)[0])
        #print idx, cnt
        xmlgrid.domains().set_marker((cell.index(), idx), 3)
        cnt += 1
    
matid = dolfin.MeshFunction("size_t", xmlgrid, 3, xmlgrid.domains())

dolfin.File(casename+"_subDomain"+".pvd") << matid

xmlgrid, xmlfacet, xmledges = vtk_py.extractFeNiCsBiVFacet(ugrid, geometry="LV")

VQuadelem = dolfin.VectorElement("Quadrature", 
                              xmlgrid.ufl_cell(), 
                              degree=4, 
                              quad_scheme="default")
VQuadelem._quad_scheme = 'default'

fiberFS = dolfin.FunctionSpace(xmlgrid, VQuadelem)

ef, es, en, eC, eL, eR = vtk_py.addLVfiber(xmlgrid, fiberFS, casename, LVangle, -LVangle, [] , isepiflip, isendoflip)

f = dolfin.HDF5File(xmlgrid.mpi_comm(), casename+".hdf5", 'w')
f.write(xmlgrid, casename)
f.close()

f = dolfin.HDF5File(xmlgrid.mpi_comm(), casename+".hdf5", 'a') 
f.write(xmlfacet, casename+"/"+"facetboundaries") 
f.write(xmledges, casename+"/"+"edgeboundaries") 
f.write(ef, casename+"/"+"eF") 
f.write(es, casename+"/"+"eS") 
f.write(en, casename+"/"+"eN")
f.write(eC, casename+"/"+"eC") 
f.write(eL, casename+"/"+"eL") 
f.write(eR, casename+"/"+"eR") 
f.write(matid, casename+"/"+"matid")
f.close()

dolfin.File(casename+"facetboundaries.pvd")  << xmlfacet
dolfin.File(casename+"Edges.pvd")  << xmledges

'''
meshfilename = casename+".hdf5"
mesh = dolfin.Mesh()

f = dolfin.HDF5File(dolfin.mpi_comm_world(), meshfilename, 'r')
f.read(mesh, casename, False)

matid1 = dolfin.CellFunction('size_t', mesh)
f.read(matid1, casename+"/"+"matid")

print matid1[1598]
'''
'''for cellid in range(0,celldata.GetNumberOfTuples()):
		data = celldata.GetTuple(cellid)[0]
'''
'''Qelem2 = dolfin.FiniteElement("Quadrature", mesh.ufl_cell(), degree=4, quad_scheme="default")
FS = dolfin.FunctionSpace(mesh, Qelem2)
kappa = dolfin.Function(FS)
for cell in dolfin.cells(mesh):
        idx = cell.index() #int(ugrid.GetCellData().GetArray("Regionid").GetTuple(cnt)[0])
	print idx, matid1[idx], kappa.vector()[idx]
	if matid1[idx] ==1 or matid1[idx] == 6:
		kappa.vector()[idx] = 2/9.0
	else:
		kappa.vector()[idx] = 0.0 
	print idx, matid1[idx], kappa.vector()[idx]
	
        #print idx, cnt
        #xmlgrid.domains().set_marker((cell.index(), idx), 3)
        #cnt += 1
'''
