import sys

sys.path.append("/mnt/home/caichen3/lab/")
import os
import vtk_py as vtk_py
import dolfin as dolfin
isepiflip = True
isendoflip = True


casenumber=399019
casename = "tomtec" + str(casenumber)+"_refine"+"_rot"

LVangle = 60

if os.path.exists('fiber.vtu'):
   os.rename('fiber.vtu','fiber_coarse.vtu')

meshfilename = casename + ".vtk"

#cmd = "gmsh"+" -3 ellipsoidal_baseline.geo -o " + meshfilename
#os.system(cmd)

ugrid = vtk_py.readUGrid(meshfilename)

xmlgrid = vtk_py.convertUGridToXMLMesh(ugrid)

xmlgrid, xmlfacet, xmledges = vtk_py.extractFeNiCsBiVFacet(ugrid, geometry="LV")

VQuadelem = dolfin.VectorElement("Quadrature", 
                              xmlgrid.ufl_cell(), 
                              degree=4, 
                              quad_scheme="default")
VQuadelem._quad_scheme = 'default'

fiberFS = dolfin.FunctionSpace(xmlgrid, VQuadelem)

ef, es, en, eC, eL, eR = vtk_py.addLVfiber(xmlgrid, fiberFS, casename, LVangle, -LVangle, [] , isepiflip, isendoflip)
print("4")
f = dolfin.HDF5File(xmlgrid.mpi_comm(), casename+".hdf5", 'w')
f.write(xmlgrid, casename)
f.close()

f = dolfin.HDF5File(xmlgrid.mpi_comm(), casename+".hdf5", 'a') 
f.write(xmlfacet, casename+"/"+"facetboundaries") 
f.write(xmledges, casename+"/"+"edgeboundaries") 
print('3')
f.write(ef, casename+"/"+"eF") 
f.write(es, casename+"/"+"eS") 
f.write(en, casename+"/"+"eN")
f.write(eC, casename+"/"+"eC") 
f.write(eL, casename+"/"+"eL") 
f.write(eR, casename+"/"+"eR") 
f.close()
print('2')
dolfin.File(casename+"facetboundaries.pvd")  << xmlfacet
dolfin.File(casename+"Edges.pvd")  << xmledges
print('1')
if os.path.exists('fiber.vtu'):
   os.rename('fiber.vtu','fiber_refine.vtu')
