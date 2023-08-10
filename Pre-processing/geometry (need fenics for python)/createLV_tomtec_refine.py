import sys
sys.path.append("/mnt/home/caichen3/lab/")

import vtk_py as vtk_py
import os
import numpy as np
import vtk
from pyquaternion import Quaternion

x_data=-0.963
y_data=0.187
z_data=0.194

Laxis = np.array([x_data,y_data,z_data])


refinemesh=0.8

zoffset = 0.4 #0.38

angle = np.arccos(Laxis[2])


raxis = np.cross(Laxis, np.array([0,0,1]))
raxis = raxis / np.linalg.norm(raxis)




print (angle)
volumedatadir = './'#+pigname+'/'+week+'/Volume_data/'
casenumber=399019
endostr = str(casenumber)+"endo"
epistr = str(casenumber)+"epi"
#endofiles_tp = sorted([int(f[len(f)-6:len(f)-4]) for f in os.listdir(volumedatadir) if (f[len(f)-4:len(f)] == '.stl') and (endostr in f)])

#dataname = pigname+"_Aortic_Banding_"+week+"_LV"
#thicknessfilename = volumedatadir + dataname + "_thickness.txt"




endofilename = volumedatadir + endostr +  '.stl'
endopdata = vtk_py.readSTL(endofilename)

epifilename = volumedatadir + epistr + '.stl'
epipdata = vtk_py.readSTL(epifilename)

rotatedepipdata = vtk_py.rotatePData_w_axis(epipdata, angle, raxis)
rotatedendopdata = vtk_py.rotatePData_w_axis(endopdata, angle, raxis) #untitled

height =  min(rotatedendopdata.GetBounds()[5], rotatedepipdata.GetBounds()[5]) - zoffset

clipped_endo, clipped_epi = vtk_py.clipSurfacesForCutLVMesh(rotatedendopdata, rotatedepipdata, height, verbose=True)



vtk_py.writeSTL(clipped_epi, "clipped_epi_tmp.stl")
vtk_py.writeSTL(clipped_endo, "clipped_endo_tmp.stl")



filename = './tomtec'+str(casenumber)+'_refine'

vtk_py.createLVmesh(filename, refinemesh, "clipped_epi_tmp.stl", "clipped_endo_tmp.stl")


os.remove("clipped_epi_tmp.stl")
os.remove("clipped_endo_tmp.stl")

os.remove("LVtemp.geo.bak")

ugrid = vtk_py.readUGrid(filename +  '.vtk')
rotmat = vtk.vtkMatrix4x4()
rotmat.SetElement(0,0,1)
rotmat.SetElement(1,1,1)
rotmat.SetElement(2,2,1)
rotmat.SetElement(2,3,-ugrid.GetBounds()[5]*1)
rotmat.SetElement(3,3,1)

#rotmat = np.array([[0.1, 0, 0, 0], [0, 0.1, 0, 0], [0, 0, 0.1, 0], [0, 0, 0, 1]])
ugrid_transform = vtk_py.transform_mesh_w_4x4mat(filename + '.vtk', filename + '_rot.vtk', rotmat)
#ugrid_transform = vtk_py.readUGrid(filename +  '_rot.vtk')
print ("ugrid=trans", ugrid_transform)

