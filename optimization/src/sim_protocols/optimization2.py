from dolfin import *
from dolfin_adjoint import *
import numpy as np
import os
from pyadjoint.placeholder import Placeholder
#from ..optimization.MyReducedFunctional3 import MyReducedFunctional
#from pyadjoint.reduced_functional_numpy import ReducedFunctionalNumPy 
#from ..sim_protocols.heArtsolver import heArtsolver
from ..sim_protocols.heArtsolver_clean2 import heArtsolver_clean
from collections import OrderedDict
from mpi4py import MPI as pyMPI
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
from matplotlib import pylab as plt

def GetTmax(Tmax_):

    comm = MPI.comm_world
    rank = comm.Get_rank()

    try:
        val_local = Tmax_.vector()[0]
    except IndexError:
        val_local = 0.0

    Tmax = MPI.sum(comm, val_local)

    return Tmax


def printTmax(Tmax_array):

    comm = MPI.comm_world
    rank = comm.Get_rank()

    if(rank ==0):
        print("Tmax = ")

    for m1 in Tmax_array:

        Tmax = GetTmax(m1)
        
        if(rank == 0):
            print("m = %15.10f" % (Tmax), flush=True)



def eval_cb(j, m):

    comm = MPI.comm_world
    rank = comm.Get_rank()

    print("j = %10.7f" %(j), flush=True)
    printTmax(m)

def readmesh(parameters):

    comm = MPI.comm_world
    rank = comm.Get_rank()

    mesh = Mesh() 
    directory = parameters["directory_me"]
    casename = parameters["casename"]
    outputfolder = parameters["outputfolder"]
    folderName = parameters["folderName"]
    
    meshfilename = directory + casename + ".hdf5"
    f = HDF5File(MPI.comm_world, meshfilename, 'r')
    f.read(mesh, casename, False)
    
    facetboundaries = MeshFunction("size_t", mesh, 2)
    f.read(facetboundaries, casename+"/"+"facetboundaries")
    
    edgeboundaries = MeshFunction("size_t", mesh, 1)
    f.read(edgeboundaries, casename+"/"+"edgeboundaries")
    
    deg = 4
    VQuadelem = VectorElement("Quadrature", 
                              mesh.ufl_cell(), 
                              degree=deg, 
                              quad_scheme="default")
    VQuadelem._quad_scheme = 'default'
    
    fiberFS = FunctionSpace(mesh, VQuadelem)
    
    f0 = Function(fiberFS)
    s0 = Function(fiberFS)
    n0 = Function(fiberFS)
    
    f.read(f0, casename+"/"+"eF")
    f.read(s0, casename+"/"+"eS")
    f.read(n0, casename+"/"+"eN")
		
    f0 = f0/sqrt(inner(f0, f0))
    s0 = s0/sqrt(inner(s0, s0))
    n0 = n0/sqrt(inner(n0, n0))



    return mesh, facetboundaries, edgeboundaries, f0, s0, n0

def optimization(IODet, SimDet):

    tape = Tape()
    set_working_tape(tape)

    mesh, facetboundaries, edgeboundaries, f0, s0, n0 = readmesh(IODet)

    PVloopdatafilename = IODet["PVloop_data_file"]
    LVPtargets = np.load(PVloopdatafilename)["LVP"]#[0:10]
    LVVtargets = np.load(PVloopdatafilename)["LVV"]#[0:10]
    ndatapts = len(LVPtargets)

    Rspace = FunctionSpace(mesh, "R", 0)
    Tmax_ctrls = OrderedDict()
    for p in range(0,ndatapts):
        Tmax_ctrls[p] = Function(Rspace)

    Tmax_bds = [[0.0]*ndatapts, [500e3]*ndatapts]

    w, J, LVP_array, LVV_array = heArtsolver_clean(Tmax_ctrls, LVVtargets, LVPtargets, mesh, facetboundaries, edgeboundaries, f0, s0, n0, SimDet, IODet)

    m1 = [Control(Tmax_ctrl) for Tmax_ctrl in Tmax_ctrls.values()]
    rf = ReducedFunctional(J, m1, eval_cb_post = eval_cb)
    f_opt = minimize(rf, bounds = Tmax_bds, options={"maxiter": 50})

    printTmax(f_opt)

    #[print("m = %15.10f" % (float(m1.vector()[0]))) for m1 in f_opt]

    w, J, LVP_array, LVV_array = heArtsolver_clean(f_opt, LVVtargets, LVPtargets, mesh, facetboundaries, edgeboundaries, f0, s0, n0, SimDet, IODet, isannotate=False)

    print(LVP_array)
    print(LVV_array)


    plt.figure(2)
    Tmax_array = []
    for T in f_opt:
        Tmax_ = GetTmax(T)
        Tmax_array.append(Tmax_)
    plt.plot(np.arange(0,len(Tmax_array)), Tmax_array)
    plt.xlabel("Time point")
    plt.ylabel("Tmax (Pa)")
    plt.savefig("Tmax2.png")

    plt.figure(1)
    plt.plot(LVV_array, LVP_array)
    plt.plot(LVVtargets, LVPtargets, '*')
    plt.xlabel("LVV (ml)")
    plt.ylabel("LVP (mmHg)")
    plt.savefig("PVloop2.png")
 


