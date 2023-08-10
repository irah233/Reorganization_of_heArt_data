from dolfin import *

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

    eC0 = Function(fiberFS)
    eL0 = Function(fiberFS)
    eR0 = Function(fiberFS)
    
    f.read(f0, casename+"/"+"eF")
    f.read(s0, casename+"/"+"eS")
    f.read(n0, casename+"/"+"eN")

    f.read(eC0, casename+"/"+"eC")
    f.read(eL0, casename+"/"+"eL")
    f.read(eR0, casename+"/"+"eR")
		
    f0 = f0/sqrt(inner(f0, f0))
    s0 = s0/sqrt(inner(s0, s0))
    n0 = n0/sqrt(inner(n0, n0))

    eC0 = eC0/sqrt(inner(eC0, eC0))
    eL0 = eL0/sqrt(inner(eL0, eL0))
    eR0 = eR0/sqrt(inner(eR0, eR0))

    matid = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)
    f.read(matid, casename+"/"+"matid")

    return mesh, facetboundaries, edgeboundaries, f0, s0, n0, eC0, eL0, eR0, matid

directory = "./data/data_PIG393485/" # Directory of the data
IODetails = {"casename" : "tomtec393485_rot",
             "directory_me" : directory+'/mesh/',
             "directory_ep" : './', 
             "outputfolder" : '/mnt/Research/optimization_heart/src/demo/outputs_optimize',
             "folderName" : '',
             "caseID" : '',
             "isLV" : True} 

QUAD_DEG = 4
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["quadrature_degree"] = QUAD_DEG

# Display material ID, material directions and fiber directions
mesh, facetboundaries, edgeboundaries, f0, s0, n0, eC0, eL0, eR0, matid = readmesh(IODetails)
File(directory+"matid.pvd") << matid
File(directory+"eC0.pvd") << project(eC0, VectorFunctionSpace(mesh, "DG", 0))
File(directory+"eL0.pvd") << project(eL0, VectorFunctionSpace(mesh, "DG", 0))
File(directory+"eR0.pvd") << project(eR0, VectorFunctionSpace(mesh, "DG", 0))
File(directory+"ef.pvd") << project(f0, VectorFunctionSpace(mesh, "DG", 0))
File(directory+"es.pvd") << project(s0, VectorFunctionSpace(mesh, "DG", 0))
File(directory+"en.pvd") << project(n0, VectorFunctionSpace(mesh, "DG", 0))
