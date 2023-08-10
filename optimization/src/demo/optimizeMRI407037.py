import sys, pdb
from dolfin import * 
from dolfin_adjoint import *

#sys.path.append("/mnt/Research")
sys.path.append("/mnt/home/caichen3/lab")
from heArt_optimization.src.sim_protocols.optimization_o import optimization as optimization 

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
IODetails = {"casename" : "UnloadMesh",
             "directory_me" : '../data/MRI407037/mesh/',
             "PVloop_data_file" : '../data/MRI407037/PVloop_Base.npz',
             "Strain_data_file" : '../data/MRI407037/Strain_Base.npz',
             "outputfolder" : './outputs_optimize407037_MRI/',
             "Initfile": "Log407037_MRI.hdf5",
             } 

GuccioneParams = {"ParamsSpecified" : True, 
     		  "Passive model": {"Name": "Guccione"},
     		  "Passive params": {"Cparam": Constant(100*1.2),	#starting point in the end of nLoadSteps 
		                     "bff"  : Constant(29.0*0.4),
		                     "bfx"  : Constant(13.3*0.4),
                  		     "bxx"  : Constant(26.6*0.4),
				   },
	         "Active model": {"Name": "Scalar-contractility"},
     	         "Active params": {"tau" : 25, "t_trans" : 300, "B" : 4.75,  "t0" : 275,  "l0" : 1.08, \
				   "Tmax" : 0e3, "Ca0" : 4.35, "Ca0max" : 4.35, "lr" : 1.85},   #1.58 to 1
                 "HomogenousActivation": True,
		 "deg" : 4, 
                 "Kappa": 1e5,
                 "incompressible" : True,
                 }

SimDetails = {    
                  "GiccioneParams" : GuccioneParams, 
                  "nLoadSteps": 18, 
                  "UpperBd": 2000e3,
                  "LowerBd": 0.0,
                  "Nintermediate_stp": 30,#20,#10,#80,#70,#60,#50,#40,#30,#20,#10,#18, #8,#14, #num of intermedia steps
                  "Kreg": 0e-9,
                  "Kstr": 0,#3e2,#3e3,#3e1, #3e-2,#1e-2,	#try larger
                  "Kell": 0,#5e-1,#2,#8e-1,#1,#0,#1e-2,
                  "Kecc": 0,#5e-1,#2,#8e-1,#1,#0,#1e-1,#1e-2,
                  "Klvp": 1e4,#6e3,#1e5,#8e4,#4e4,#1e4,#7e3,#5e5,#5e4,#1e6, #1e5,	#try larger
                  "init_opt_val" : 2e2,#4e2,#6e2,#1e3,#4e3,#1e3,#2e2,#6e2,#8e3,#1e3,#4e2,#2e2,#5e1,#5e3,  	#shear max(changing distance of points in calculating)
                  "topid" : 4,
                  "LVendoid" : 2,
                  "RVendoid" : 0,
                  "epiid" : 1,
		  "abs_tol" : 1e-8,
		  "rel_tol" : 5e-7,
		  "max_newton_iter" : 40,
                 }

# Run Simulation
optimization(IODet=IODetails, SimDet=SimDetails)
#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
    
