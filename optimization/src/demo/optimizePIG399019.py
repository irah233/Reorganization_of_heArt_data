import sys, pdb
from dolfin import * 
from dolfin_adjoint import *

#sys.path.append("/mnt/Research")
sys.path.append("/mnt/home/caichen3/lab")
from heArt_optimization.src.sim_protocols.optimization3 import optimization as optimization 

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
IODetails = {"casename" : "UnloadMesh",
             "directory_me" : '../data/data_PIG399019/mesh/',
             "PVloop_data_file" : '../data/data_PIG399019/PVloop_Base.npz',
             "Strain_data_file" : '../data/data_PIG399019/Strain_Base.npz',
             "outputfolder" : './outputs_optimize399019/',
             "Initfile": "Log399019.hdf5",
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
				   "Tmax" : 0e3, "Ca0" : 4.35, "Ca0max" : 4.35, "lr" : 1.85}, #1.58 to 1.38
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
                  "Nintermediate_stp": 30,#20,#10,#20,#70,#60,#40,#20,#80,#60,#40,#20,#12,#18, #8,#14, #num of intermedia steps
                  "Kreg": 0e-9,
                  "Kstr": 0,#3e1,#0,#3e3,#3e1, #3e-2,#1e-2,	#try larger
                  "Kell": 0,#1,#8e-1,#1,#3,#3,#3,#2,#5,#8e-1,#0,#4e1,#8e1,#1e2,#1e1,#1,#8e-1,#1e-1,#0,#1e-1,
                  "Kecc": 0,#1,#8e-1,#0,#6e1,#8e1,#1e1,#1,#1e-1,
                  "Klvp": 3e4,#1e4,#6e3,#4e3,#4e3,#3e3,#1e3,#3e3,#1e5,#4e,	#try larger
                  "init_opt_val" : 1e4,#2e2,#1e4,#4e4,#1e3,#8e2,#4e4,#5e3,  	#shear max(changing distance of points in calculating)
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
    
