import sys, pdb
from dolfin import * 
from dolfin_adjoint import *

#sys.path.append("/mnt/Research")
sys.path.append("/mnt/home/caichen3/lab")
from heArt_optimization.src.sim_protocols.optimization3 import optimization as optimization 

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
IODetails = {"casename" : "tomtec393485_refine_rot",
             "directory_me" : '../data/data_PIG393485/mesh/',
             "PVloop_data_file" : '../data/data_PIG393485/PVloop_Base.npz',
             "Strain_data_file" : '../data/data_PIG393485/Strain_Base.npz',
             "outputfolder" : './outputs_optimize393485/',
             #"Initfile": "Log393485.hdf5",
             } 

GuccioneParams = {"ParamsSpecified" : True, 
     		  "Passive model": {"Name": "Guccione"},
     		  "Passive params": {"Cparam": Constant(119*8.7),	#starting point in the end of nLoadSteps 
		                     "bff"  : Constant(29.0*1.5),
		                     "bfx"  : Constant(13.3*1.5),
                  		     "bxx"  : Constant(26.6*1.5),
				   },
	         "Active model": {"Name": "Scalar-contractility"},
     	         "Active params": {"tau" : 25, "t_trans" : 300, "B" : 4.75,  "t0" : 275,  "l0" : 1.58, \
				   "Tmax" : 0e3, "Ca0" : 4.35, "Ca0max" : 4.35, "lr" : 1.85},
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
                  "Nintermediate_stp": 40,#20,#12,#18, #8,#14, #num of intermedia steps
                  "Kreg": 0e-9,
                  "Kstr": 0,#3e2,#3e3,#3e1, #3e-2,#1e-2,	#try larger
		  "Kell": 0,
		  "Kecc": 0,
                  "Klvp": 6e7,#5e5,#5e4,#5e2,#5e6,#1e6, #1e5,	#try larger
                  "init_opt_val" : 2e1,#4e4,#2e3,#2e1,#5e1,#5e3,  	#Tmax(changing distance of points in calculating)
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
    
