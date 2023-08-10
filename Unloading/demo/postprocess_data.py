import sys, pdb
from dolfin import *

sys.path.append("/mnt/home/fanlei1")
from heArt.src.sim_protocols.run_BiV_ClosedLoop import run_BiV_ClosedLoop as run_BiV_ClosedLoop
from heArt.src.postprocessing.postprocessdataBiV2 import postprocessdata as postprocessdata

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
IODetails = {"casename" : "biv_idealized",
             "directory_me" : '../BiVMesh/',
             "directory_ep" : '../BiVMesh/',
             "outputfolder" : './outputs_BiVelectromechanics/',
             "folderName" : '',
             "caseID" : 'BiVelectromechanics'}

Circparam = {     "stop_iter" : 5,
                  };

SimDetails = {
                  "HeartBeatLength": 800.0,
                  "dt": 1.0,
                  "writeStep": 10.0,
                  #"GiccioneParams" : GuccioneParams,
                  "nLoadSteps": 5,
                  "DTI_EP": False,
                  "DTI_ME": False,
                  "d_iso": 1.5*0.01,
                  "d_ani_factor": 4.0,
                  "ploc": [[-0.083, 5.6,-1.16, 2.0, 1.0]],
                  "pacing_timing": [[4.0, 20.0]],
                  "closedloopparam": Circparam,
                  "Ischemia": False,
                  "isLV" : False,
                  "topid" : 4,
                  "LVendoid" : 2,
                  "RVendoid" : 3,
                  "epiid" : 1,
                  "abs_tol": 1e-9,
                  "rel_tol": 1e-9,
                 }



# Postprocessing
postprocessdata(IODet=IODetails, SimDet=SimDetails)

#  - - - - - - - - - - -
