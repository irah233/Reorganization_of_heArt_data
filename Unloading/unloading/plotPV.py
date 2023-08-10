import csv
import glob
import numpy as np
import csv
import math
import os
import matplotlib 
matplotlib.use('Agg')
from matplotlib import pylab as plt

def extract_unloadPV(filename):

	reader = csv.reader(open(filename), delimiter=" ")
	tpt_array = []
	LVP_array = []
	LVV_array = []
	for row in reader:
		tpt_array.append(float(row[0]))
		LVP_array.append(float(row[1]))
		LVV_array.append(float(row[2]))
	
	
	tpt_array = np.array(tpt_array)
	LVP_array = np.array(LVP_array)
	LVV_array = np.array(LVV_array)
	
        tpt_new_array = []
        LVP_new_array = []
        LVV_new_array = []
        
        ncyc = int(np.max(tpt_array))
        for cyc in range(0, ncyc):
	    ind = np.where(tpt_array == cyc)
            tpt_new_array.append(tpt_array[ind])
            LVP_new_array.append(LVP_array[ind])
            LVV_new_array.append(LVV_array[ind])
	
	#return tpt_array[ind], LVP_array[ind], LVV_array[ind], Qmv_array[ind]
	return tpt_new_array, LVP_new_array, LVV_new_array

def readEchoPV(echoPVfilename):

    LVP = np.load(echoPVfilename)["LVP"]
    LVV = np.load(echoPVfilename)["LVV"]

    return LVP, LVV

directory = "./outputs_407037_MRI/tomtec407037/deformation_unloadED"
filename = directory+"/"+"BiV_unloadPV.txt"
echofilename = "407037_PVloop_MRI.npz"

tpt_arr, LVP_arr, LVV_arr = extract_unloadPV(filename)
LVP_exp, LVV_exp = readEchoPV(echofilename)

plt.figure(1)
cnt = 0
for (LVV, LVP) in zip(LVV_arr, LVP_arr):
    cnt += 1

    if(cnt == len(LVV_arr)):
        plt.plot(LVV, LVP, '-k', linewidth=3.0)

    else:
        plt.plot(LVV, LVP, '--', linewidth=3.0)


plt.plot(LVV_exp, LVP_exp, '*')
plt.savefig(directory+"unloadedPV_ECHO.png")
 
