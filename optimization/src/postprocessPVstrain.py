import os
import pandas as pd
import numpy as np
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
from matplotlib import pylab as plt
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator


res = 1 # Resolution of the interpolation (if 1: use the original data)
ms = 3 # Marker size for plotting
directory = "./data_PIG2829/" # Directory of the data


# Post process PV data ##########################################################
PVdata = pd.read_csv(directory +"RA.csv")
tpt = np.array(PVdata.values[:,0])
LVP = np.array(PVdata.values[:,1])
LVV = np.array(PVdata.values[:,2])

# Interpolate PV data
LVPfunc = interp1d(tpt - tpt[0], LVP, kind='linear', fill_value="extrapolate")
LVVfunc = interp1d(tpt - tpt[0], LVV, kind='linear', fill_value="extrapolate")
maxtpt = max(tpt) - tpt[0]
PV_tpt_arr = np.linspace(0, maxtpt, res*len(LVP))
LVP_fit_arr = LVPfunc(PV_tpt_arr)
LVV_fit_arr = LVVfunc(PV_tpt_arr)

print("Number of fitted data points = ", len(PV_tpt_arr))

plt.figure(1)
plt.plot(LVV, LVP, '-*', label="Expt", markersize=ms)
plt.plot(LVV_fit_arr, LVP_fit_arr, 'o', label="Fitted", markersize=ms)
plt.savefig(directory+"PVloop.png")

plt.figure(2)
plt.plot(tpt, LVP, '-*', label="Expt", markersize=ms)
plt.plot(PV_tpt_arr, LVP_fit_arr, 'o', label="Fitted", markersize=ms)
plt.savefig(directory+"Pwave.png")

plt.figure(3)
plt.plot(tpt, LVV, '-*', label="Expt", markersize=ms)
plt.plot(PV_tpt_arr, LVV_fit_arr, 'o', label="Fitted", markersize=ms)
plt.savefig(directory+"Vwave.png")

# Output interpolated PV loop to npz for optimization
np.savez(directory + "PVloop.npz", \
        LVP = LVP_fit_arr,\
        LVV = LVV_fit_arr)

#################################################################


# Post process strain data ##########################################################
Ecc = np.load(directory+"E_cc_RAP100.npz")["strain"]
Ecc_t = np.load(directory+"E_cc_RAP100.npz")["t"]

Ell = np.load(directory+"E_ll_RAP100.npz")["strain"]
Ell_t = np.load(directory+"E_ll_RAP100.npz")["t"]

# Determine ED time point at which strain = 0 #)
min_Ecc_tpt_arr = [np.argmin(abs(Ecc_)) for Ecc_ in Ecc]
min_Ell_tpt_arr = [np.argmin(abs(Ell_)) for Ell_ in Ell]
min_tpt_arr = np.concatenate((min_Ecc_tpt_arr, min_Ell_tpt_arr))
EDtpt = np.bincount(min_tpt_arr).argmax()

# Shift strain waveform so that Ecc and Ell = 0 at ED
Ecc_new = np.concatenate((Ecc[:,EDtpt:len(Ecc[0])], Ecc[:,0:EDtpt]), axis=1)
Ell_new = np.concatenate((Ell[:,EDtpt:len(Ell[0])], Ell[:,0:EDtpt]), axis=1)

# Interpolate strain data so that the number of time point is equal to that from the PV loop (i.e., len(LVP)
Eccfunc = interp1d(Ecc_t - Ecc_t[0], Ecc_new, kind='linear', fill_value="extrapolate")
Ellfunc = interp1d(Ell_t - Ell_t[0], Ell_new, kind='linear', fill_value="extrapolate")
maxtpt = max(Ecc_t) - Ecc_t[0]
E_tpt_arr = np.linspace(0, maxtpt, res*len(LVP))
Ecc_fit_arr = Eccfunc(E_tpt_arr)
Ell_fit_arr = Ellfunc(E_tpt_arr)



plt.figure(4)
fig, axs = plt.subplots(4, 4, sharex=True, sharey=True)
cnt = 0
for p in range(0, len(Ecc_new)):
    i = int(cnt/4)
    j = cnt%4
    cnt += 1
    axs[i,j].plot(Ecc_t, Ecc_new[p], '-')
    axs[i,j].plot(E_tpt_arr, Ecc_fit_arr[p], '*')
    axs[i,j].set_title("Sector = " + str(cnt), fontsize=8)
fig.tight_layout()
plt.savefig(directory+"Ecc.png")


plt.figure(5)
fig, axs = plt.subplots(4, 4, sharex=True, sharey=True)
cnt = 0
for p in range(0, len(Ell_new)):
    i = int(cnt/4)
    j = cnt%4
    cnt += 1
    axs[i,j].plot(Ell_t, Ell_new[p], '-')
    axs[i,j].plot(E_tpt_arr, Ell_fit_arr[p], '*')
    axs[i,j].set_title("Sector = " + str(cnt), fontsize=8)
fig.tight_layout()
plt.savefig(directory+"Ell.png")

# Output interpolated Strain waveforms to npz for optimization
np.savez(directory+"Strain.npz", \
        tpt = E_tpt_arr,\
        Ecc = Ecc_fit_arr,\
        Ell = Ell_fit_arr)

