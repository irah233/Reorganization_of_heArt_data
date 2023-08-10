import os
import pandas as pd
import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
plt.rcParams["font.family"] = "Times New Roman"
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt

font1 = {'weight' : 'normal',
'fontname' : 'Times New Roman',
'size' : 24,}
ms=3

case = [396669,396670,399018,399019,407034,407036,407037]
for i in case:
    pig_ = i
    if pig_== 396669:
        total_t = 0.66925           
    if pig_== 396670:
        total_t = 0.4845
    if pig_== 399018:
        total_t = 0.5383
    if pig_== 399019:
        total_t = 0.5383
    if pig_== 407034:
        total_t = 0.5383
    if pig_== 407036:
        total_t = 0.48445
    if pig_== 407037:
        total_t = 0.5921
    path = input('Solve for ECHO or MRI: ')
    if path.lower().strip() == 'echo':
        path = str(pig_)+'echo'
    if path.lower().strip() == 'mri':
        path = str(pig_)+'mri'
    PVdata = pd.read_csv(path+'/PVLoop.csv')       
    PVbdata = pd.read_csv(path+'/'+str(pig_)+'.csv')
    t = np.array(PVbdata.values[:,0]) 
    V = np.array(PVbdata.values[:,2])
    P = np.array(PVbdata.values[:,1])                   
    LVP = np.array(PVdata.values[:,1])           
    LVV = np.array(PVdata.values[:,0])
    V_t = t
    tt=[]
    m=1
    preloading_point = (len(LVP)-len(t))
    while m < preloading_point+1 :
        t=np.append(t,((max(t)+1)*total_t))
        m+=1;
    res = 1
    LVPfunc = interp1d(t - t[0], LVP, kind='linear', fill_value="extrapolate")
    LVVfunc = interp1d(t - t[0], LVV, kind='linear', fill_value="extrapolate")
    maxtpt = max(t) - t[0]
    
    Vfunc = interp1d(V_t - V_t[0], V, kind='linear', fill_value="extrapolate")
    maxt = max(V_t) - V_t[0]
    V_tpt_arr = np.linspace(0, maxt, res*len(V))
    V_fit_arr = Vfunc(V_tpt_arr)
    
    PV_tpt_arr = np.linspace(0, maxtpt, res*len(V))
    LVP_fit_arr = LVPfunc(PV_tpt_arr)
    LVV_fit_arr = LVVfunc(PV_tpt_arr)
    
    plt.figure(1)
    plt.plot(LVV, LVP, '-', label="Simulation", markersize=ms, color='b', linewidth=2)
    plt.plot(V, P, 'o', label="Measurement", markersize=ms, color='darkorange', linewidth=2)
    plt.legend(loc='center',frameon=False, ncol=1, fontsize=24)
    plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})
    plt.xlabel(r"Volume (mL)", font1)
    plt.ylabel(r'Pressure (mmHg)', font1)
    plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
    plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
    plt.tight_layout()
    plt.savefig(path+"PVLoop"+".png")
    plt.close('all')
        
    
    sdata = pd.read_csv(path+'/total_Strain.csv')
    n=1
    t = np.array(sdata.values[:,0])
    ell=[]
    ecc=[]
    while n<33:
        ell.append(np.array(sdata.values[:,n]))
        n+=1
    cnt=0
    plt.figure(2)
    plt.figure(figsize=(14,9))
    for p in range(0, len(ell)//2):
        plt.subplot(4, 4, cnt+1)
        cnt += 1
        plt.plot(t, ell[2*p], '*', label = 'measurement')
        plt.plot(t, ell[2*p+1], '-', label = 'simulation')
        mini = min(ell[2*p])
        if min(ell[2*p+1]) < mini:
            mini = min(ell[2*p+1])
        maxi = max(ell[2*p])
        if max(ell[2*p+1]) > maxi:
            maxi = max(ell[2*p+1])
        #plt.set_title("Sector = " + str(cnt), fontsize=24)
        y_ticks = np.linspace(mini, maxi, 2) # two values for y-axis for the largest value
        x_ticks = np.linspace(min(t), max(t), 2) # two values for x-axis for the largest value
        plt.yticks(y_ticks)
        plt.xticks(x_ticks)
        plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
        plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
###############################################################################              
    plt.xlabel("Time (ms)",font1)
    plt.ylabel("Ell",font1)
    plt.legend(loc='upper right',  handletextpad = 0.05, borderaxespad= 5.5, frameon=False, ncol=2, numpoints =4,fontsize=24)
    plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})
    plt.tight_layout()
    plt.subplots_adjust(left=0.1,
                    bottom=0.2,
                    right=0.9,
                    top=0.9,
                    wspace=0.9,
                    hspace=0.9)
    plt.savefig(path+"Ell.png")
    plt.close('all')
        
    while n<65:
        ecc.append(np.array(sdata.values[:,n]))
        n+=1
    cnt=0
    plt.figure(3)
    plt.figure(figsize=(14,9))
    for p in range(0, len(ecc)//2):
        plt.subplot(4, 4, cnt+1)
        cnt += 1
        plt.plot(t, ecc[2*p], '*', label = 'measurement')
        plt.plot(t, ecc[2*p+1], '-', label = 'simulation')
        mini = min(ecc[2*p])
        if min(ecc[2*p+1]) < mini:
            mini = min(ecc[2*p+1])
        maxi = max(ecc[2*p])
        if max(ecc[2*p+1]) > maxi:
            maxi = max(ecc[2*p+1])
        #plt.set_title("Sector = " + str(cnt), fontsize=24)
        y_ticks = np.linspace(mini, maxi, 2) # two values for y-axis for the largest value
        x_ticks = np.linspace(min(t), max(t), 2) # two values for x-axis for the largest value
        plt.yticks(y_ticks)
        plt.xticks(x_ticks)
        plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
        plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
###############################################################################              
    plt.xlabel("Time (ms)",font1)
    plt.ylabel("Ecc",font1)
    plt.legend(loc='upper right',  handletextpad = 0.05, borderaxespad= 5.5, frameon=False, ncol=2, numpoints =4,fontsize=24)
    plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})
    plt.tight_layout()
    plt.subplots_adjust(left=0.1,
                    bottom=0.2,
                    right=0.9,
                    top=0.9,
                    wspace=0.9,
                    hspace=0.9)
    plt.savefig(path+"Ecc.png")
    plt.close('all')
