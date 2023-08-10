import os
import pandas as pd
import numpy as np
import matplotlib as mpl
import statsmodels.api as sm
import math
import random
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
plt.rcParams["font.family"] = "Times New Roman"
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt

font3 = {'weight' : 'normal',
'fontname' : 'Times New Roman',
'size' : 24,}

font1 = {'weight' : 'normal',
'fontname' : 'Times New Roman',
'size' : 48,}

font2 = {'weight' : 'normal',
'fontname' : 'Times New Roman',
'size' : 64,}
all_color=['indianred','maroon','red','sienna','chocolate','goldenrod',
       'yellowgreen','forestgreen','lime','lightseagreen','cyan','deepskyblue',
       'navy','blueviolet','indigo','magenta']
path=''
ms=3
case = [396669,396670,399018,399019,407034,407036,407037]
echo_Plst=[];mri_Plst=[];echo_Plmt=[];mri_Plmt=[];
echo_Vlst=[];mri_Vlst=[];echo_Vlmt=[];mri_Vlmt=[];
ellE=[];ellM=[];ellES=[];ellMS=[]
eccE=[];eccM=[];eccES=[];eccMS=[]
time_m = []; time_e = []
Tmax_m=[];Tmax_e=[];all_time=[]
for i in case:
    pig_ = i
    if pig_== 396669:
        total_t = 0.66925           
    elif pig_== 396670:
        total_t = 0.4845
    elif pig_== 399018:
        total_t = 0.5383
    elif pig_== 399019:
        total_t = 0.5383
    elif pig_== 407034:
        total_t = 0.5383
    elif pig_== 407036:
        total_t = 0.48445
    elif pig_== 407037:
        total_t = 0.5921
    kind=['echo','mri']
    for na in kind:
        if len(path) > 4:
            lastl = path[-1]
            if lastl == 'o':
                path = 'echo'
            else:
                path = 'mri'
        if na.lower().strip() == 'echo':
            path = str(pig_)+na
        elif na.lower().strip() == 'mri':
            path = str(pig_)+na
        PVdata = pd.read_csv(path+'/PVLoop.csv')       
        PVbdata = pd.read_csv(path+'/'+str(pig_)+'.csv')
        t = np.array(PVbdata.values[:,0][:]) 
        V = np.array(PVbdata.values[:,2][:])
        P = np.array(PVbdata.values[:,1][:])                   
        LVP = np.array(PVdata.values[:,1][17:])           
        LVV = np.array(PVdata.values[:,0][17:])
        V_t = t
        
        m=1
        preloading_point = (len(LVP)-len(t))
        while m < preloading_point+1 :
            t=np.append(t,((max(t))*total_t))
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
        
        if 'echo' in path:
            echo_Plst.append(LVP_fit_arr)
            echo_Vlst.append(LVV_fit_arr)
            echo_Plmt.append(P)
            echo_Vlmt.append(V)
            time_e.append(t)
        else:
            mri_Plst.append(LVP_fit_arr)
            mri_Vlst.append(LVV_fit_arr)
            mri_Plmt.append(P)
            mri_Vlmt.append(V)
            '''if '396670' in path:
                t = np.linspace(0,0.6692499999999981,29)
            elif '399018' in  path:
                t = np.linspace(0,0.6692499999999981,29)
            elif'399019' in path:'''
            t = np.linspace(0,0.6692499999999981,31)
                #0.6692499999999981
            time_m.append(t)
        sdata = pd.read_csv(path+'/total_Strain.csv')
        n=1
        t = np.array(sdata.values[:,0])
        while n<65:
            if n<33:
                if n%2==0:
                    if 'echo' in path:
                        ellES.extend(np.array(sdata.values[:,n]))
                    else:
                        ellMS.extend(np.array(sdata.values[:,n]))
                else:
                    if 'echo' in path:
                        ellE.extend(np.array(sdata.values[:,n]))
                    else:
                        ellM.extend(np.array(sdata.values[:,n]))
            else:
                if n%2==0:
                    if 'echo' in path:
                        eccES.extend(np.array(sdata.values[:,n]))
                    else:
                        eccMS.extend(np.array(sdata.values[:,n]))
                else:
                    if 'echo' in path:
                        eccE.extend(np.array(sdata.values[:,n]))
                    else:
                        eccM.extend(np.array(sdata.values[:,n]))
            n+=1
    csre=pd.read_csv(str(i)+'echo/'+'Tmax.csv')
    csrm=pd.read_csv(str(i)+'mri/'+'Tmax.csv')
    time = np.array(csre.values[:,0])
    ECHO=np.array(csre.values[:,1])
    MRI=np.array(csrm.values[:,1])
    
    all_time.append(time*total_t/len(time))
    Tmax_m.append(MRI/1000)
    Tmax_e.append(ECHO/1000)

fig = plt.figure(figsize=(36,12))
gs = fig.add_gridspec(1,7, wspace=0)
axs = gs.subplots(sharey=True)
cnt=0
for p in range(len(mri_Plst)):
    cnt += 1
    axs[p].plot(all_time[p], Tmax_e[p], '-', lw=8, label = 'ECHO',color='blue')
    axs[p].plot(all_time[p], Tmax_m[p], '-', lw=8, label = 'MRI',color='darkorange')
    mini = min(echo_Plst[p])
    if min(mri_Plst[p]) < mini:
        mini = min(mri_Plst[p])
    maxi = max(echo_Plst[p])
    if max(mri_Plst[p]) > maxi:
        maxi = max(mri_Plst[p])
    y_ticks = np.linspace(0, 40, 6) # two values for y-axis for the largest value
    x_ticks = np.linspace(0.1, 0.5, 2)
    #plt.yticks(y_ticks)
    plt.yticks(fontsize=48, rotation=0, fontname="Times New Roman")
    axs[p].set_xlim([0, 0.6])
    axs[p].set_xticks(np.linspace(0.1, 0.5, 2))
    axs[p].set_xticklabels(x_ticks, fontsize = 48)
##############################################################################
axs[0].set_yticklabels(y_ticks, fontsize = 48)
plt.legend(bbox_to_anchor=(-1, -0.025),  frameon=False, ncol=4, numpoints =4, fontsize=48)
plt.xlabel("Time (s)",font1, loc = 'right')
axs[0].set_ylabel("Contractility (kPa)",font1)
#plt.legend(loc='lower right',  handletextpad = 0.05, borderaxespad= 5.5, frameon=False, ncol=2, numpoints =4,fontsize=24)
plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})
plt.tight_layout()
plt.savefig("Tmax.png")
plt.close('all')
      
fig = plt.figure(figsize=(50,14))
gs = fig.add_gridspec(1,7, wspace=0)
axs = gs.subplots(sharey=True)
cnt=0
for p in range(len(mri_Plst)):
    cnt += 1
    axs[p].plot(time_e[p], echo_Plst[p], '-', lw=8,label = 'ECHO Model',color='blue')
    axs[p].plot(time_m[p], mri_Plst[p], '-', lw=8,label = 'MRI Model',color='darkorange')
    #axs[p].plot(time_e[p], echo_Plmt[p], '*', label = 'ECHO Expe.',color='blue')
    axs[p].plot(time_m[p], mri_Plmt[p], '*', markersize=20,label = 'Expt.',color='blue')
    mini = min(echo_Plst[p])
    if min(mri_Plst[p]) < mini:
        mini = min(mri_Plst[p])
    maxi = max(echo_Plst[p])
    if max(mri_Plst[p]) > maxi:
        maxi = max(mri_Plst[p])
    y_ticks = np.linspace(0, 100, 5) # two values for y-axis for the largest value
    x_ticks = np.linspace(0.1, 0.5, 2)
    plt.yticks(y_ticks)
    plt.yticks(fontsize=48, rotation=0, fontname="Times New Roman")
    axs[p].set_xlim([0, 0.6])
    axs[p].set_xticks(np.linspace(0.1, 0.5, 2))
    axs[p].set_xticklabels(x_ticks, fontsize = 48)
##############################################################################
axs[0].set_yticklabels(y_ticks, fontsize = 48)
plt.legend(bbox_to_anchor=(-0.5, -0.05),  frameon=False, ncol=4, numpoints =4, fontsize=48)
plt.xlabel("Time (s)",font1, loc = 'right')
axs[0].set_ylabel("Pressure (mmHg)",font1)
#plt.legend(loc='lower right',  handletextpad = 0.05, borderaxespad= 5.5, frameon=False, ncol=2, numpoints =4,fontsize=24)
plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})
plt.tight_layout()
plt.savefig("Pressure_sub.png")
plt.close('all')
    
fig = plt.figure(figsize=(72,16))
gs = fig.add_gridspec(1,7, wspace=0)
axs = gs.subplots(sharey=True)
cnt=0
for p in range(len(mri_Plst)):
    cnt += 1
    axs[p].plot(echo_Vlst[p], echo_Plst[p], '-', lw=8, label = 'ECHO Model',color='blue')
    axs[p].plot(mri_Vlst[p], mri_Plst[p], '-', lw=8, label = 'MRI Model',color='darkorange')
    axs[p].plot(echo_Vlmt[p], echo_Plmt[p], 'o', markersize=20, label = 'ECHO Expt.',color='blue')
    axs[p].plot(mri_Vlmt[p], mri_Plmt[p], '>', markersize=20, label = 'MRI Expt.',color='darkorange')
    mini = min(echo_Plst[p])
    if min(mri_Plst[p]) < mini:
        mini = min(mri_Plst[p])
    maxi = max(echo_Plst[p])
    if max(mri_Plst[p]) > maxi:
        maxi = max(mri_Plst[p])
    y_ticks = np.linspace(0, 100, 5) # two values for y-axis for the largest value
    x_ticks = np.linspace(40, 120, 2)
    plt.yticks(y_ticks)
    plt.yticks(fontsize=64, rotation=0, fontname="Times New Roman")
    axs[p].set_xlim([0, 180])
    axs[p].set_xticks(np.linspace(40, 120, 2))
    axs[p].set_xticklabels(x_ticks, fontsize = 64)
##############################################################################
axs[0].set_yticklabels(y_ticks, fontsize = 64)
plt.legend(bbox_to_anchor=(0.2, -0.025),  frameon=False, ncol=4, numpoints =4, fontsize=64)
plt.xlabel("Volume (mL)",font2, loc = 'right')
axs[0].set_ylabel("Pressure (mmHg)",font2)
#plt.legend(loc='lower right',  handletextpad = 0.05, borderaxespad= 5.5, frameon=False, ncol=2, numpoints =4,fontsize=24)
plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})
plt.tight_layout()
plt.savefig("PV_sub.png")
plt.close('all')   

def func(x, a):
    return a * x
  
def oneline(echo,mri,type_,unit,name_p):
    filename1=[]
    max_x = max(echo)
    max_y = max(mri)
    min_x = min(echo)
    min_y = min(mri)
    maxrange = max_x
    minrange = min_x
    if max_x < max_y:
        maxrange = max_y
    if min_x > min_y:
        minrange = min_y
    filename1.append(type_)
    color1 = random.choice(all_color)
    plt.figure(figsize=(8,6))
    res = sm.OLS((mri), (echo)).fit()
    x = np.linspace(minrange, maxrange, 10)          # process identical data point for pressure, volume, and strain
    plt.plot(x, func(x,res.params),'-', label="MRI = %.3fECHO" % res.params[0], color='k')
    i=0
    for i in range(0,7):
        plt.plot(echo[len(echo)//7*i:len(echo)//7*(i+1)-1], mri[len(echo)//7*i:len(echo)//7*(i+1)-1],'o', label='case '+str(i+1), markersize=2, markeredgewidth=2 ,color=all_color[i])
    plt.text(minrange+(maxrange-minrange)*3/5, maxrange-minrange*-6/6, "R = %.3f" % float(math.sqrt(res.rsquared)),fontdict = font3)
    #plt.text(minrange+(maxrange-minrange)/6, maxrange-minrange*2/7, "R^2 = %.3f" % float((res.rsquared)),fontdict = font1)
    plt.legend(loc='upper left', frameon=False,ncol=2, numpoints=3,\
               handletextpad=0.05,prop={'family' : 'Times New Roman', 'size' : '24'})
    plt.ylabel(name_p+r' MRI ('+unit+')', font3)
    plt.xlabel(name_p+r' ECHO ('+unit+')', font3)
    y_ticks = np.linspace(minrange, maxrange+60, 3) # three values for y-axis for the largest value
    x_ticks = np.linspace(minrange, maxrange, 3) # three values for x-axis for the largest value
    plt.yticks(y_ticks)
    plt.xticks(x_ticks)
    plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
    plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
    plt.tight_layout()
    plt.savefig(filename1[0])
    plt.close('all')
expell=[]
simll=[]
expecc=[]
simcc=[]
oneline((ellE),(ellES),'Ell ECHO','%','Ell ')
oneline((eccE),(eccES),'Ecc ECHO','%','Ecc ')
oneline((ellM),(ellMS),'Ell MRI','%','Ell ')
oneline((eccM),(eccMS),'Ecc MRI','%','Ecc ')    
oneline((ellES),(ellMS),'Ell','%','Ell')
oneline((eccES),(eccMS),'Ecc','%','Ecc ')
expell.extend(ellE)
expell.extend(ellM)
simll.extend(ellES)
simll.extend(ellMS)
oneline((expell),(simll),'Ell ALL','%','$E_{ll} $')
expecc.extend(eccE)
expecc.extend(eccM)
simcc.extend(eccES)
simcc.extend(eccMS)
oneline((expecc),(simcc),'Ecc ALL','%','$E_{cc} $')


