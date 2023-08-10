import os
import pandas as pd
import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt

font1 = {'weight' : 'normal',
'fontname' : 'Times New Roman',
'size' : 24,}
ms=3
#path = input('Solve for ECHO or MRI: ')
case = [396669,396670,399018,399019,407034,407036,407037]#396669,396670,399018,399019,407034,407036,407037
Pe1=[];Ve1=[];filename1=[];Pe2=[];Ve2=[]
echo_ellt_error=[];echo_ecct_error=[];mri_ellt_error=[];mri_ecct_error=[]
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
    kind=['echo','mri']#'echo','mri'
    for na in kind:
        '''if len(path) > 4:
            lastl = path[-1]
            if lastl == 'o':
                path = 'echo'
            else:
                path = 'mri'''
        if na.lower().strip() == 'echo':
            path = str(pig_)+'echo'
        elif na.lower().strip() == 'mri':
            path = str(pig_)+'mri'
        PVdata = pd.read_csv(path+'/PVLoop.csv')       
        PVbdata = pd.read_csv(path+'/'+str(pig_)+'.csv')
        t = np.array(PVbdata.values[:,0][:30]) 
        V = np.array(PVbdata.values[:,2][:30])
        P = np.array(PVbdata.values[:,1][:30])                   
        LVP = np.array(PVdata.values[:,1][18:])           
        LVV = np.array(PVdata.values[:,0][18:])
        V_t = t
        Vsum=Psum=0
        for i in range(len(LVV)):
            Vsum += abs(V[i] - LVV[i])/V[i]*100
            Psum += abs(P[i] - LVP[i])/abs(P[i])*100
        Verror = Vsum/len(LVV)
        Perror = Psum/len(LVP)
        if na.lower().strip() == 'echo':
            Pe1.append(Perror)
            Ve1.append(Verror*1e10)
        if na.lower().strip() == 'mri':
            Pe2.append(Perror)
            Ve2.append(Verror*1e10)
        sdata = pd.read_csv(path+'/'+'total_Strain.csv')
        data_1=np.array(sdata.values)
        m=1; sum_error=0; ell_error=[]; ecc_error=[]
        while m != 33:
            S = np.array(sdata.values[:,m])
            num = m-1
            if (m%2)==0:
                S1 =  np.array(sdata.values[:,num])
                n=0;maxi=[]
                while n<len(S1):
                    maxi.append(max(abs(S1)))
                    n+=1
                maxi =np.array(maxi)
                SS = abs(S-S1)/abs(maxi)*100
                for ele in SS:
                    sum_error += ele
                sum_error = sum_error/len(SS)
                ell_error.append(sum_error)
                sum_error=0
            m+=1
        while m != 65:
            S = np.array(sdata.values[:,m])
            num = m-1
            if (m%2)==0:
                S1 =  np.array(sdata.values[:,num])
                n=0;maxi=[]
                while n<len(S1):
                    maxi.append(max(abs(S1)))
                    n+=1
                maxi =np.array(maxi)
                SS = abs(S-S1)/abs(maxi)*100
                for ele in SS:
                    sum_error += ele
                sum_error = sum_error/len(SS)
                ecc_error.append(sum_error)
                sum_error=0
            m+=1
        if na.lower().strip() == 'echo':
            echo_ellt_error.append(ell_error)
            echo_ecct_error.append(ecc_error)
        if na.lower().strip() == 'mri':
            mri_ellt_error.append(ell_error)
            mri_ecct_error.append(ecc_error)
        clst=[]
        filename1=[]
        part = []; err = []
        color=['indianred','maroon','red','sienna','chocolate','goldenrod',
               'yellowgreen','forestgreen','lime','lightseagreen','cyan','deepskyblue',
               'navy','blueviolet','indigo','magenta']
        m=0
        while (m<(16)):
            part.append(m+1)
            m += 1
        filename1.append('Ell Error')
        plt.figure(figsize=(7,5))
        plt.bar(part,ell_error,color=color)
        plt.xlabel(r"AHA segment (#)", font1)
        plt.ylabel(r'Ell error (%)', font1)
        x_ticks = np.linspace(0, 16, 5)  # 5 values for x-axis between the largest value to 0;
        plt.xticks(x_ticks)
        plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
        y_ticks = np.linspace(0, 150, 3)  # 3 values for x-axis between the largest value to 0;
        plt.yticks(y_ticks)
        plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
        plt.tight_layout()
        plt.savefig(path+filename1[0])
        plt.close('all')
        part = []; err = []
        m=0
        while  (m<(16)):
            part.append(m+1)
            m += 1
        filename1.append('Ecc Error')
        plt.figure(figsize=(7,5))
        plt.bar(part,ecc_error,color=color)
        plt.xlabel(r"AHA segment (#)", font1)
        plt.ylabel(r'Ecc error (%)', font1)
        x_ticks = np.linspace(0, 16, 5)  
        plt.xticks(x_ticks)
        plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
        y_ticks = np.linspace(0, 50, 3)  
        plt.yticks(y_ticks)
        plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
        plt.tight_layout()
        plt.savefig(path+filename1[1])
        plt.close('all')
        part = ['Pressure','Volume']; err = []  
filename1=[]
x_axis = np.arange(len(case))
max_x = max(Pe1)
max_y = max(Pe2)
maxrange = max(max_x,max_y)
filename1.append('Pressure Abs Error')
plt.figure(figsize=(8,6))
plt.bar(x_axis-0.2, Pe1, width = 0.4, label = 'ECHO', color='blue')
plt.bar(x_axis+0.2, Pe2, width = 0.4, label = 'MRI', color='darkorange')
plt.legend(loc='lower center', borderaxespad= -4, frameon=False, ncol=2, numpoints=4,\
           handletextpad=0.05,prop={'family' : 'Times New Roman', 'size' : '24'})
plt.ylabel(r'Pressure (%)', font1)
y_ticks = np.linspace(0, maxrange, 3)
plt.yticks(y_ticks)
plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.tight_layout()
plt.savefig(filename1[0])
plt.close('all')

maxrange = max(max(Ve1),max(Ve2))
filename1.append('Volume Abs Error')
evalue= 10**10
plt.figure(figsize=(8,6))
plt.bar(x_axis-0.2, Ve1, width = 0.4, label = 'ECHO', color='blue')
plt.bar(x_axis+0.2, Ve2, width = 0.4, label = 'MRI', color='darkorange')
plt.legend(loc='lower center', borderaxespad= -4, frameon=False, ncol=2, numpoints=4,\
           handletextpad=0.05,prop={'family' : 'Times New Roman', 'size' : '24'})
plt.ylabel(r'Volume (%)(10^(-10))', font1)
y_ticks = np.linspace(0, maxrange, 3)
plt.yticks(y_ticks)
plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.tight_layout()
plt.savefig(filename1[1])
plt.close('all')