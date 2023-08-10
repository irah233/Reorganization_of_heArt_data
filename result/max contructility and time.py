import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

font1 = {'weight' : 'normal',
'fontname' : 'Times New Roman',
'size' : 24,}

def cal_mean_and_diff(datalist):
    mean = []
    diff = []
    for n in range(0,16):
        lst=[]
        for m in range(0,6):
            if datalist[m][n] < 400 and datalist[m][n] > 100:
                lst.append(datalist[m][n]/0.0075/1000)
        mean1 = np.mean(np.array(lst))
        diff1 = np.std(np.array(lst))
        mean.append(mean1)
        diff.append(diff1)
    return mean, diff        
            

case = [396669,396670,399018,399019,407034,407036,407037]
totalecho=[];totalmri=[];ttecho=[];ttmri=[]
echowall=[];mriwall=[]
for swine in case:
    echol = []
    mril = []
    echot = []
    mrit = []
    pig_ = swine
    if pig_ == 396669:
            total_t = 0.66925           
    elif pig_ == 396670:
            total_t = 0.4845
    elif pig_ == 399018:
            total_t = 0.5383
    elif pig_ == 399019:
            total_t = 0.5383
    elif pig_ == 407034:
            total_t = 0.5383
    elif pig_ == 407036:
            total_t = 0.48445
    elif pig_ == 407037:
            total_t = 0.5921
    kind=['echo','mri']
    for na in kind:
        '''file =open(str(swine)+na+'/regional_Tmax.csv')
        cscr= csv.reader(file)
        n=0
        while n < 16:
            print(len(cscr[0]))
            n+=1'''
        csr=pd.read_csv(str(swine)+na+'/regional_Tmax.csv')
        m=0
        res=[]
        while m < 16:
            findmax = []
            value = np.array(csr.values[:,m][:30])
            value1 = value.tolist()
            maxf = max(value)
            index = ((value1.index(maxf))+1)*total_t*1000/30
            findmax.append(index)
            findmax.append(maxf)
            if na == 'echo':
                echol.append((maxf))
                echot.append((index))
            if na == 'mri':
                mril.append((maxf))
                mrit.append((index))
            m+=1
            res.append(m)
        if na == 'echo':
            totalecho.append(np.array(echol))
            ttecho.append(np.array(echot))
        if na == 'mri':
                totalmri.append(np.array(mril))
                ttmri.append(np.array(mrit))
    filename1 = []
    echo = np.array(echot)
    mri = np.array(mrit)
    res = np.array(res)
    max_x = max(echot)
    max_y = max(mrit)
    echotl=echo.tolist()
    mritl=mri.tolist()
    maxrange = max(max_x,max_y)
    filename1.append(str(swine)+'Time to Peak')
    plt.figure(figsize=(8,6))
    plt.bar(res-0.2, echo/0.0075, width = 0.4, label = 'ECHO', color='blue')
    plt.bar(res+0.2, mri/0.0075, width = 0.4, label = 'MRI', color='darkorange')
    plt.legend(loc='lower center', borderaxespad= -4, frameon=False, ncol=2, numpoints=4,\
        handletextpad=0.05,prop={'family' : 'Times New Roman', 'size' : '24'})
    plt.ylabel(r'Time (ms)', font1)
    y_ticks = np.linspace(0, maxrange, 3)
    plt.yticks(y_ticks)
    plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
    plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
    #plt.title("Percentage Error of Volume", fontsize=24, fontname="Times New Roman")
    plt.tight_layout()
    plt.savefig(filename1[0])
    plt.close('all')

mecho, stdecho=cal_mean_and_diff(totalecho)
mmri, stdmri=cal_mean_and_diff(totalmri)
print((np.array(mmri)-np.array(mecho))/(np.array(mecho))*100)
plt.figure(figsize=(8,6))
plt.bar(res-0.2, mecho, width = 0.4, label = 'ECHO', color='blue')
plt.bar(res+0.2, mmri, width = 0.4, label = 'MRI', color='darkorange')
plt.errorbar(res-0.2, mecho, yerr=stdecho, fmt="o", color="r")
plt.errorbar(res+0.2, mmri, yerr=stdmri, fmt="o", color="r")
plt.legend(loc='lower center', borderaxespad= -4, frameon=False, ncol=2, numpoints=4,\
        handletextpad=0.05,prop={'family' : 'Times New Roman', 'size' : '24'})
plt.ylabel(r'Contractility (kPa)', font1)
y_ticks = np.linspace(1, max(max(mecho+abs(np.array(stdecho))),max(mmri+abs(np.array(stdmri)))), 3)
x_ticks = np.linspace(1, 16, 16)
plt.yticks(y_ticks)
plt.xticks(x_ticks)
plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.tight_layout()  
plt.savefig('Max Contractility')
plt.close('all')

mecho=np.mean(ttecho, axis = 0)
mmri=np.mean(ttmri, axis = 0)
stdecho=np.std(ttecho, axis = 0)
stdmri=np.std(ttmri, axis = 0)
plt.figure(figsize=(8,6))
plt.bar(res-0.2, mecho, width = 0.4, label = 'ECHO', color='blue')
plt.bar(res+0.2, mmri, width = 0.4, label = 'MRI', color='darkorange')
plt.errorbar(res-0.2, mecho, yerr=stdecho, fmt="o", color="r")
plt.errorbar(res+0.2, mmri, yerr=stdmri, fmt="o", color="r")
plt.legend(loc='lower center', borderaxespad= -4, frameon=False, ncol=2, numpoints=4,\
        handletextpad=0.05,prop={'family' : 'Times New Roman', 'size' : '24'})
plt.ylabel(r'Time (s)', font1)
y_ticks = np.linspace(0, max(max(mecho+stdecho),max(mmri+stdmri)), 3)
x_ticks = np.linspace(1, 16, 16)
plt.yticks(y_ticks)
plt.xticks(x_ticks)
plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.tight_layout()  
plt.savefig('Time to peak')
plt.close('all')
