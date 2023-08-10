import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

font1 = {'weight' : 'normal',
'fontname' : 'Times New Roman',
'size' : 24,}

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

    file =pd.read_csv(str(swine)+'wallthickness.csv')
    wecho = np.array(file.values[:,0])
    wmri = np.array(file.values[:,1])
    echowall.append(wecho)
    mriwall.append(wmri)


mecho = np.mean(echowall, axis = 0)
mmri = np.mean(mriwall, axis = 0)
stdecho = np.std(echowall, axis = 0)
stdmri = np.std(mriwall, axis = 0)
res = np.arange(1,17)
plt.figure(figsize=(8,6))
plt.bar(res-0.2, mecho, width=0.4, label='ECHO',color='blue')
plt.bar(res+0.2, mmri, width=0.4, label='MRI', color='darkorange')
plt.errorbar(res-0.2, mecho, yerr=stdecho, fmt="o", color="r")
plt.errorbar(res+0.2, mmri, yerr=stdmri, fmt="o", color="r")
plt.legend(loc='lower center', borderaxespad= -4, frameon=False, ncol=2, numpoints=4,\
        handletextpad=0.05,prop={'family' : "Times New Roman", 'size' : '24'})
plt.ylabel(r'Thickness (cm)', font1)
y_ticks = np.linspace(0, max(max(mecho+stdecho),max(mmri+stdmri)), 3)
x_ticks = np.linspace(1, 17, 17)
plt.yticks(y_ticks)
plt.xticks(x_ticks)
plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.tight_layout()  
plt.savefig('total wall')
plt.close('all')
print(mecho,mmri)
print((np.mean(mecho[:5]),np.mean(mmri[:5])))

