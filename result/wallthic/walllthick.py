import pandas as pd
import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

font1 = {'weight' : 'normal',
'fontname' : 'Times New Roman',
'size' : 24,}
case = [396669,396670,399018,399019,407034,407036,407037]
echo=[]
mri=[]
totale=[]
totalm=[]
for swine in case:
    csr=pd.read_csv(str(swine)+'wallthickness.csv')
    
    ECHO=np.array(csr.values[:,0])
    MRI=np.array(csr.values[:,1])
    mean=np.mean([ECHO,MRI],axis=0)
    diff=(ECHO-MRI)
    md=np.mean(diff)
    sd=np.std(diff,axis=0)
    filename1=[]
    md1=np.mean(diff/mean)*100
    sd1=np.std(diff/mean,axis=0)*100
    totale.append(ECHO)
    totalm.append(MRI)
    echo.append(np.mean(ECHO))
    mri.append(np.mean(MRI))


    filename1.append(str(swine)+'wallthickness mean')
    plt.figure(figsize=(8,6))
    plt.scatter(mean,diff/mean*100,color="b")
    plt.axhline(md1,linestyle='-')
    plt.axhline(md1+1.96*sd1,linestyle='--')
    plt.axhline(md1-1.96*sd1,linestyle='--')
    plt.ylabel(r'Difference/Mean (%)', font1)
    plt.xlabel(r'Mean (cm)', font1)
    plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
    plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
    y_ticks = np.linspace(md1-1.96*sd1, md1+1.96*sd1, 3) 
    plt.yticks(y_ticks)
    plt.tight_layout()
    plt.savefig(filename1[0])
    plt.close('all')
    
    
    filename1.append(str(swine)+'wallthickness')
    plt.figure(figsize=(8,6))
    plt.scatter(mean,diff,color="darkorange")
    plt.axhline(md,linestyle='-')
    plt.axhline(md+1.96*sd,linestyle='--')
    plt.axhline(md-1.96*sd,linestyle='--')
    plt.ylabel(r'Difference (cm)', font1)
    plt.xlabel(r'Mean (cm)', font1)
    plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
    plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
    y_ticks = np.linspace(md-1.96*sd, md+1.96*sd, 3) 
    plt.yticks(y_ticks)
    plt.tight_layout()
    plt.savefig(filename1[1])
    plt.close('all')
    
    filename1.append(str(swine)+'wallthicknessEM')
    plt.figure(figsize=(8,6))
    max_x = max(ECHO)
    max_y = max(MRI)
    maxrange = max_x
    if max_x < max_y:
        maxrange = max_y
    min_x = min(ECHO)
    min_y = min(MRI)
    minrange = min_x
    if min_x > min_y:
        minrange = min_y
    z = np.polyfit(ECHO, MRI, 1)
    p = np.poly1d(z)
    func= ("MRI = %.3fECHO+(%.3f)"%(z[0],z[1]))
    x = np.linspace(minrange, max_x, 10)          # process identical data point for pressure, volume, and strain
    plt.plot(x, p(x),'-', label=func)
    plt.plot(ECHO,MRI,'o', markersize=5, markeredgewidth=5,color='blue')
    plt.legend(loc='upper center', frameon=False, numpoints=1,\
               handletextpad=0.05,prop={'family' : 'Times New Roman', 'size' : '24'})
    plt.ylabel(r'MRI (cm)', font1)
    plt.xlabel(r'ECHO (cm)', font1)
    y_ticks = np.linspace(minrange, maxrange, 3) # three values for y-axis for the largest value
    x_ticks = np.linspace(minrange, maxrange, 3) # three values for x-axis for the largest value
    plt.yticks(y_ticks)
    plt.xticks(x_ticks)
    plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
    plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
    plt.tight_layout()
    plt.savefig(filename1[2])
    plt.close('all')  
    
print(echo[6])
print(mri[6])
print(np.std(np.array(echo)-np.array(mri)))
print(np.std(np.array(echo)))
print(np.std(np.array(mri)))

'''
print(np.std(totale))
print(np.std(totalm))
print(avgecho,avgmri)
'''
