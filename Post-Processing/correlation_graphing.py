import numpy as np
import matplotlib.pyplot as plt

aT_echo = [275.0495843, 189.0816201, 215.4755936, 262.9338237, 189.1226136, 281.9576813, 218.0716974]   #active tension in ECHO
aT_mri = [229.7009363, 160.7534176, 176.8056527, 139.7658522, 158.1578118, 194.4538006, 173.672695]
EDV_echo = [140.06, 52.95, 84.23, 72.45, 75.02, 59.04, 68.46]
EDV_mri = [97.407492, 75.187678, 63.06153, 61.969764, 61.784404, 64.88807, 75.352568]
ESV_echo = [53.35732143, 19.95, 32.17670384, 21.98026022, 35.57664808, 22.64004647, 29.06628324]
ESV_mri = [42.92906, 37.15827647, 26.69678, 21.134951, 33.9853134, 40.9584486, 42.613271]
Perror_echo = [0.175178682, 0.120147454, 0.216753395, 0.144770112, 0.17070662, 0.260716651, 0.145156943]
Perror_mri = [0.209379757, 0.227028281, 0.176200786, 0.20833536, 0.215350689, 0.199961955, 0.203534967]
Verror_echo = [0.06619639, 0.061814603, 0.066989418, 0.082223709, 0.058540396, 0.060403981, 0.055864997]
Verror_mri = [0.058483123, 0.053436365, 0.082623658, 0.089350635, 0.040527157, 0.039236701, 0.050048113]
case = ['1','2','3','4','5','6','7']    #case = ['396669','396670','399018','399019','407034','407036', '407037']
part = ['ECHO','MRI']; filename1 = []
trendlinex = [] ; trendliney = []
font1 = {'weight' : 'normal',
'fontname' : 'Times New Roman',
'size' : 24,}
####################correlation of active tension between ECHO and MRI
max_x = max(aT_echo)
max_y = max(aT_mri)
maxrange = max_x
if max_x < max_y:
    maxrange = max_y
filename1.append('Active Tension')
plt.figure(figsize=(8,6))
z = np.polyfit(aT_echo, aT_mri, 1)
p = np.poly1d(z)
func= ("MRI = %.3fECHO+(%.3f)"%(z[0],z[1]))
x = np.linspace(130, max_x, 10)          # process identical data point for pressure, volume, and strain
plt.plot(x, p(x),'-', label=func)
plt.plot(aT_echo, aT_mri,'o', markersize=5, markeredgewidth=5 ,color='blue')
plt.legend(loc='lower center', borderaxespad= -4, frameon=False, numpoints=1,\
           handletextpad=0.05,prop={'family' : 'Times New Roman', 'size' : '24'})
plt.ylabel(r'MRI (mmHg)', font1)
plt.xlabel(r'ECHO (mmHg)', font1)
y_ticks = np.linspace(130, maxrange, 3) # three values for y-axis for the largest value
x_ticks = np.linspace(130, maxrange, 3) # three values for x-axis for the largest value
plt.yticks(y_ticks)
plt.xticks(x_ticks)
plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.title("Active Tension in ECHO and MRI", fontsize=24, fontname="Times New Roman")
plt.tight_layout()
plt.savefig(filename1[0])
####################correlation between EDV and active tension in ECHO
max_x1 = max(EDV_echo)
max_x2 = max(ESV_echo)
max_y1 = max(EDV_mri)
max_y2 = max(ESV_mri)
max_x = max_x1
max_y = max_y1
if max_x2 > max_x1:
    max_x = max_x2
if max_y2 > max_y1:
    max_y = max_y2
filename1.append('EDV&ESV')
plt.figure(figsize=(8,6))
z1 = np.polyfit(EDV_echo, EDV_mri, 1)
z2 = np.polyfit(ESV_echo, ESV_mri, 1)
p1 = np.poly1d(z1)
p2 = np.poly1d(z2)
func1= ("MRI = %.3fECHO+(%.3f)"%(z1[0],z1[1]))
func2= ("MRI = %.3fECHO+(%.3f)"%(z2[0],z2[1]))
x1 = np.linspace(0, max_x, 10)     
x2 = np.linspace(0, max_x, 10)     
plt.plot(x1, p1(x1),'-', label=func1, color='blue')
plt.plot(x2, p2(x2),'-', label=func2, color='darkorange')
plt.plot(EDV_echo, EDV_mri,'o', markersize=5, markeredgewidth=5 , label = 'EDV', color='blue')
plt.plot(ESV_echo, ESV_mri,'o', markersize=5, markeredgewidth=5 , label = 'ESV', color='darkorange')
plt.legend(loc=(0, 0.57),  frameon=False, numpoints=1,\
           handletextpad=0.05,prop={'family' : 'Times New Roman', 'size' : '24'})
plt.ylabel(r'MRI Volume (mL)', font1)
plt.xlabel(r'ECHO Volume (mL)', font1)
y_ticks = np.linspace(0, max_y, 3)  
x_ticks = np.linspace(0, max_x, 3)
plt.yticks(y_ticks)
plt.xticks(x_ticks)
plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.title("EDV & ESV in ECHO and MRI", fontsize=24, fontname="Times New Roman")
plt.tight_layout()
plt.savefig(filename1[1])
############################Comparison of error in pressure
m=0
while (m!=len(Perror_echo)):
    Perror_echo[m]=Perror_echo[m]*100
    Perror_mri[m]=Perror_mri[m]*100
    Verror_echo[m]=Verror_echo[m]*100
    Verror_mri[m]=Verror_mri[m]*100
    m+=1
max_x1 = max(Perror_echo)
max_y1 = max(Perror_mri)
max_x2 = max(Verror_echo)
max_y2 = max(Verror_mri)
maxrange = max(max_x1,max_x2,max_y1,max_y2)
filename1.append('PVLoop Error')
plt.figure(figsize=(8,6))
plt.plot(case, Perror_echo,'o', markersize=5, markeredgewidth=5, label = 'ECHO', color='blue')
plt.plot(case, Perror_mri,'+', markersize=5, markeredgewidth=5, label = 'MRI', color='blue')
plt.plot(case, Verror_echo,'o', markersize=5, markeredgewidth=5, label = 'ECHO', color='darkorange')
plt.plot(case, Verror_mri,'+', markersize=5, markeredgewidth=5, label = 'MRI', color='darkorange')
plt.legend(loc='lower center', borderaxespad= -4, frameon=False, ncol=2, numpoints=4,\
           handletextpad=0.05,prop={'family' : 'Times New Roman', 'size' : '24'})
plt.ylabel(r'Pressure & Volume (%)', font1)
y_ticks = np.linspace(0, maxrange, 3)
plt.yticks(y_ticks)
plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.tight_layout()
plt.savefig(filename1[2])
##########################Separation from 'PVLoop Error'
##############
x_axis = np.arange(len(case))
max_x = max(Perror_echo)
max_y = max(Perror_mri)
maxrange = max(max_x,max_y)
filename1.append('Pressure Error')
plt.figure(figsize=(8,6))
plt.bar(x_axis-0.2, Perror_echo, width = 0.4, label = 'ECHO', color='blue')
plt.bar(x_axis+0.2, Perror_mri, width = 0.4, label = 'MRI', color='darkorange')
plt.legend(loc='lower center', borderaxespad= -4, frameon=False, ncol=2, numpoints=4,\
           handletextpad=0.05,prop={'family' : 'Times New Roman', 'size' : '24'})
plt.ylabel(r'Pressure (%)', font1)
y_ticks = np.linspace(0, maxrange, 3)
plt.yticks(y_ticks)
plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.title("Percentage Error of Pressure", fontsize=24, fontname="Times New Roman")
plt.tight_layout()
plt.savefig(filename1[3])
##############
x_axis = np.arange(len(case))
max_x = max(Verror_echo)
max_y = max(Verror_mri)
maxrange = max(max_x,max_y)
filename1.append('Volume Error')
plt.figure(figsize=(8,6))
plt.bar(x_axis-0.2, Verror_echo, width = 0.4, label = 'ECHO', color='blue')
plt.bar(x_axis+0.2, Verror_mri, width = 0.4, label = 'MRI', color='darkorange')
plt.legend(loc='lower center', borderaxespad= -4, frameon=False, ncol=2, numpoints=4,\
           handletextpad=0.05,prop={'family' : 'Times New Roman', 'size' : '24'})
plt.ylabel(r'Volume (%)', font1)
y_ticks = np.linspace(0, maxrange, 3)
plt.yticks(y_ticks)
plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
plt.title("Percentage Error of Volume", fontsize=24, fontname="Times New Roman")
plt.tight_layout()
plt.savefig(filename1[4])