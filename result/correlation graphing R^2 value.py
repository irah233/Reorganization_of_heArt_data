import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm
import math
import random
import scipy.stats as stats

aT_echo = []   #active tension in ECHO
aT_mri = [] 
EDV_mriwo =[97.407492, 75.187678, 63.06153, 61.969764, 61.784404, 64.88807, 75.352568]
EDV_echo = [140.06, 52.95, 84.23, 72.45, 75.02, 59.04, 68.46]
EDV_mriw = [156.977492, 94.587678, 85.61153, 84.519764, 73.964404, 86.72807, 97.852568]#[97.407492, 75.187678, 63.06153, 61.969764, 61.784404, 64.88807, 75.352568]
ESV_echo = [53.35732143, 19.95, 32.17670384, 21.98026022, 35.57664808, 22.76061958, 28.22]
ESV_mriwo = [42.92906, 35.14480307, 26.69678, 21.134951, 35.636926, 37.343256, 42.613271]
ESV_mriw = [65.44906, 24.18480307, 16.96678, 11.404951, 42.766926, 44.493256, 53.003271]
EDV_echo_t = [126.224, 40.5, 76.34, 62.589, 59.028, 49.811, 57.404]
ESV_echo_t = [39.52132143, 11.4, 24.28670384, 13.18950186, 19.58464808, 13.53161958, 21.394]
EDV_mri_t = [76.342492, 53.673678, 57.01253, 57.009764, 55.912404, 58.04707, 55.612568]
ESV_mri_t = [21.86406, 15.64427647, 20.64778, 16.174951, 29.764926, 30.502256, 22.873271]

EDV_echo=np.array(EDV_echo)
ESV_echo=np.array(ESV_echo)
EDV_mriw=np.array(EDV_mriw)
ESV_mriw=np.array(ESV_mriw)
EDV_mriwo=np.array(EDV_mriwo)
ESV_mriwo=np.array(ESV_mriwo)
EDV_echo_t=np.array(EDV_echo_t)
ESV_echo_t=np.array(ESV_echo_t)
EDV_mri_t=np.array(EDV_mri_t)
ESV_mri_t=np.array(ESV_mri_t)

case = ['1','2','3','4','5','6','7']    #case = ['396669','396670','399018','399019','407034','407036', '407037']
cases = ['396669','396670','399018','399019','407034','407036', '407037']
part = ['ECHO','MRI']; filename1 = []; filenames=[]
trendlinex = [] ; trendliney = []
font1 = {'weight' : 'normal',
'fontname' : 'Times New Roman',
'size' : 24,}
font2 = {'weight' : 'normal',
'fontname' : 'Times New Roman',
'size' : 30,}
ECHO_t=[];MRI_t=[]
Pe1=[];Ve1=[];Pe2=[];Ve2=[]
time_e=[];time_m=[]
all_color=['indianred','maroon','red','sienna','chocolate','goldenrod',
       'yellowgreen','forestgreen','lime','lightseagreen','cyan','deepskyblue',
       'navy','blueviolet','indigo','magenta']


def ttest(data_group1,data_group2):

     # Perform the two sample t-test with equal variances
    retuslt = stats.ttest_ind(a=data_group1, b=data_group2, equal_var=True)
    print(retuslt)

def pact_info(swine,time):
    if swine == '396669':
        total_t = 0.66925
        lo = 1.08
    elif swine == '396670':
        total_t = 0.4845
        lo = 1.08
    elif swine == '399018':
        total_t = 0.5383
        lo = 1.18
    elif swine == '399019':
        total_t = 0.5383
        lo = 0.88
    elif swine == '407034':
        total_t = 0.5383
        lo = 1.08
    elif swine == '407036':
        total_t = 0.48445
        lo = 1.08
    elif swine == '407037':
        total_t = 0.5921
        lo = 1.18
    return total_t, lo

def ActiveLength(lamda, lr, lo):
    ls = lamda*lr;
    if ls<=lo:
        lso = 0.002;
    else:
        lso = ls-lo
    return lso

def sqrt(num):
    return math.sqrt(num)

def exp(num):
    return math.exp(num)

def cos(num):
    return math.cos(num)

def find_pact(pig_,ECHO_t,MRI_t,tt,lo,tl):
    pi = math.pi
    lr = 1.85;
    CaoMax = 4.35;
    Cao = 4.35;
    B = 4.75;
    to = 275;
    ttrans = 300;
    tau = 25;
    BCL = 750;
    lambda_ = 1.5
    techo = ECHO_t
    tmri = MRI_t
    ls_Pm=[];
    ls_Pe=[];
    lso = ActiveLength(lambda_, lr, lo)
    deno = sqrt(exp((B*lso)-1))
    ECa50 = CaoMax/deno
    for i in range(0,len(tl)):
        at1=[]
        at2=[]
        for j in range(0,BCL):
            ta = j
            if ta<ttrans:
                 Ct = 0.5*(1-cos(pi*ta/to))
            else:
                 Ct = 0.5*(1-cos(pi*ttrans/to))*exp(-((ta-ttrans)/tau))
            CaTerm = (Cao**2)/(Cao**2 + ECa50**2)
            Pacte = techo[i]*CaTerm*Ct
            Pactm = tmri[i]*CaTerm*Ct
            at1.append(Pacte)
            at2.append(Pactm)
        ls_Pe.append(max(at1))
        ls_Pm.append(max(at2))
    ls_Pe = np.array(ls_Pe)
    ls_Pm = np.array(ls_Pm)
    filename1 = []
    filename1.append(str(pig_)+'Active Tension')
    plt.figure(figsize=(7,5))
    plt.plot(tl, ls_Pe,'-', markersize=5, markeredgewidth=5,label='ECHO' ,color='blue')
    plt.plot(tl, ls_Pm,'--', markersize=5, markeredgewidth=5,label ='MRI' ,color='darkorange')
    plt.legend(loc='upper right',  frameon=False, ncol=1, fontsize=24)# borderaxespad= -3, bbox_to_anchor=(0.3, 0),
    plt.ylabel(r'Active Tension (mmHg)', font1)
    plt.xlabel(r'Time (ms)', font1)
    y_ticks = np.linspace(0, max(max(ls_Pm),max(ls_Pe)), 3)
    x_ticks = np.linspace(0, max(tl), 3)
    plt.yticks(y_ticks)
    plt.xticks(x_ticks)
    plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
    plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
    plt.tight_layout()
    plt.savefig(filename1[0])
    plt.close('all')
    return ls_Pe, ls_Pm

def func(x, a):
    return a * x
####################Altman-Bland Graph
def al_bl(mean,diff,md1,md,sd1,sd,type_,unit,name_p):
    filenames = []
    filenames.append(type_+'_difff mean')
    plt.figure(figsize=(9,7))
    plt.scatter(mean,diff/mean*100,color="b")
    plt.axhline(md1,linestyle='-')
    plt.axhline(md1+1.96*sd1,linestyle='--')
    plt.axhline(md1-1.96*sd1,linestyle='--')
    plt.ylabel(name_p+r' (ECHO-MRI)/Mean (%)', font2)
    plt.xlabel(name_p+r' Mean '+unit, font2)
    plt.yticks(fontsize=30, rotation=0, fontname="Times New Roman")
    plt.xticks(fontsize=30, rotation=0, fontname="Times New Roman")
    y_ticks = np.linspace(md1-1.96*sd1, md1+1.96*sd1, 3) 
    plt.yticks(y_ticks)
    x_ticks = np.linspace(min(mean), max(mean), 3) 
    plt.xticks(x_ticks)
    plt.tight_layout()
    plt.savefig(filenames[0])
    plt.close('all')
    filenames.append(type_+'_diff')
    plt.figure(figsize=(9,7))
    plt.scatter(mean,diff,color="b")
    plt.axhline(md,linestyle='-')
    plt.axhline(md+1.96*sd,linestyle='--')
    plt.axhline(md-1.96*sd,linestyle='--')
    plt.ylabel(name_p+r' ECHO - MRI '+unit, font2)
    plt.xlabel(name_p+r' Mean '+unit, font2)
    plt.yticks(fontsize=30, rotation=0, fontname="Times New Roman")
    plt.xticks(fontsize=30, rotation=0, fontname="Times New Roman")
    y_ticks = np.linspace(md-1.96*sd, md+1.96*sd, 3) 
    plt.yticks(y_ticks)
    x_ticks = np.linspace(min(mean), max(mean), 3) 
    plt.xticks(x_ticks)
    plt.tight_layout()
    plt.savefig(filenames[1])
    plt.close('all')
    
def albl2(echo,mri,mean1,mean2,diff,type_,unit,e_std1,m_std1):
    filenames = []
    filenames.append(type_+'_diff ECHO and MRI')
    plt.scatter(mean1,echo,color="b")
    plt.scatter(mean2,mri,'<',color="b")
    plt.axhline(mean1,linestyle='-')
    plt.axhline(mean1+1.96*e_std1,linestyle='--')
    plt.axhline(mean1-1.96*e_std1,linestyle='--')
    plt.axhline(mean2,linestyle='-')
    plt.axhline(mean2+1.96*m_std1,linestyle='--')
    plt.axhline(mean2-1.96*m_std1,linestyle='--')
    plt.ylabel(type_+r' ECHO - MRI '+unit, font1)
    plt.xlabel(type_+r' Mean '+unit, font1)
    plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
    plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
    y_ticks = np.linspace(max(mean1+1.96*e_std1,mean2+1.96*m_std1), min(mean1-1.96*e_std1,mean2-1.96*m_std1), 3) 
    plt.yticks(y_ticks)
    plt.tight_layout()
    plt.savefig(filenames[0])
    plt.close('all')
####################correlation for only one valariable for both 3D ECHO and MRI
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
    plt.plot(x, func(x,res.params),'-',  lw = 10, label="MRI = %.3fECHO" % res.params[0], color=color1)
    plt.plot(echo, mri,'o', markersize=5, markeredgewidth=5 ,color=color1)
    #plt.text(minrange+(maxrange-minrange)/100, maxrange-minrange*3/6, "R = %.3f" % float(math.sqrt(res.rsquared)),fontdict = font1)
    plt.text(minrange+(maxrange-minrange)/100, maxrange-minrange*1/6, r"$R^2$ = %.3f" % float((res.rsquared)),fontdict = font1)
    plt.legend(loc='upper left', frameon=False, numpoints=1,\
               handletextpad=0.05,prop={'family' : 'Times New Roman', 'size' : '24'})
    plt.ylabel(name_p+r' MRI ('+unit+')', font1)
    plt.xlabel(name_p+r' ECHO ('+unit+')', font1)
    y_ticks = np.linspace(minrange, maxrange, 3) # three values for y-axis for the largest value
    x_ticks = np.linspace(minrange, maxrange, 3) # three values for x-axis for the largest value
    plt.yticks(y_ticks)
    plt.xticks(x_ticks)
    plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
    plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
    plt.tight_layout()
    plt.savefig(filename1[0])
    plt.close('all')

####################correlation for both EDV and ESV
def twoline(EDV_echo,ESV_echo,EDV_mriw,ESV_mriw,type_):
    filename1 = []
    max_x1 = max(EDV_echo)
    max_x2 = max(ESV_echo)
    max_y1 = max(EDV_mriw)
    max_y2 = max(ESV_mriw)
    max_x = max_x1
    max_y = max_y1
    if max_x2 > max_x1:
        max_x = max_x2
    if max_y2 > max_y1:
        max_y = max_y2
    max_xy = max(max_x,max_y)
    filename1.append(type_)
    plt.figure(figsize=(8,6))
    res_d = sm.OLS(EDV_mriw, EDV_echo).fit()
    
    res_s = sm.OLS(ESV_mriw, ESV_echo).fit()
    
    x1 = np.linspace(max_x2, max_y+50, 10)     
    x2 = np.linspace(0, max_x2, 10)     
    plt.plot(x1, func(x1,res_d.params),'-', lw=8,label="MRI = %.3fECHO" % res_d.params[0], color='blue')
    plt.plot(x2, func(x2,res_s.params),'-', lw=8,label="MRI = %.3fECHO" % res_s.params[0], color='darkorange')
    plt.plot(EDV_echo, EDV_mriw,'o', markersize=8, markeredgewidth=5 , label = 'EDV', color='blue')
    plt.plot(ESV_echo, ESV_mriw,'o', markersize=8, markeredgewidth=5 , label = 'ESV', color='darkorange')
    #plt.text(minrange+(maxrange-minrange)/100, maxrange-minrange*3/6, "R^2 = %.3f" % float((res.rsquared)),fontdict = font1)
    plt.text(80, 60, r"$R^2$ = %.3f" % float((res_d.rsquared)),fontdict = font1, color = 'blue')
    #plt.text(minrange+(maxrange-minrange)/100, maxrange-minrange*3/6, "R^2 = %.3f" % float((res.rsquared)),fontdict = font1)
    plt.text(80, 48, r"$R^2$ = %.3f" % float((res_s.rsquared)),fontdict = font1, color = 'darkorange')
    plt.legend(loc=(0, 0.56),  frameon=False, numpoints=1,\
               handletextpad=0.05,prop={'family' : 'Times New Roman', 'size' : '24'})
    plt.ylabel(r'MRI Volume (mL)', font1)
    plt.xlabel(r'ECHO Volume (mL)', font1)
    y_ticks = np.linspace(min(ESV_mriw), max_xy, 3)  
    x_ticks = np.linspace(0, max_x, 3)
    plt.yticks(y_ticks)
    plt.xticks(x_ticks)
    plt.yticks(fontsize=24, rotation=0, fontname="Times New Roman")
    plt.xticks(fontsize=24, rotation=0, fontname="Times New Roman")
    plt.tight_layout()
    plt.savefig(filename1[0])
    plt.close('all')
##########################Separation from 'PVLoop Error'
def error_bar(echo,mri,type_):
    filename1=[]
    x_axis = np.arange(len(case))
    max_x = max(echo)
    max_y = max(mri)
    maxrange = max(max_x,max_y)
    filename1.append(type_)
    plt.figure(figsize=(8,6))
    plt.bar(x_axis-0.2, echo, label = 'ECHO', color='blue')
    plt.bar(x_axis+0.2, mri, label = 'MRI', color='darkorange')
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
    
def mean_diff(item_e, item_m):
    echot=np.array(item_e)
    mrit=np.array(item_m)
    mean1=np.mean(echot)
    mean2=np.mean(mrit)
    e_std1=np.std(echot)
    m_std1=np.std(mrit)
    mean=np.mean([item_e,item_m],axis=0)
    diff=(echot-mrit)
    md=np.mean(diff)
    sd=np.std(diff,axis=0)
    md1=np.mean(diff/mean)*100
    sd1=np.std(diff/mean,axis=0)*100
    return mean1,mean2,mean,diff,md1,md,sd1,sd,e_std1,m_std1
    
for i in cases:
    pig_ = i
    csre=pd.read_csv(str(i)+'echo/'+'Tmax.csv')
    csrm=pd.read_csv(str(i)+'mri/'+'Tmax.csv')
    time = np.array(csre.values[:,0])
    ECHO=np.array(csre.values[:,1])
    MRI=np.array(csrm.values[:,1])
    tt, lo = pact_info(pig_,time)
    tl = time*tt/len(time)*1000
    ECHO_t.append(max(ECHO))
    MRI_t.append(max(MRI))
    ls_Pe, ls_Pm = find_pact(pig_,ECHO*0.0075,MRI*0.0075,tt,lo,tl)
    e1 = max(ECHO); m1 = max(MRI)
    
    inde = (ECHO.tolist()).index(e1);
    indm = (MRI.tolist()).index(m1);
    time_e.append(time[inde]*tt/len(time)*1000)
    time_m.append(time[indm]*tt/len(time)*1000)
    
    aT_echo.append(max(ls_Pe))
    aT_mri.append(max(ls_Pm))
    kind=['echo','mri']
    for na in kind:
        if na.lower().strip() == 'echo':
            path = (pig_)+'echo'
        elif na.lower().strip() == 'mri':
            path = (pig_)+'mri'
        PVdata = pd.read_csv(path+'/PVLoop.csv')       
        PVbdata = pd.read_csv(path+'/'+(pig_)+'.csv')
        t = np.array(PVbdata.values[:,0][:30]) 
        V = np.array(PVbdata.values[:,2][:30])
        P = np.array(PVbdata.values[:,1][:30])                   
        LVP = np.array(PVdata.values[:,1][18:])           
        LVV = np.array(PVdata.values[:,0][18:])
        V_t = t
        Vsum=Psum=0
        for i in range(len(LVV)):
            Vsum += abs(V[i] - LVV[i])/max(V)*100
            Psum += abs(P[i] - LVP[i])/max(P)*100
        Verror = Vsum/len(LVV)
        Perror = Psum/len(LVP)
        if na.lower().strip() == 'echo':
            Pe1.append(Perror)
            Ve1.append(Verror)
        if na.lower().strip() == 'mri':
            Pe2.append(Perror)
            Ve2.append(Verror)
#print(np.array(ECHO_t)/1000)
#print(np.array(MRI_t)/1000)
print(ECHO_t.pop(3))
print(MRI_t.pop(3))
EF_e = (EDV_echo_t-ESV_echo_t)/EDV_echo_t*100
EF_m = (EDV_mri_t-ESV_mri_t)/EDV_mri_t*100
EF_m1 = (EDV_mriwo-ESV_mriwo)/EDV_mriwo*100
time_e = np.array(time_e)
time_m = np.array(time_m)
####################contractility
mean1,mean2,mean,diff,md1,md,sd1,sd,e_std1,m_std1 = mean_diff(np.array(ECHO_t)/1000,np.array(MRI_t)/1000)
unit = '(kPa)'
al_bl(mean,diff,md1,md,sd1,sd,'Tmax',unit,'Tmax')
'''albl2(ECHO_t,MRI_t,mean1,mean2,diff,'Tmax',unit,e_std1,m_std1)'''
####################active tension
mean1,mean2,mean,diff,md1,md,sd1,sd,e_std1,m_std1 = mean_diff(np.array(aT_echo),np.array(aT_mri))
unit = '(kPa)'
al_bl(mean,diff,md1,md,sd1,sd,'Pact',unit,'Pact')
####################EDVw
mean1,mean2,mean,diff,md1,md,sd1,sd,e_std1,m_std1 = mean_diff(EDV_echo,EDV_mriw)
unit = '(mL)'
al_bl(mean,diff,md1,md,sd1,sd,'EDV with CSO',unit,'EDV')
####################EDVwo
mean1,mean2,mean,diff,md1,md,sd1,sd,e_std1,m_std1 = mean_diff(EDV_echo,EDV_mriwo)
unit = '(mL)'
al_bl(mean,diff,md1,md,sd1,sd,'EDV without CSO',unit,'EDV')
####################ESVw
mean1,mean2,mean,diff,md1,md,sd1,sd,e_std1,m_std1 = mean_diff(ESV_echo,ESV_mriw)
unit = '(mL)'
al_bl(mean,diff,md1,md,sd1,sd,'ESV with CSO',unit,'ESV')
####################ESVwo
mean1,mean2,mean,diff,md1,md,sd1,sd,e_std1,m_std1 = mean_diff(ESV_echo,ESV_mriwo)
unit = '(mL)'
al_bl(mean,diff,md1,md,sd1,sd,'ESV without CSO',unit,'ESV')
####################SV
mean1,mean2,mean,diff,md1,md,sd1,sd,e_std1,m_std1 = mean_diff(EDV_echo_t-ESV_echo_t,EDV_mri_t-ESV_mri_t)
unit = '(mL)'
al_bl(mean,diff,md1,md,sd1,sd,'SV truncated',unit,'SV')
####################EF
mean1,mean2,mean,diff,md1,md,sd1,sd,e_std1,m_std1 = mean_diff(EF_e,EF_m)
unit = '(%)'
al_bl(mean,diff,md1,md,sd1,sd,'EF truncated',unit,'EF')
####################Time to Peak Tmax
mean1,mean2,mean,diff,md1,md,sd1,sd,e_std1,m_std1 = mean_diff(time_e,time_m)
unit = '(ms)'
al_bl(mean,diff,md1,md,sd1,sd,'Tmax time',unit,'Tmax time')
####################EDV_t
mean1,mean2,mean,diff,md1,md,sd1,sd,e_std1,m_std1 = mean_diff(EDV_echo_t[1:],EDV_mri_t[1:])
unit = '(mL)'
al_bl(mean,diff,md1,md,sd1,sd,'EDV truncated',unit,'EDV')
####################ESV_T
mean1,mean2,mean,diff,md1,md,sd1,sd,e_std1,m_std1 = mean_diff(ESV_echo_t,ESV_mri_t)
unit = '(mL)'
al_bl(mean,diff,md1,md,sd1,sd,'ESV truncated',unit,'ESV')
########################################one correlate line
oneline(np.array(ECHO_t)/1000,np.array(MRI_t)/1000,'Contractility','kPa','Tmax')
oneline(np.array(aT_echo)/1000,np.array(aT_mri)/1000,'Active Tension','kPa','Pact')
oneline(EDV_echo-ESV_echo,EDV_mri_t-ESV_mri_t,'Stroke Volume truncated','mL','SV')
oneline(EF_e,EF_m,'Ejection Fraction truncated','%','EF')
oneline(EF_e,EF_m1,'Ejection Fraction wo CSO','%','EF')
oneline(time_e,time_m,'global time to peak','ms','Time to peak Tmax')
########################################two correlate line
twoline(EDV_echo,ESV_echo,EDV_mriw,ESV_mriw,'EDV&ESV with CSO')
twoline(EDV_echo,ESV_echo,EDV_mriwo,ESV_mriwo,'EDV&ESV without CSO')
twoline(EDV_echo_t[1:],ESV_echo_t,EDV_mri_t[1:],ESV_mri_t,'EDV&ESV truncated')
########################################error bar
error_bar(np.array(Pe1),np.array(Pe2),'Pressure Error')
ECHO_t = np.array(ECHO_t)/1000,
MRI_t = np.array(MRI_t)/1000

ttest(EDV_echo,EDV_mri_t)
ttest(ESV_echo,ESV_mri_t)
ttest(EDV_echo-ESV_echo,EDV_mri_t-ESV_mri_t)
ttest(EF_e,EF_m)

