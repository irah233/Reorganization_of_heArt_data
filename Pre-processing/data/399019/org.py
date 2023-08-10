import csv
import pandas as pd
import numpy as np
from matplotlib import pylab as plt
import os as os
import scipy.io as sio
import scipy.io
from scipy.interpolate import interp1d

##  test process
pignumber=399019           #running case
data_number = 31         #process data numbers for simulation
P_filename = 'LVpressure'       #name for pressure file
txt_filename = '399019'
ms=3
try:
    print((txt_filename+'.csv') == True)
    open(txt_filename+'.csv')
except: 
    extra='_tomtec'
    ## Delete corresponding terms if no such condition in the case
    case_name_1=['Base','RAP','RVP','S75','S75+CRT','S100','S100+CRT']
    case_name_2=['Base','RAP','RVP','S75','S75+CRT']
    case_name_3=['Base','RAP','RVP','S100','S100+CRT']
    case_name_4=['RAP_100','RVP_100']
    case_name_5=['Base','RAP100','RAP140','RVP100','RVP140']
    case_name_6=['Base','RAP100','RAP140','RVP100']
    case_name_7=['Base','LCX+RAP+RVP_100',\
                 'LCX+RAP+RVP+LVP(Endo_Mid_Anterior_Wall)_100',\
                 'LCX+RAP+RVP+LVP(Endo_Mid_Lateral_Wall)_100',\
                 'LCX+RAP+RVP+LVP(Epi_Mid_Anterior_Wall)_100',\
                 'RAP_100','RAP+RVP_100',\
                 'RAP+RVP+LVP(Endo_Mid_Anterior_Wall)_100',\
                 'RAP+RVP+LVP(Endo_Mid_Lateral_Wall)_100',\
                 'RAP+RVP+LVP(Epi_Mid_Anterior_Wall)_100']
    print(case_name_1,case_name_2,case_name_3,case_name_4,case_name_5,case_name_6,case_name_7)
    case_name=(input('what is the case name for file (if none or not known just type none): '))
    if case_name.lower().strip(' ') == 'none':
        case_name=['']
    case_type = input('it is ECHO (w/testprocess) or MRI(w/otestprocess) (no space): ')
    if case_type.lower().strip() == 'echo':
            ## Check saved file name whether is '_' or ' '
            blank='_'
            basedirectory = "./"+str(pignumber)+extra+"/"
            outdirectory = "echo_volume_measurements/"
            '''filename = "PIG"+blank+str(pignumber)+"_PIG "+str(pignumber)+"_" '''
            
            filename = "PIG "+txt_filename      #the txt file name and be careful about the space in the the txt file name
            xlsfilename=str(pignumber)+'_strain.xlsx'
            
            if not os.path.exists(outdirectory):
                os.makedirs(outdirectory)
            
            for typename in case_name:
            	data=[]
            	with open(filename+typename+".txt") as f:
            	    for line in f:
            	
            	
            	        a1_test4 = [elt.strip() for elt in line.split(',')]
            	        data.append(a1_test4)
                    
            	t = np.asfarray(data[832][1:-1],float)
            	Vol=np.asfarray(data[833][1:-1],float)
            	Mass=np.asfarray(data[38][1],float)
            	SDI=np.asfarray(data[39][1],float)
            	
            	#LL_results=np.ones([16,2])
                #print(data[662][1])
            	#for ii in range(0,6):
            		#LL_results[ii][0]=np.asfarray(data[662+ii][1],float)
            
            	tan_results=np.ones([16,2])
            	for ii in range(0,6):
            		tan_results[ii][0]=np.asfarray(data[490+ii][1],float)
            		tan_results[ii][1]=np.asfarray(data[519+ii][1],float)
            	for ii in range(6,12):
            		tan_results[ii][0]=np.asfarray(data[491+ii][1],float)
            		tan_results[ii][1]=np.asfarray(data[520+ii][1],float)
            	for ii in range(12,16):
            		tan_results[ii][0]=np.asfarray(data[492+ii][1],float)
            		tan_results[ii][1]=np.asfarray(data[521+ii][1],float)
                
                
            	CC_results=np.ones([16,2])
            	for ii in range(0,6):
            		CC_results[ii][0]=np.asfarray(data[576+ii][1],float)
            		CC_results[ii][1]=np.asfarray(data[605+ii][1],float)
            	for ii in range(6,12):
            		CC_results[ii][0]=np.asfarray(data[577+ii][1],float)
            		CC_results[ii][1]=np.asfarray(data[606+ii][1],float)
            	for ii in range(12,16):
            		CC_results[ii][0]=np.asfarray(data[578+ii][1],float)
            		CC_results[ii][1]=np.asfarray(data[607+ii][1],float)
            	
            	LL_results=np.ones([16,2])
            	for ii in range(0,6):
            		LL_results[ii][0]=np.asfarray(data[662+ii][1],float)
            		LL_results[ii][1]=np.asfarray(data[691+ii][1],float)
            	for ii in range(6,12):
            		LL_results[ii][0]=np.asfarray(data[663+ii][1],float)
            		LL_results[ii][1]=np.asfarray(data[692+ii][1],float)
            	for ii in range(12,16):
            		LL_results[ii][0]=np.asfarray(data[664+ii][1],float)
            		LL_results[ii][1]=np.asfarray(data[693+ii][1],float)
            
            	RR_results=np.ones([16,2])
            	for ii in range(0,6):
            		RR_results[ii][0]=np.asfarray(data[748+ii][1],float)
            		RR_results[ii][1]=np.asfarray(data[777+ii][1],float)
            	for ii in range(6,12):
            		RR_results[ii][0]=np.asfarray(data[749+ii][1],float)
            		RR_results[ii][1]=np.asfarray(data[778+ii][1],float)
            	for ii in range(12,16):
            		RR_results[ii][0]=np.asfarray(data[750+ii][1],float)
            		RR_results[ii][1]=np.asfarray(data[779+ii][1],float)
            	
            	
            	np.savetxt(outdirectory+'_Tangential'+'_data.txt',tan_results,delimiter=',')
            	np.savetxt(outdirectory+'_cc'+'_data.txt',CC_results,delimiter=',')
            	np.savetxt(outdirectory+'_rr'+'_data.txt',RR_results,delimiter=',')
            	np.savetxt(outdirectory+'_ll'+'_data.txt',LL_results,delimiter=',')
            	np.savetxt(outdirectory+'_mass_SDI_data.txt',(Mass,SDI),delimiter=',')
            
            	E_tang=np.ones([16,len(t)])
            	for i in range(0,16):
            		E_tang[i][:]=(np.asfarray(data[980+i][1:-1],float))
            	
            	
            	E_cc=np.ones([16,len(t)])
            	for i in range(0,16):
            		E_cc[i][:]=(np.asfarray(data[1005+i][1:-1],float))
            	
            	E_ll=np.ones([16,len(t)])
            	for i in range(0,16):
            		E_ll[i][:]=(np.asfarray(data[1030+i][1:-1],float))
            	
            	E_rr=np.ones([16,len(t)])
            	for i in range(0,16):
            		E_rr[i][:]=(np.asfarray(data[1055+i][1:-1],float))
            	
            	
            
            	plt.figure(figsize=(15,10))
            	plt.plot(t, Vol)
            	plt.xlabel('Time(ms)')
            	plt.ylabel('Volume(ml)')
            	plt.savefig(outdirectory+str(pignumber)+"_volume_"+"_"+typename+".png")
            	np.savez(outdirectory+"Vol_"+typename+".npz", t=t, volume_array=Vol)
            	sio.savemat(outdirectory+"Vol_"+typename+".mat", {'t':t, 'volume_array':Vol})
            
            
            	##tangential
            	plt.figure(figsize=(15,10))
            	for j in range(16):
            		plt.plot(t,E_tang[j][:],label=str(j+1))
            		plt.legend()
            	plt.xlabel('Time(ms)')
            	plt.ylabel('Strain(%)')
            	plt.savefig(outdirectory+str(pignumber)+"_E_tangential_"+"_"+typename+".png")
            	np.savez(outdirectory+"E_tangential_"+typename+".npz", t=t, strain=E_tang)
            	sio.savemat(outdirectory+"E_tangential_"+typename+".mat", {'t':t, 'strain':E_tang})
            	
            	
            	##circumferential
            	plt.figure(figsize=(15,10))
            	for j in range(16):
            		plt.plot(t,E_cc[j][:],label=str(j+1))
            		plt.legend()
            	plt.xlabel('Time(ms)')
            	plt.ylabel('Strain(%)')
            	plt.savefig(outdirectory+str(pignumber)+"_E_cc_"+"_"+typename+".png")
            	np.savez(outdirectory+"E_cc_"+typename+".npz", t=t, strain=E_cc)
            	sio.savemat(outdirectory+"E_cc_"+typename+".mat", {'t':t, 'strain':E_cc})
            	
            	
            	#longitudinal
            	plt.figure(figsize=(15,10))
            	for j in range(16):
            		plt.plot(t,E_ll[j][:],label=str(j+1))
            		plt.legend()
            	plt.xlabel('Time(ms)')
            	plt.ylabel('Strain(%)')
            	plt.savefig(outdirectory+str(pignumber)+"_E_ll_"+"_"+typename+".png")
            	np.savez(outdirectory+"E_ll_"+typename+".npz", t=t, strain=E_ll)
            	sio.savemat(outdirectory+"E_ll_"+typename+".mat", {'t':t, 'strain':E_ll})
            	
            	
            	##Radial
            	plt.figure(figsize=(15,10))
            	for j in range(16):
            		plt.plot(t,E_rr[j][:],label=str(j+1))
            		plt.legend()
            	plt.xlabel('Time(ms)')
            	plt.ylabel('Strain(%)')
            	plt.savefig(outdirectory+str(pignumber)+"_E_rr_"+"_"+typename+".png")
            	np.savez(outdirectory+"E_rr_"+typename+".npz", t=t, strain=E_rr)
            	sio.savemat(outdirectory+"E_rr_"+typename+".mat", {'t':t, 'strain':E_rr})
            
            	print (str(SDI)+'%')
            	print (str(Mass)+'g')
    
    ################################################################
    #interpolation and output Strain_base file
            Ecc = np.load(outdirectory+"E_cc_.npz")["strain"]           #E_cc_strain
            Ecc_t = np.load(outdirectory+"E_cc_.npz")["t"]
            Ell = np.load(outdirectory+"E_ll_.npz")["strain"]           #E_ll_strain
            Ell_t = np.load(outdirectory+"E_ll_.npz")["t"]
            
                
            min_Ecc_tpt_arr = [np.argmin(abs(Ecc_)) for Ecc_ in Ecc]
            min_Ell_tpt_arr = [np.argmin(abs(Ell_)) for Ell_ in Ell]
            
            min_tpt_arr = np.concatenate((min_Ecc_tpt_arr, min_Ell_tpt_arr))
            EDtpt = np.bincount(min_tpt_arr).argmax()
            Ecc_new = np.concatenate((Ecc[:,EDtpt:len(Ecc[0])], Ecc[:,0:EDtpt]), axis=1)
            Ell_new = np.concatenate((Ell[:,EDtpt:len(Ell[0])], Ell[:,0:EDtpt]), axis=1)
            Eccfunc = interp1d(Ecc_t - Ecc_t[0], Ecc_new, kind='linear', fill_value="extrapolate")
            Ellfunc = interp1d(Ell_t - Ell_t[0], Ell_new, kind='linear', fill_value="extrapolate")
            
            maxtpt = max(Ecc_t) - Ecc_t[0]
            E_tpt_arr = np.linspace(0, maxtpt, data_number)     # process identical data point for pressure, volume, and strain
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
            plt.savefig("Ecc.png")          #image for Ecc_strain
            
            
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
            plt.savefig("Ell.png")          #image for Ell_strain
            
            data_1 = scipy.io.loadmat('echo_volume_measurements/'+ 'Vol_')
            x_1=data_1['t']
            y_1=data_1['volume_array']
            t_1 = np.array(x_1[0])/1000                 #changing time from millisecond to second
            volume_array_1 = np.array(y_1[0])
            
            # Output interpolated strain to npz for optimization
            np.savez("Strain_Base.npz", \
                     tpt = E_tpt_arr,\
                     Ecc = Ecc_fit_arr,\
                     Ell = Ell_fit_arr)
    if case_type.lower().strip() == 'mri':
        data=[]
        x_1=[]
        y_1=[]
        VInterp=[]
        with open("LV_Volume.txt") as f:
            for line in f:
                if line.strip() != '':
                    a1_test4 = [elt.strip() for elt in line.split(',')]
                    if (a1_test4[0][1]) != '.':
                        x_1.append(int(a1_test4[0][0:2])*0.04408)   #time in second
                        y_1.append(a1_test4[0][9:])
                    else:
                        x_1.append(int(a1_test4[0][0])*0.04408)     #time in second
                        y_1.append(a1_test4[0][9:])
        t_1 = x_1
        volume_array_1 =(y_1)
        
    data_0 = scipy.io.loadmat(P_filename)
    
    x_0=data_0['t']
    y_0=data_0['LVP']
    
    t_0 = np.array(x_0[0])
    LVP_0 = np.array(y_0[0])

    print('The pressure mat points: t='+ str(len(t_0)) + ' and LVP='+ str(len(LVP_0)))
    print('The volume mat points: t='+ str(len(t_1)) + ' and Volume='+ str(len(volume_array_1)))
    
    ######interpolation
    if max(t_0)>= max(t_1):                     #the size of pressrue and volume is not equal, we need to match them in same time point
        max_t = max(t_1)
    if max(t_0)<=max(t_1):
        max_t = max(t_0)
    #intepolation for pressure and volume
    Pfx = interp1d(t_0, LVP_0, kind = 'linear',fill_value="extrapolate")
    Vfx = interp1d(t_1, volume_array_1, kind = 'linear',fill_value="extrapolate")
    xInterp = np.linspace(0, max_t, data_number)    # process identical data point for pressure, volume, and strain
    xInterp_1 = np.linspace(0, max(t_0), data_number)      # process identical data point for pressure, volume, and strain
    yInterp = Vfx(xInterp)
    LVP_1 = Pfx(xInterp_1)
    t_1 = np.linspace(0,t_1[len(t_1)-1],data_number)
    
    plt.plot(xInterp,yInterp,'o-', markersize=ms)
    plt.title("NEW"+ str(pignumber))
    plt.xlabel("Time")
    plt.ylabel("LVP")
    plt.show()
    print('The pressure mat points after interpolation: t='+str(len(xInterp))+' and LVV='+str(len(yInterp)))
    t = xInterp
    P = LVP_1
    yInterp = np.array(yInterp)
    yInterp -= 4.96
    yInterp = yInterp.tolist()
    V = yInterp
    #print(len(t_1),len(P),len(V))
    n=0
    with open(txt_filename+'.csv','w',newline='') as f:       #CSV file name
        fieldnames = ['t','P','V']
        thewriter = csv.DictWriter(f, fieldnames=fieldnames)
        while n!= int(data_number):
            thewriter.writerow({'t':t_1[n],'P':P[n],'V':V[n]})  #check length for each list or numpy   -13.217
            n+=1



#################################################################### 
#postprocess
try:
    tpt = xInterp
    LVP = LVP_1
    LVV = yInterp
    V_t = xInterp
    V = yInterp
except:
    PVdata = pd.read_csv(txt_filename+'.csv')
    tpt = np.array(PVdata.values[:,0])
    LVP = np.array(PVdata.values[:,1])
    LVV = np.array(PVdata.values[:,2])
    V_t = np.array(PVdata.values[:,0])
    V = np.array(PVdata.values[:,2])
res = 1
# Interpolate PV data
LVPfunc = interp1d(tpt - tpt[0], LVP, kind='linear', fill_value="extrapolate")
LVVfunc = interp1d(tpt - tpt[0], LVV, kind='linear', fill_value="extrapolate")
maxtpt = max(tpt) - tpt[0]

Vfunc = interp1d(V_t - V_t[0], V, kind='linear', fill_value="extrapolate")
maxt = max(V_t) - V_t[0]
V_tpt_arr = np.linspace(0, maxt, res*len(V))
V_fit_arr = Vfunc(V_tpt_arr)

PV_tpt_arr = np.linspace(0, maxtpt, res*len(V))
LVP_fit_arr = LVPfunc(PV_tpt_arr)
LVV_fit_arr = LVVfunc(PV_tpt_arr)


print("Number of fitted data points = ", len(PV_tpt_arr))

plt.figure(1)
plt.plot(LVV, LVP, '-*', label="Expt", markersize=ms)
plt.plot(LVV_fit_arr, LVP_fit_arr, 'o', label="Fitted", markersize=ms)
plt.savefig("PVloop.png")

plt.figure(2)
plt.plot(tpt, LVP, '-*', label="Expt", markersize=ms)
plt.plot(PV_tpt_arr, LVP_fit_arr, 'o', label="Fitted", markersize=ms)
plt.savefig("Pwave.png")

plt.figure(3)
plt.plot(tpt, LVV, '-*', label="Expt", markersize=ms)
plt.plot(PV_tpt_arr, LVV_fit_arr, 'o', label="Fitted", markersize=ms)
plt.savefig("Vwave.png")


# Output interpolated PV loop to npz for optimization
np.savez( "PVloop_RV100.npz", \
        LVP = LVP_fit_arr,\
        LVV = LVV_fit_arr,\
        V = V_fit_arr)
