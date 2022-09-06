#Reorganization of Pressure Volume Loop and Strain Graphing from HeArt Project 

The org.py is processing the pressure, volume, and strain value. The example is provided. There are two process methods 3D ECHO (with strain value) (in special format) and MRI (Magnetic Resonance Images) (no strain value). The ‘PIG 393485.txt’ is the echo file and ‘LV_Volume.txt’ is the MRI file. Also, the ‘LVpressure.mat’ is pressure data in a full cycle. When running the echo parts, the program will generate a folder with volume and multiple different strains in it.  

After the first successful running, it can generate a csv file which has pressure, volume, and time in it. If you would like to make some data changes, you can change inside, and the code will use the data in the csv file to run in second time.  

Recent edited 294-297. It compares the time of pressure and volume.Then, it takes the shorter period for interpolation. Shown here:

if max(t_0)>= max(t_1): 
     max_t = max(t_1)
if max(t_0)<max(t_1):
     max_t = max(t_0)
