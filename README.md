#Reorganization of Pressure Volume Loop and Strain Graphing from HeArt Project 

Simulator of mechanics in the heart based on [FEniCS](https://fenicsproject.org/) library.

<!-- TOC -->
  - [Installation and Running the code](#installation-and-running-the-code)
  - [Organization of the code](#organization-of-the-code)
  - [Simulation protocols](#simulation-protocols)

<!-- /TOC -->


### Installation and Running the code
A singularity "build" [file](./SingularitY/Singularity_fenics2017_msucompbiomechlab) is provided that will install all the libraries required to run the code.

1. Install singularity by following the instruction in [here](https://sylabs.io/guides/3.5/admin-guide/installation.html)

2. Build a singularity container using the "build" [file](./SingularitY/Singularity_fenics2017_msucompbiomechlab) with
```
sudo singularity build <container_name>.img Singularity_fenics2017_msucompbiomechlab
```

3. Once the container is built, you can launch the Singularity container by
```
 singularity run <container_name>.img
```

4. The code can be run within the singularity container. For example, for the code [createLV_refine](./ed_mesh_create/Patient_1/createLV_refine.py)  
```
python createLV_refine.py
```
or in parallel
```
mpirun.mpich -np <# processors> python createLV_refine.py
```
###
The org.py is processing the pressure, volume, and strain value. The example is provided. There are two process methods 3D ECHO (with strain value) (in special format) and MRI (Magnetic Resonance Images) (no strain value). The ‘PIG 393485.txt’ is the echo file and ‘LV_Volume.txt’ is the MRI file. Also, the ‘LVpressure.mat’ is pressure data in a full cycle. When running the echo parts, the program will generate a folder with volume and multiple different strains in it.  

After the first successful running, it can generate a csv file which has pressure, volume, and time in it. If you would like to make some data changes, you can change inside, and the code will use the data in the csv file to run in second time.  


