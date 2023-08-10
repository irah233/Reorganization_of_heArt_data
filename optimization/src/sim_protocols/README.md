# Simulation protocols

The simulation protocols folder contains 2 scripts for the forward run [`heArtsolver_clean3.py`](./heArtsolver_clean3.py) and the optimization [`optimization3.py`](./optimization3.py).

## Forward Run
The forward run script solve the forward problem of cardiac mechanics to determine the LV pressures (`Ptargets`) given  prescribed LV volumes (`Vtargets`) and contractility (`Tmax_ctrls`). The contractility is given as a list of function in the `DG0` space (i.e., constant in each element) that has length as that in `Ptargets` and `Vtargets`.

The forward run will also compute the objective functional `j` based on the absolute difference between the computed LV pressures (`pendo_`) and the target LV pressures (`Ptargets`), and the target strains in % (`eCCtargets` and `eLLtargets`) if present . The objective functional `j` will be used for optimization.

The objective function is given by
<!--
j_1 = K_{lvp}\sum_{i=0}^{N} \left( \frac{(LVP_i - LVP_{target,i})}{LVP_{target,i}} \right )^2 \\
\\
j_2 = K_{reg}\sum_{i=0}^{N} \frac{1}{V}\int_{\Omega} \left(\frac{(T_{max,i} - \bar{T}_{max,i})}{\bar{T}_{max,i}} \right )^2 \hspace{3pt} dV
\mbox{where} \hspace{0.2cm} \bar{T}_{max,i} = \int_{Omega}T_{max,i} \hspace{3pt}dV \\
\\
j_3 = K_{str} \sum_{i=0}^{N} \sum_{p=0}^{N_{mat}}\frac{1}{V_{p}} \left(\int_{\Omega_{p}} e_{cc,i} - e_{cc,target,i}\right )^2 \\
\\
j_4 = K_{str} \sum_{i=0}^{N} \sum_{p=0}^{N_{mat}}\frac{1}{V_{p}} \left(\int_{\Omega_{p}} e_{ll,i} - e_{ll,target,i}\right )^2 \\
\\
\\
j = j_1 + j_2 + j_3 + j_4
-->

![equation](https://latex.codecogs.com/gif.latex?j_1%20%3D%20K_%7Blvp%7D%5Csum_%7Bi%3D0%7D%5E%7BN%7D%20%5Cleft%28%20%5Cfrac%7B%28LVP_i%20-%20LVP_%7Btarget%2Ci%7D%29%7D%7BLVP_%7Btarget%2Ci%7D%7D%20%5Cright%20%29%5E2%20%5C%5C%20%5C%5C%20j_2%20%3D%20K_%7Breg%7D%5Csum_%7Bi%3D0%7D%5E%7BN%7D%20%5Cfrac%7B1%7D%7BV%7D%5Cint_%7B%5COmega%7D%20%5Cleft%28%5Cfrac%7B%28T_%7Bmax%2Ci%7D%20-%20%5Cbar%7BT%7D_%7Bmax%2Ci%7D%29%7D%7B%5Cbar%7BT%7D_%7Bmax%2Ci%7D%7D%20%5Cright%20%29%5E2%20%5Chspace%7B3pt%7D%20dV%20%5Cmbox%7Bwhere%7D%20%5Chspace%7B0.2cm%7D%20%5Cbar%7BT%7D_%7Bmax%2Ci%7D%20%3D%20%5Cint_%7BOmega%7DT_%7Bmax%2Ci%7D%20%5Chspace%7B3pt%7DdV%20%5C%5C%20%5C%5C%20j_3%20%3D%20K_%7Bstr%7D%20%5Csum_%7Bi%3D0%7D%5E%7BN%7D%20%5Csum_%7Bp%3D0%7D%5E%7BN_%7Bmat%7D%7D%5Cfrac%7B1%7D%7BV_%7Bp%7D%7D%20%5Cleft%28%5Cint_%7B%5COmega_%7Bp%7D%7D%20e_%7Bcc%2Ci%7D%20-%20e_%7Bcc%2Ctarget%2Ci%7D%5Cright%20%29%5E2%20%5C%5C%20%5C%5C%20j_4%20%3D%20K_%7Bstr%7D%20%5Csum_%7Bi%3D0%7D%5E%7BN%7D%20%5Csum_%7Bp%3D0%7D%5E%7BN_%7Bmat%7D%7D%5Cfrac%7B1%7D%7BV_%7Bp%7D%7D%20%5Cleft%28%5Cint_%7B%5COmega_%7Bp%7D%7D%20e_%7Bll%2Ci%7D%20-%20e_%7Bll%2Ctarget%2Ci%7D%5Cright%20%29%5E2%20%5C%5C%20%5C%5C%20%5C%5C%20j%20%3D%20j_1%20&plus;%20j_2%20&plus;%20j_3%20&plus;%20j_4)

Inputs are given as follows:
- `Tmax_ctrls`:   List of contractility given as function in the real space
- `Vtargets`: Target volumes
- `Ptargets`: Target pressures for computing the objective functional
- `mesh`: LV mesh
- `eCCtargets`: (Optional) Target regional circumferential strains for computing the objective functional
- `eLLtargets`: (Optional) Target regional longitudinal strains for computing the objective functional
- `facetboundaries`: facet labels for imposing BCs
- `edgeboundaries`: edge labels for imposing BCs
- `f0`: fiber unit vectors
- `s0`: sheet unit vectors
- `n0`: sheet-normal unit vectors
- `eC0`: circumferential unit vectors
- `eL0`: longitudinal unit vectors
- `eR0`: radial unit vectors
- `matid`: material IDs given as a DG0 function
- `SimDet`: Simulation details containing material parameters and other options.
- `IODet`: I/O details

## Optimization

The optimization script is used to find the list of contractility `Tmax_ctrls` that best matches the PV loop given by `Vtargets`, `Ptargets`, `eCCtargets` (if available) and `eLLtargets` (if available).

It takes the parameter sets given in `IODet` and `SimDet`.

In `IODet`, the following parameters are used:

- `casename`: case name associated with the HDF5 mesh files.
- `directory_me`: path of folder containing mesh.
- `outputfolder`: user prescribed output directory for the simulation
- `folderName`: user prescribed output folder for the simulation
- `caseID`: user prescribed case ID for the simulation
- `PVloop_data_file`: data file containing the pressures and volumes as variable `LVP` and `LVV` in .npz format.
- `Strain_data_file`: data file containing the eCC and eLL as variable `Ecc` and `Ell` in .npz format. These variables are given as a numpy 2D array where row is associated with region and column is associated with time point.
- `outputfolder`: output folder for dumping of optimization results
- `Initfile`: (optional) initialization of the optimization variable `Tmax`. If not given, the code will assume a homogeneous value specified by `init_opt_val` in `SimDet`.

In `SimDet`, the following parameters are used:

- `GiccioneParams`: Parameters for the mechanics model (See [`GiccioneParams`](#giccioneparams))
- `nLoadSteps`: Number of load steps during the preloading from unloaded geometry to EDV (given as first entry of `Vtargets`).
- `UpperBd`: Upper bound for the optimization variables
- `LowerBd`: Lower bound for the optimization variables
- `Nintermediate_stp`: Intermediate steps to run between each time point in the data set
- `Kreg`: Penalty factor associated with the regularization functional
- `Kstr`: Penalty factor associated with the strain functional
- `Klvp`: Penalty factor associated with the PV functional
- `init_opt_val`: Initialization value for optimization variables
- `topid`: Basal facet ID
- `LVendoid`: LV endocardium facet ID
- `RVendoid`: RV endocardium facet ID (not used)
- `epiid`: Epicardium facet ID
- `abs_tol`: Absolute tolerance for newton iteration
- `rel_tol`: Relative tolerance for newton iteration
- `max_newton_iter`: Maximum number of newton iteration

#### `GiccioneParams`
- `Passive model`: Name of the passive model given in python dictionary `{Name: "<model name>"}`. The model implemented so far include only Guccione (`Guccione`).
- `Passive params`: Parameters of the passive model given in python dictionary.
- `Active model`: Name of the active model given in python dictionary `{Name: "<model name>"}`. The model implemented so far include Time-varying (`Time-varying`) with time dependency removed.
- `Active params`: Parameters of the active model given in python dictionary.

An demo example is given in [`optimizePIG2455_2.py`](../demo/optimizePIG2455_2.py) where it uses the mesh, PV and strain data given in the [`data directory`](../data/data_PIG2455). The script optimize the time-course of `Tmax` (given as a function in DG0 space) to fit the PV loops and strain. The optimized solution is given in [`Optimal_Klvp1e5_Kstr1e-2.hdf5`](../demo/Optimal_Klvp1e5_Kstr1e-2.hdf5). Specifying this file as an initial guess using the parameter `Initfile` will result in the optimization converging in 1 step.

**Note**: The code will generate a `Log.hdf5` file in the output directory at each optimization step. This file can be used as a restart for the optimization by specifying its path in `Initfile` in `IODet`. This is useful in the case when the newton iterations do not converge because the change in `Tmax` is too large in each load step. In that case, one can restart the optimization using more intermediate loading step specified by `Nintermediate_stp`.

*Optimized Tmax time course*

![Tmax_time_course](../../figures/Tmax.png)

*Optimized PV loop (dots: data line: fitted results)*

![PV loop](../../figures/PVloop.png)

*Optimized eCC (dots: data line: fitted results)*

![PV loop](../../figures/Ecc.png)

*Optimized eLL (dots: data line: fitted results)*

![PV loop](../../figures/Ell.png)
