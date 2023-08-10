# Preparation of Data

## Postprocessing strain and PV data
The pressure, volume and strain data are postprocessed for optimization using the file
[`postprocessPVstrain.py`](postprocessPVstrain.py).

The python script converts the pressure, volume waveforms into `PVloop.npz` and strains into `Strain.npz` that can be used as input into the optimization.

In that script, pressure and volume data are read and can be interpolated to create a higher resolution data. The script will also read the strain and interpolate it so that it has the **same resolution** as that of the pressure volume data.

## Previewing of fiber orientation and strain basis orientation
A python script [`postprocessmesh.py`](postprocessmesh.py) is given to output the circumferential, longitudinal, fiber, cross-fiber and sheet directions, as well as material regions, as `.pvd` file.

An example is given using a data set in directory [`data_PIG2455`](data_PIG2455).
