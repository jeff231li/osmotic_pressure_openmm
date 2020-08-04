# Info
This repository contains OpenMM scripts to calculate the osmotic pressure of NaCl. The method was developed by Luo and Roux and is described in the article below:

`Luo, Y. & Roux, B. Simulation of Osmotic Pressure in Concentrated Aqueous Salt Solutions. J. Phys. Chem. Lett. 1, 183–189 (2010).`

Basically, once we have the system (i.e. psf and pdb files) we run an initial minimization and equilibration. Then we perform an ensemble of 10-20 simulations to calculate the osmotic pressure.

For the NaCl system, the bash scripts below controls the flow of the simulation. The simulation part can be parallelized over the different systems (i.e. over concentration).
```
sh 01-create_system.sh
sh 02-run_simulation.sh
sh 03-analyse-Osmotic.sh
```
The coordinates of the ions are stored in a NetCDF file, which is analysed with the `analyse-Osmotic.py` python script in the common folder.

These scripts can be used for other molecules and the parameters & options are controlled in the input scripts `omm_equil.inp` and `omm_osmotic.inp`.
