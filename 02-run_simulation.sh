#!/bin/bash

for i in 0.5 1.0 2.0 3.0 4.0 5.0
do
	cd ${i}M/
	# Equilibration
	python ../common/openmm_simulation.py omm_equil.inp
	
	# Production
	python ../common/openmm_simulation.py omm_osmotic.inp
	cd ../
done
