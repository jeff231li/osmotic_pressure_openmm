#!/bin/bash

for conc in 0.5 1.0 2.0 3.0 4.0 5.0
do
	cd ${i}M/
	python ../common/analyse-Osmotic.py osmotic-ions.nc
	cd ../
done
