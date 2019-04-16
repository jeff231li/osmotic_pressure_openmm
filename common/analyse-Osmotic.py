import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
import sys

ncfile    = netcdf.Dataset(sys.argv[1])

Nframes   = ncfile.dimensions['frames'].size
Nensemble = ncfile.dimensions['ensemble'].size
Positions = ncfile.variables['positions'][:,:,:]

Xlen      = ncfile.variables['xlen'][:][0]
Ylen      = ncfile.variables['ylen'][:][0]
Zmax      = ncfile.variables['zmax'][:][0]
Zmin      = ncfile.variables['zmin'][:][0]
Kwall     = ncfile.variables['kz'][:][0]

Area      = Xlen*Ylen*1e-20
Osmotic   = np.zeros((Nensemble, 1), dtype=np.float32)

for i in range(Nensemble):
	MatZ       = Positions[:,:,i].flatten()*10
	
	Izmax      = MatZ >= Zmax
	Izmin      = MatZ <= Zmin
	
	UpperWall  = Kwall*abs(MatZ[Izmax] - Zmax)
	LowerWall  = Kwall*abs(MatZ[Izmin] - Zmin)
	
	Fupper     = np.sum(UpperWall)/Nframes
	Flower     = np.sum(LowerWall)/Nframes
	Fwall      = (Fupper+Flower)/2.0 * 69.5e-12
	
	Osmotic[i] = Fwall/Area/1.01e5
	
	#print("Osmotic Pressure: %.1f Bar" % (Osmotic[i]))

ave = np.mean(Osmotic)
std = np.std(Osmotic)/np.sqrt(Nensemble)
print("%.1f +- %.1f Bar" % (ave, std))
