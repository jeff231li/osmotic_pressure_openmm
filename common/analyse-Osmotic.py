import sys
import numpy as np
import netCDF4 as netcdf

ncfile = netcdf.Dataset(sys.argv[1])

n_frames = ncfile.dimensions['frames'].size
n_ensemble = ncfile.dimensions['ensemble'].size
coordinates = ncfile.variables['positions'][:, :, :]

x_len = ncfile.variables['x_len'][:][0]
y_len = ncfile.variables['y_len'][:][0]
z_max = ncfile.variables['z_max'][:][0]
z_min = ncfile.variables['z_min'][:][0]
k_wall = ncfile.variables['k_z'][:][0]

Area = x_len * y_len * 1e-20
Osmotic = np.zeros((n_ensemble, 1), dtype=np.float32)

for i in range(n_ensemble):
    pos_z = coordinates[:, :, i].flatten() * 10

    i_z_max = pos_z >= z_max
    i_z_min = pos_z <= z_min

    UpperWall = k_wall * abs(pos_z[i_z_max] - z_max)
    LowerWall = k_wall * abs(pos_z[i_z_min] - z_min)

    Fupper = np.sum(UpperWall) / n_frames
    Flower = np.sum(LowerWall) / n_frames
    Fwall = (Fupper + Flower) / 2.0 * 69.5e-12

    Osmotic[i] = Fwall / Area / 1.01e5
    print(f"Ensemble {i+1} -- Pressure: {Osmotic[i]:.1f} Bar")
print('-'*50)

ave = np.mean(Osmotic)
std = np.std(Osmotic) / np.sqrt(n_ensemble)
print(f"Average: {ave:.1f} +- f{std:.1f} Bar")
