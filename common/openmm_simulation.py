import simtk.unit as unit
import simtk.openmm as openmm
import simtk.openmm.app as app
from common.openmm_readinputs import read_config, create_system
import netCDF4 as netcdf
import time
import sys

# User Input
inp = read_config(sys.argv[1])
inp.initialize()
sys.stdout = open('%s.log' % inp.output_name, 'w')

# Hardware Properties
platform = openmm.Platform.getPlatformByName(inp.Hardware_type)
properties = {}
if inp.Hardware_type == 'CUDA':
    properties = {'CudaDeviceIndex': inp.Hardware_gpu_idx,
                  'CudaPrecision'  : inp.Hardware_gpu_prec}

# Create OpenMM System
system, topology, coordinates = create_system(inp)

# Wall restraints
if inp.POT_wall == 'on':
    print("Flat-bottom wall setup")
    sys.stdout.flush()
    z_wall = openmm.CustomExternalForce("(k_z/2.0)*(max(0, z-z_max)^2 + min(0, z-z_min)^2)")
    z_wall.addGlobalParameter("k_z",   inp.POT_kz * unit.kilocalories_per_mole/unit.angstroms**2)
    z_wall.addGlobalParameter("z_max", inp.POT_zmax * unit.angstroms)
    z_wall.addGlobalParameter("z_min", inp.POT_zmin * unit.angstroms)
    idx = 0
    ion_index = []
    for line in open(inp.POT_atoms):
        if line.startswith("ATOM"):
            dummy = line.split()
            tags = float(line[inp.POT_colA:inp.POT_colB])
            if tags == inp.POT_tag:
                z_wall.addParticle(idx, [])
                ion_index.append(idx)
            idx += 1
    n_ions = len(ion_index)
    system.addForce(z_wall)
    z_wall.setForceGroup(10)

# Thermostat and Barostat
integrator = openmm.LangevinIntegrator(
    inp.temperature * unit.kelvin,
    inp.fric_coeff / unit.picosecond,
    inp.dt * unit.femtosecond,
)

if inp.npt == 'on':
    print("Running with anisotropic Monte Carlo barostat")
    sys.stdout.flush()
    barostat = openmm.MonteCarloAnisotropicBarostat(
        (inp.npt_val, inp.npt_val, inp.npt_val) * unit.bar,
        inp.temperature * unit.kelvin,
        False,
        False,
        True,
        inp.npt_freq
    )
    system.addForce(barostat)

# Simulation Object
print("Creating Simulation object")
sys.stdout.flush()
simulation = app.Simulation(topology.topology, system, integrator, platform, properties)
simulation.context.setPositions(coordinates.positions)

if inp.restart_name:
    print("Reading positions and velocities from %s" % inp.restart_name+'.xml')
    sys.stdout.flush()
    with open(inp.restart_name+'.xml', 'r') as f:
        simulation.context.setState(openmm.XmlSerializer.deserialize(f.read()))

# Minimization Phase
if inp.min_run == 'on':
    print("\nMinimizing Energy")
    sys.stdout.flush()
    simulation.minimizeEnergy(
        maxIterations=inp.min_steps,
        tolerance=inp.min_tol * unit.kilocalories_per_mole
    )
    simulation.context.setVelocitiesToTemperature(inp.temperature * unit.kelvin)

    print("Energy after minimization")
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

# Equilibration Phase
if inp.MD_run == 'on':
    print("\nEquilibrating for %i steps" % inp.MD_equil)
    sys.stdout.flush()
    if inp.dcdreporter == 'on':
        simulation.reporters.append(app.DCDReporter('%s.dcd' % inp.output_name, inp.nstdcd))

    if inp.reporter == 'on':
        simulation.reporters.append(
            app.StateDataReporter(
                '%s_reporter.log' % inp.output_name,
                inp.nstout,
                step=True,
                kineticEnergy=True,
                potentialEnergy=True,
                totalEnergy=True,
                temperature=True,
                totalSteps=inp.MD_equil,
                progress=True,
                remainingTime=True,
                speed=True,
                separator="\t",
            )
        )

    # Equilibration
    t1 = time.time()
    simulation.step(inp.MD_equil)
    t2 = time.time()

    # Housekeeping
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
    print("\nEquilibration completed in %.2f" % (t2-t1))

    print("\nSaving Coordinates")
    sys.stdout.flush()
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(simulation.topology, positions, open('%s.pdb' % inp.output_name, 'w'))
    simulation.saveCheckpoint('%s.chk' % inp.output_name)

    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open('%s.xml' % inp.output_name, 'w') as f:
        f.write(openmm.XmlSerializer.serialize(state))

# Production - Ensemble
if inp.OP_run == 'on':
    # NetCDF Output file
    ncfile = netcdf.Dataset(inp.OP_output, 'w', format='NETCDF4')
    ncfile.createDimension('harmonic', None)
    ncfile.createDimension('frames', 0)
    ncfile.createDimension('atoms', n_ions)
    ncfile.createDimension('ensemble', inp.OP_ensembles)
    ncfile.createVariable('positions', 'f', ('frames', 'atoms', 'ensemble'))
    ncfile.createVariable('x_len', 'f', ('harmonic',))
    ncfile.createVariable('y_len', 'f', ('harmonic',))
    ncfile.createVariable('z_max', 'f', ('harmonic',))
    ncfile.createVariable('z_min', 'f', ('harmonic',))
    ncfile.createVariable('k_z',   'f', ('harmonic',))

    ncfile.variables['x_len'][:] = inp.Box_xlen
    ncfile.variables['y_len'][:] = inp.Box_ylen
    ncfile.variables['z_max'][:] = inp.POT_zmax
    ncfile.variables['z_min'][:] = inp.POT_zmin
    ncfile.variables['k_z'][:] = inp.POT_kz

    if inp.dcdreporter is 'on' or inp.dcdreporter is 'yes':
        simulation.reporters.append(app.DCDReporter('%s.dcd' % inp.output_name, inp.nstdcd))

    if inp.reporter is 'on' or inp.reporter is 'yes':
        simulation.reporters.append(
            app.StateDataReporter(
                '%s_reporter.log' % inp.output_name,
                inp.nstout,
                step=True,
                kineticEnergy=True,
                potentialEnergy=True,
                totalEnergy=True,
                temperature=True,
                totalSteps=inp.OP_total,
                progress=True,
                remainingTime=True,
                speed=True,
                separator="\t",
              )
        )

    for i in range(inp.OP_ensembles):
        print("\nRunning ensemble simulation: %i" % (i+1))
        start = time.time()

        # Initialize
        print("\tLoading initial structure...")
        simulation.loadCheckpoint('%s.chk' % inp.restart_name)
        print("\tInitial system energy ", simulation.context.getState(getEnergy=True).getPotentialEnergy())
        print("\tSetting random velocities...")
        simulation.context.setVelocitiesToTemperature(inp.temperature * unit.kelvin)
        sys.stdout.flush()

        # Equilibration
        print("\tEquilibrating for %i steps" % inp.OP_equil)
        sys.stdout.flush()
        simulation.step(inp.OP_equil)

        # Production
        print("\tRunning Production for %i steps" % inp.OP_production)
        sys.stdout.flush()
        for j in range(inp.OP_nblocks):
            simulation.step(inp.OP_freq)
            positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
            ncfile.variables['positions'][j, :, i] = positions[n_ions, -1]

        # Finalize
        end = time.time()
        print("\tTotal time taken: %.2f sec" % (end - start))

    ncfile.close()
