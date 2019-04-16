from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from openmm_readinputs import *
import netCDF4 as netcdf
import time
import sys


# User Input
inp = read_config(sys.argv[1])
inp.initialize()
sys.stdout = open('%s.log' % inp.output_name, 'w')

# Hardware
platform = Platform.getPlatformByName(inp.Hardware_type)
if inp.Hardware_type == 'CUDA':
    properties = {'CudaDeviceIndex': inp.Hardware_gpu_idx,
                  'CudaPrecision'  : inp.Hardware_gpu_prec}

# Load CHARMM files
print("Loading structures")
sys.stdout.flush()
psf, pdb, params = load_structure(inp)

# Create System Object
print("Creating System")
sys.stdout.flush()
system = psf.createSystem(params,
                          nonbondedMethod=PME,
                          nonbondedCutoff=inp.cutoff,
                          switchDistance=inp.switch,
                          constraints=inp.constraint,
                          rigidWater=True)

## WALL RESTRAINTS
#-------------------------------------------------------------------------------#
if inp.POT_wall == 'on':
    print("Flat-bottom restraint setup")
    sys.stdout.flush()
    zForce = CustomExternalForce("(kz/2.0)*(max(0, z-zmax)^2 + min(0,z-zmin)^2)")
    zForce.addGlobalParameter("kz",   inp.POT_kz*kilocalories_per_mole/angstroms**2)
    zForce.addGlobalParameter("zmax", inp.POT_zmax*angstroms)
    zForce.addGlobalParameter("zmin", inp.POT_zmin*angstroms)
    I = 0
    IonsI = []
    for line in open(inp.POT_atoms):
        if line.startswith("ATOM"):
            dummy = line.split()
            tags = float(line[inp.POT_colA:inp.POT_colB])
            if tags == inp.POT_tag:
                zForce.addParticle(I, [])
                IonsI.append(I)
            I += 1
    Nions = len(IonsI)
    system.addForce(zForce)

## THERMOSTAT & BAROSTAT
#-------------------------------------------------------------------------------#
integrator = LangevinIntegrator(inp.temperature*kelvin,
                                inp.fric_coeff/picosecond,
                                inp.dt*femtoseconds)

if inp.npt == 'on':
    print("Running with Monte Carlo barostat")
    sys.stdout.flush()
    barostat = MonteCarloAnisotropicBarostat(
                               (inp.npt_val, inp.npt_val, inp.npt_val)*bar,
                               inp.temperature*kelvin,
                               False,
                               False,
                               True,
                               inp.npt_freq)
    system.addForce(barostat)

## SIMULATION OBJECT
#-------------------------------------------------------------------------------#
print("Creating Simulation object")
sys.stdout.flush()
simulation = Simulation(psf.topology, system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)

if inp.restart_name:
    print("Reading positions and velocities from %s" % inp.restart_name+'.rst')
    sys.stdout.flush()
    with open(inp.restart_name+'.rst', 'r') as f:
        simulation.context.setState(XmlSerializer.deserialize(f.read()))

## MINIMIZATION
# -------------------------------------------------------------------------------#
if inp.min_run == 'on':
    print("\nMinimizing Energy")
    sys.stdout.flush()
    simulation.minimizeEnergy(maxIterations=inp.min_nstep,
                              tolerance=inp.min_tol*kilocalories_per_mole)
    simulation.context.setVelocitiesToTemperature(inp.temperature*kelvin)
    print("Energy after minimization")
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

## EQUILIBRATION
# -------------------------------------------------------------------------------#
if inp.MD_run == 'on':
    print("\nEquilibrating for %i steps" % inp.MD_equil)
    sys.stdout.flush()
    if inp.dcdreporter == 'on':
        simulation.reporters.append(DCDReporter('%s.dcd' % inp.output_name, inp.nstdcd))

    if inp.reporter == 'on':
        simulation.reporters.append(StateDataReporter('%s_reporter.log' % inp.output_name,
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
                                    separator="\t"))
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
    PDBFile.writeFile(simulation.topology, positions, open('%s.pdb' % inp.output_name, 'w'))
    simulation.saveCheckpoint('%s.chk' % inp.output_name)

    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open('%s.rst' % inp.output_name, 'w') as f:
        f.write(XmlSerializer.serialize(state))

## RUN ENSEMBLE SIMULATION
#-------------------------------------------------------------------------------#
if inp.OP_run == 'on':
    # NetCDF Output file
    ncfile = netcdf.Dataset(inp.OP_output, 'w', version='NETCDF4')
    ncfile.createDimension('harmonic', None)
    ncfile.createDimension('frames', 0)
    ncfile.createDimension('atoms', Nions)
    ncfile.createDimension('ensemble', inp.OP_ensembles)
    ncfile.createVariable('positions', 'f', ('frames', 'atoms', 'ensemble'))
    ncfile.createVariable('xlen', 'f', ('harmonic',))
    ncfile.createVariable('ylen', 'f', ('harmonic',))
    ncfile.createVariable('zmax', 'f', ('harmonic',))
    ncfile.createVariable('zmin', 'f', ('harmonic',))
    ncfile.createVariable('kz',   'f', ('harmonic',))

    ncfile.variables['xlen'][:] = inp.Box_xlen
    ncfile.variables['ylen'][:] = inp.Box_ylen
    ncfile.variables['zmax'][:] = inp.POT_zmax
    ncfile.variables['zmin'][:] = inp.POT_zmin
    ncfile.variables['kz'][:] = inp.POT_kz

    if inp.dcdreporter == 'on':
        simulation.reporters.append(DCDReporter('%s.dcd' % inp.output_name, inp.nstdcd))

    if inp.reporter == 'on':
        simulation.reporters.append(StateDataReporter('%s_reporter.log' % inp.output_name,
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
                                    separator="\t"))

    for i in range(inp.OP_ensembles):
        print("\nRunning ensemble simulation: %i" % (i+1))
        start = time.time()

        # Initialize
        print("\tLoading initial structure...")
        simulation.loadCheckpoint('%s.chk' % inp.restart_name)
        print("\tInitial system energy ", simulation.context.getState(getEnergy=True).getPotentialEnergy())
        print("\tSetting random velocities...")
        simulation.context.setVelocitiesToTemperature(inp.temperature*kelvin)
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
            ncfile.variables['positions'][j, :, i] = positions[IonsI, -1]

        # Finalize
        end = time.time()
        print("\tTotal time taken: %.2f sec" % (end - start))

    ncfile.close()

