from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np


class _OpenMMOSMOConfig():
    def __init__(self):
        # Hardware
        self.Hardware_type = 'CUDA'
        self.Hardware_gpu_prec = 'mixed'
        self.Hardware_gpu_idx = '0'

        # Input/Output
        self.structure = None
        self.coordinates = None
        self.output_name = None
        self.restart_name = None
        self.param_file = None
        self.nstout = 500
        self.nstdcd = 500
        self.reporter = 'off'
        self.dcdreporter = 'off'

        # MD Integration and CutOff
        self.dt = 2.0
        self.cutoff = 12.0
        self.switch = 10.0
        self.constraint = HBonds

        # NVT
        self.fric_coeff = 1.0
        self.temperature = 300

        # NPT
        self.npt = 'off'
        self.npt_val = 0.0
        self.npt_freq = 25
        self.Box_xlen = 0.0
        self.Box_ylen = 0.0

        # Minimization
        self.min_run = 'off'
        self.min_nstep = 1000
        self.min_tol = 0.1

        # MD Equilibration
        self.MD_run = 'off'
        self.MD_equil = 50000

        # Flat-bottom Potential
        self.POT_wall = 'on'
        self.POT_zmax = 0.0
        self.POT_zmin = 0.0
        self.POT_kz = 0.0
        self.POT_atoms = None
        self.POT_colA = 0
        self.POT_colB = 0
        self.POT_tag = 1.0

        # Osmotic Pressure
        self.OP_run = 'off'
        self.OP_ensembles = 10
        self.OP_equil = 50000
        self.OP_production = 50000
        self.OP_freq = 50
        self.OP_output = None
        self.OP_total = 0
        self.OP_nblocks = 0

    def initialize(self):
        self.OP_total = (self.OP_equil + self.OP_production) * self.OP_ensembles
        self.OP_nblocks = int(self.OP_production / self.OP_freq)

    def readconfig(self, inputFile):
        for line in open(inputFile, 'r'):
            if not line.startswith("#"):
                if line.find("#") >= 0:
                    line = line.split("#")[0]
                dummy = line.strip().split("=")
                if len(dummy) > 1:
                    inp_param = dummy[0].upper().strip()
                    try:
                        inp_value = dummy[1].strip()
                    except:
                        inp_value = None
                    if inp_value:
                        # Hardware
                        if inp_param == 'HARDWARE_TYPE':
                            self.Hardware_type = inp_value
                        if inp_param == 'HARDWARE_GPU_PREC':
                            self.Hardware_gpu_pre = inp_value
                        if inp_param == 'HARDWARE_GPU_IDX':
                            self.Hardware_gpu_idx = inp_value

                        # Input/Output
                        if inp_param == 'STRUCTURE':
                            self.structure = inp_value
                        if inp_param == 'COORDINATES':
                            self.coordinates = inp_value
                        if inp_param == 'OUTPUT_NAME':
                            self.output_name = inp_value
                        if inp_param == 'RESTART_NAME':
                            self.restart_name = inp_value
                        if inp_param == 'PARAM_FILE':
                            self.param_file = inp_value
                        if inp_param == 'NSTOUT':
                            self.nstout = int(inp_value)
                        if inp_param == 'NSTDCD':
                            self.nstdcd = int(inp_value)
                        if inp_param == 'REPORTER':
                            self.reporter = inp_value
                        if inp_param == 'DCDREPORTER':
                            self.dcdreporter = inp_value

                        # MD Integration and CutOff
                        if inp_param == 'TIME_STEP':
                            self.dt = float(inp_value)
                        if inp_param == 'CUTOFF':
                            self.cutoff = float(inp_value) * angstrom
                        if inp_param == 'SWITCHING':
                            self.switch = float(inp_value) * angstrom
                        if inp_param == 'CONSTRAINT':
                            if inp_value == 'NONE':
                                self.constraint = None
                            if inp_value == 'HBONDS':
                                self.constraint = HBonds
                            if inp_value == 'ALLBONDS':
                                self.constraint = AllBonds
                            if inp_value == 'HANGLES':
                                self.constraint = HAngles

                        # NVT
                        if inp_param == 'FRIC_COEFF':
                            self.fric_coeff = float(inp_value)
                        if inp_param == 'TEMPERATURE':
                            self.temperature = float(inp_value)

                        # NPT
                        if inp_param == 'PRESSURE_SIM':
                            self.npt = inp_value
                        if inp_param == 'PRESSURE_VAL':
                            self.npt_val = float(inp_value)
                        if inp_param == 'PRESSURE_FREQ':
                            self.npt_freq = int(inp_value)
                        if inp_param == 'BOX_DIM_XY':
                            temporary = dummy[-1].split()
                            if len(temporary) == 1:
                                self.Box_xlen = float(temporary[0])
                                self.Box_ylen = float(temporary[0])
                            elif len(temporary) == 2:
                                self.Box_xlen = float(temporary[0])
                                self.Box_ylen = float(temporary[1])
                            
                        # Minimization
                        if inp_param == 'MIN_RUN':
                            self.min_run = inp_value
                        if inp_param == 'MIN_NSTEPS':
                            self.min_nstep = int(inp_value)
                        if inp_param == 'MIN_TOL':
                            self.min_tol = float(inp_value)

                        # MD Equilibration
                        if inp_param == 'MD_RUN':
                            self.MD_run = inp_value
                        if inp_param == 'MD_EQUIL':
                            self.MD_equil = int(inp_value)

                        # Flat-bottom Potential
                        if inp_param == 'POT_wall':
                            self.POT_wall = inp_value
                        if inp_param == 'POT_ZMAX':
                            self.POT_zmax = float(inp_value)
                        if inp_param == 'POT_ZMIN':
                            self.POT_zmin = float(inp_value)
                        if inp_param == 'POT_KZ':
                            self.POT_kz = float(inp_value)
                        if inp_param == 'POT_ATOMSFILE':
                            self.POT_atoms = inp_value
                        if inp_param == 'POT_ATOMSCOL':
                            if inp_value == 'O':
                                self.POT_colA = 55
                                self.POT_colB = 61
                            elif inp_value == 'B':
                                self.POT_colA = 61
                                self.POT_colB = 67
                        if inp_param == 'POT_ATOMSCOLVALUE':
                            self.POT_tag = float(inp_value)

                        # Osmotic Pressure
                        if inp_param == 'OSMO_RUN':
                            self.OP_run = inp_value
                        if inp_param == 'OSMO_ENSEMBLES':
                            self.OP_ensembles = int(inp_value)
                        if inp_param == 'OSMO_EQUIL':
                            self.OP_equil = int(inp_value)
                        if inp_param == 'OSMO_PRODUTION':
                            self.OP_production = int(inp_value)
                        if inp_param == 'OSMO_FREQ':
                            self.OP_freq = int(inp_value)
                        if inp_param == 'OSMO_OUTPUT':
                            self.OP_output = inp_value

        return self


def read_config(configFile):
    return _OpenMMOSMOConfig().readconfig(configFile)


def read_psf(filename):
    return CharmmPsfFile(filename)


def read_pdb(filename):
    return PDBFile(filename)


def read_params(filename):
    charmmExt = ['rtf', 'prm', 'str']
    paramFiles = ()

    for line in open(filename, 'r'):
        if line.find("#") >= 0:
            line = line.split("#")[0]
        parfile = line.strip()
        if len(parfile) > 0:
            ext = parfile.lower().split('.')[-1]
            if ext in charmmExt:
                paramFiles += (parfile,)

    return CharmmParameterSet(*paramFiles)


def load_structure(inp):
    # Load Structure
    psf = read_psf(inp.structure)
    pdb = read_pdb(inp.coordinates)
    params = read_params(inp.param_file)

    # Determine z dimensions
    coords = pdb.positions
    min_crds = coords[0][2] / nanometer
    max_crds = coords[0][2] / nanometer
    for coord in coords:
        min_crds = min(min_crds, coord[2] / nanometer)
        max_crds = max(max_crds, coord[2] / nanometer)
        boxlz = max_crds - min_crds

    # Set PBC box
    psf.setBox(inp.Box_xlen * angstrom,
               inp.Box_ylen * angstrom,
               boxlz*nanometer)

    return psf, pdb, params

