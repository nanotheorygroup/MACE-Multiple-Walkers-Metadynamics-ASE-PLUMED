# Input Files

## ASE Input File

This tutorial uses ASE as the MD engine. Before you install ASE, you will need to ensure that your PLUMED is installed with Python wrappers (https://www.plumed.org/doc-v2.9/user-doc/html/_installation.html#installingpython).

To install the development version of ASE, which includes all the latest features, you can run the following commands:
```
git clone https://gitlab.com/ase/ase.git
cd ase
pip install -e .
```

Once ASE is installed, you can run MD simulations with PLUMED with the following example Python script:
```
from ase import Atoms, units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.bussi import Bussi
from ase.md import MDLogger
from ase.io import read, write
from ase.io.trajectory import Trajectory
import torch
from mace.calculators import MACECalculator
from ase.calculators.plumed import Plumed

# Load the trained MACE model
model_path = "./path/to/your/trained/MACE/model.model"
model = MACECalculator(model_path, device='cuda')

# Read the initial atomic configuration from an XYZ file
xyz_file_path = './path/to/your/initial/configuration.xyz'
atoms = read(xyz_file_path)

# Set the cell (volume) and periodic boundary conditions (PBC)
cell = [[12.98, 0.0, 0.0],  # x-axis vector
        [0.0, 12.98, 0.0],  # y-axis vector
        [0.0, 0.0, 12.98]]  # z-axis vector
atoms.set_cell(cell)
atoms.set_pbc([True, True, True])  # Set PBC along x, y, z

# Set up the PLUMED actions
setup = [f"PLUMED INPUT SECTION"]

# Set up the CSVR thermostat for the NVT ensemble
temperature_K = 300
timestep = 0.5  # fs
taut = 1000 # fs
atoms.calc = Plumed(calc=model, input=setup, timestep=timestep * units.fs, atoms=atoms, kT=temperature_K * units.kB, log='PLUMED.OUT')
MaxwellBoltzmannDistribution(atoms, temperature_K=300)
dyn = Bussi(atoms, timestep * units.fs, temperature_K, taut * units.fs)

# Save the trajectory
traj = Trajectory('mace_metad.traj', 'w', atoms)
dyn.attach(traj.write, interval=10)

# Set up a logger to monitor the simulation
dyn.attach(MDLogger(dyn, atoms, 'md_mace_metad.log', header=True, stress=False), interval=10)

# Run the MD simulation
nsteps = 1000000  # Number of MD steps
dyn.run(nsteps)

# Export the last frame in extxyz format for restarting the simulation
write("mace_metad_restart.xyz", atoms, format='extxyz')
```

## PLUMED Input Section

As shown above, the PLUMED actions can be inputted explicitly in the Python script for the ASE input:
```plumed
# Solution file of the first walker = solution/walker-0/MACE_MD_MetaD.py
 
# Default units in PLUMED
UNITS LENGTH=nm TIME=ps ENERGY=kj/mol

# Define the CVs and PLUMED actions
hw: GROUP ATOMS=1,2,5,6,7,9,11,13,15,16,18,19,21,23,25,26,27,29,30,31,33,34,35,36,39,40,43,45,46,47,48,51,52,54,55,57,58,59,61,62,63,65,66,67,68,71,72,73,75,76,78,79,80,81,82,84,85,89,90,91,92,95,96,98,99,100,102,105,106,107,109,110,112,113,115,116,117,119,120,121,122,123,124,126,129,130,132,134,135,136,137,138,140,142,143,144,147,149,151,152,153,155,156,157,158,159,160,161,162,166,167,170,171,173,175,176,177,179,181,182,183,184,186,187,188,189
ow: GROUP ATOMS=1-189 REMOVE=hw
o_no2: GROUP ATOMS=190,191
cv1: DISTANCE ATOMS=192,193
cv2: HBOND_MATRIX ACCEPTORS=o_no2 HYDROGENS=hw DONORS=ow SWITCH={RATIONAL R_0=0.345 NN=8 MM=16} HSWITCH={RATIONAL R_0=0.115 NN=12 MM=24} ASWITCH={RATIONAL R_0=0.167pi NN=6 MM=12} SUM
no1: DISTANCE ATOMS=192,190
no2: DISTANCE ATOMS=192,191

# Apply walls to restrict the exploration of the CVs if necessary
uwallcv1: UPPER_WALLS ARG=cv1 AT=0.40 KAPPA=10000.0
uwall1: UPPER_WALLS ARG=no1 AT=0.15 KAPPA=10000.0
uwall2: UPPER_WALLS ARG=no2 AT=0.15 KAPPA=10000.0
uwallcv2: UPPER_WALLS ARG=cv2.sum AT=5.0 KAPPA=10000.0

# Setup of the multiple walkers metadynamics simulations
metad: METAD ARG=cv1,cv2.sum SIGMA=0.01,0.05 HEIGHT=5.0 PACE=100 TEMP=300.0 BIASFACTOR=15 GRID_MIN=0.00,0.00 GRID_MAX=0.5,6.0 CALC_RCT RCT_USTRIDE=10 WALKERS_N=__FILL__ WALKERS_ID=__FILL__ WALKERS_DIR=__FILL__ WALKERS_RSTRIDE=__FILL__

# Print out the CVs and biases on the fly
PRINT ARG=cv1,uwallcv1.bias,cv2.sum,no1,uwall1.bias,no2,uwall2.bias,metad.bias,metad.rbias STRIDE=10 FILE=COLVAR
FLUSH STRIDE=100
```
To set up a metadynamics simulation with multiple walkers, you need to specify the number of walkers by WALKERS_N, the walker index by WALKERS_ID (which starts from 0), the directory which stores the bias potentials by WALKERS_DIR, and the frequency of the walkers reading the bias potential files by WALKERS_RSTRIDE.
