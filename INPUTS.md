# Input Files

## ASE Input File

This tutorial uses ASE as the MD engine. To install the development version of ASE, which includes all the latest features, you can run the following commands:
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
model_path = "path/to/your/trained/MACE/model.model"
model = MACECalculator(model_path, device='cuda')

# Read the initial atomic configuration from an XYZ file
xyz_file_path = 'path/to/your/initial/configuration.xyz'
atoms = read(xyz_file_path)

# Set the cell (volume) and periodic boundary conditions (PBC)
cell = [[xx.xx, 0.0, 0.0],  # x-axis vector
        [0.0, yy.yy, 0.0],  # y-axis vector
        [0.0, 0.0, zz.zz]]  # z-axis vector
atoms.set_cell(cell)
atoms.set_pbc([True, True, True])  # Set PBC along x, y, z
```
