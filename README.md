# Workflow of Running Multiple Walkers Metadynamics Simulations with ASE and PLUMED Using a MACE Potential

The following is the workflow to run multiple walkers metadynamics simulation with MACE as the interatomic potential and ASE as the MD engine:

1. Train the MACE Potential and Prepare the System:
   - Follow the external tutorials to train your MACE potential
   - Equilibrate the system at its local minimum structure from the training set
  
2. Setup the Multiple Walkers Metadynamics Simulations:
   - Decide the number of walkers and 
