# Analyze the Results and Post-Processing

Besides the diffusion of the CVs, you can perform a block analysis to assess the error and check the simulation's convergence. The PLUMED [error analysis](https://www.plumed-tutorials.org/lessons/21/002/data/NAVIGATION.html) and [metadynamics](https://www.plumed-tutorials.org/lessons/21/004/data/NAVIGATION.html) tutorials discuss these in more detail. This section provides some guidelines on reweighing the simulation to obtain the free energy surface.

## Reweighing Multiple Walkers Metadynamics Simulations

You may have already been introduced to the two main routes of reweighing a metadynamics simulation from the PLUMED metadynamics tutorial: [Calculate the weights from the time-dependence of the metadynamics bias potential](https://pubs.acs.org/doi/10.1021/jp504920s) or [calculate the weights from the metadynamics bias potential obtained at the end of the simulation](https://pubs.acs.org/doi/full/10.1021/ct3002464). This tutorial uses the former way as an example.

To reweigh the frames from multiple walkers metadynamics simulations, you can perform the reweighing for each walker independently and then combine all the histograms to convert the free energy surface. You can construct a `plumed_reweight.dat` file for each walker as follows:
```plumed
# Read COLVAR file
cv1: READ FILE=COLVAR IGNORE_TIME VALUES=cv1
cv2: READ FILE=COLVAR IGNORE_TIME VALUES=cv2
metad: READ FILE=COLVAR IGNORE_TIME VALUES=metad.rbias

# Define weights
weights: REWEIGHT_METAD TEMP=300

# Calculate histograms
hh: HISTOGRAM ARG=cv1,cv2 GRID_MIN=0.1,0.3 GRID_MAX=0.4,5.0 GRID_BIN=255,255 BANDWIDTH=0.01,0.15 LOGWEIGHTS=weights NORMALIZATION=true CLEAR=__FILL__
DUMPGRID GRID=hh FILE=histogram-walker-0.dat STRIDE=__FILL__
```
Note that you will need to determine the values of CLEAR in HISTOGRAM and STRIDE in DUMPGRID from the converged block size of your block analysis. After constructing this `plumed_reweight.dat` file, you can perform the reweighing using PLUMED tools:
```
plumed driver --plumed plumed_reweight.dat --noatoms
```
