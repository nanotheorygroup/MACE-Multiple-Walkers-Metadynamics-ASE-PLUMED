# Analyze the Results and Post-Processing

Besides the diffusion of the CVs, you can perform a block analysis to assess the error and check the simulation's convergence. The PLUMED [error analysis](https://www.plumed-tutorials.org/lessons/21/002/data/NAVIGATION.html) and [metadynamics](https://www.plumed-tutorials.org/lessons/21/004/data/NAVIGATION.html) tutorials discuss these in more detail. This section provides some guidelines on reweighing the simulation to obtain the free energy surface.

## Reweighing Multiple Walkers Metadynamics Simulations

You may have already been introduced to the two main routes of reweighing a metadynamics simulation from the PLUMED metadynamics tutorial: [Calculate the weights from the time-dependence of the metadynamics bias potential](https://pubs.acs.org/doi/10.1021/jp504920s) or [calculate the weights from the metadynamics bias potential obtained at the end of the simulation](https://pubs.acs.org/doi/full/10.1021/ct3002464). This tutorial uses the former way as an example.

To reweigh the frames from multiple walkers metadynamics simulations, you can perform the reweighing for each walker independently and then combine all the histograms to convert the free energy surface. You can construct a `plumed_reweight.dat` file for each walker as follows:
```plumed
# Read the COLVAR file
cv1: READ FILE=COLVAR IGNORE_TIME VALUES=cv1
cv2: READ FILE=COLVAR IGNORE_TIME VALUES=cv2
metad: READ FILE=COLVAR IGNORE_TIME VALUES=metad.rbias

# Define weights
weights: REWEIGHT_METAD TEMP=300

# Calculate histograms
hh: HISTOGRAM ARG=cv1,cv2 GRID_MIN=0.1,0.3 GRID_MAX=0.4,5.0 GRID_BIN=255,255 BANDWIDTH=0.01,0.15 LOGWEIGHTS=weights NORMALIZATION=true CLEAR=__FILL__
DUMPGRID GRID=hh FILE=histogram.dat STRIDE=__FILL__
```
Note that you will need to determine the values of CLEAR in HISTOGRAM and STRIDE in DUMPGRID for the block size in your block analysis. You can perform the reweighing using PLUMED tools:
```
plumed driver --plumed plumed_reweight.dat --noatoms
```

After reweighing all the walkers, you can combine all the histograms to calculate the free energy and the associated error. You can either do it using PLUMED or with the following Python script:
```
import numpy as np
import matplotlib.pyplot as plt
import glob
import scipy.stats

# This is a script to calculate the free energy surface from the histograms obtained from multiple walkers metadynamics simulations with two collective variables (CV)
# This script uses two walkers as an example

# Function to read in histogram data, normalization, and number of bins for the two CVs
def readhistogram( fname ) :
        # Read in the histogram data
        data = np.loadtxt( fname )
        with open( fname, "r" ) as myfile :
                for line in myfile :
                        if line.startswith("#! SET normalisation") : norm = line.split()[3]
                        if line.startswith("#! SET nbins_cv1") : nbin_cv1 = int(line.split()[3]) + 1
                        if line.startswith("#! SET nbins_cv2") : nbin_cv2 = int(line.split()[3]) + 1
        return float(norm), data, nbin_cv1, nbin_cv2

# Calculate weighted average
norm0, hist0, binx0, biny0 = readhistogram( "./path/to/your/walker-0/directory/histogram.dat" )
N, average = norm0, norm0*hist0[:,2]
norm1, hist1, binx1, biny1 = readhistogram( "./path/to/your/walker-1/directory/histogram.dat" )
N, average = N + norm1, average + norm1*hist1[:,2]
for filen in glob.glob( "./*/analysis.*.histogram.dat" ) :
    norm, histn, binxn, binyn = readhistogram(filen)
    N, average = N + norm, average + norm*histn[:,2]
average = average / N

# Calculate errors
error = norm0*norm0*(hist0[:,2]-average[:])**2
error = error + norm1*norm1*(hist0[:,2]-average[:])**2
for filen in glob.glob( "./*/analysis.*.histogram.dat" ) :
    norm, histn, binxn, binyn = readhistogram(filen)
    error = error + norm*norm*(histn[:,2]-average[:])**2
error = np.sqrt( (error / (N*N)) )

# Convert to free energy
fes = -2.494*np.log( average )
fes -= np.amin(fes)

# Convert to error in fes
ferr = 2.494*error / average
with open('fes-ferr-block-average.dat', 'w') as fout:
    for i in range(len(fes)):
        fout.write("%f  %f  %f  %f\n" % (hist0[:,0][i]*10, hist0[:,1][i], fes[i], ferr[I]))
```

