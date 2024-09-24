# Analyze the Results and Post-Processing

Besides the diffusion of the CVs, you can perform a block analysis to assess the error and check the simulation's convergence. The PLUMED [error analysis](https://www.plumed-tutorials.org/lessons/21/002/data/NAVIGATION.html) and [metadynamics](https://www.plumed-tutorials.org/lessons/21/004/data/NAVIGATION.html) tutorials discuss these in more detail. This section provides some guidelines on reweighing the simulation to obtain the free energy surface.

## Reweighing Multiple Walkers Metadynamics Simulations

You may have already been introduced to the two main routes of reweighing a metadynamics simulation from the PLUMED metadynamics tutorial: [Calculate the weights from the time-dependence of the metadynamics bias potential](https://pubs.acs.org/doi/10.1021/jp504920s) or [calculate the weights from the metadynamics bias potential obtained at the end of the simulation](https://pubs.acs.org/doi/full/10.1021/ct3002464). This tutorial uses the former way as an example.

