# Model System - The Dissociation of Nitrate Anion at its Triplet State in an Aqueous Solution

Nitrate anion ($\mathrm{NO_3^-}$) in aerosol particles is an essential sink of nitrogen oxide species (NOx). Its photodissociation is a 'renoxification' process, 
which converts nitrate anion solvated in water or deposited on surfaces back into NOx to the atmosphere. The dissociation of nitrate anion at its triplet state can follow two channels:

```math
\mathrm{NO_3^-} \overset{\textit{h}ν}{\rightarrow} \mathrm{NO_2} + \mathrm{O^{\cdot -}}
```
```math
\mathrm{NO_3^-} \overset{\textit{h}ν}{\rightarrow} \mathrm{NO_2^-} + \mathrm{O(^{3}P)}
```

Despite the well-studied macroscopic kinetics of the two channels, the microscopic details and their connection to the thermodynamics and kinetics of the dissociation are still inconclusive.
Previous experiments have shown that nitrate anion photodissociation in aqueous solutions has a low quantum yield of ~ 1%.
We previously employed _ab initio_ molecular dynamics simulations at the level of density functional theory with multiple walkers metadynamics to explore the two channels in an aqueous solution to unravel the atomistic and electronic structure details.
We located a solvation cage complex that can be identified as a meta-stable state that requires additional thermal energy to complete the dissociation of the N-O bond at the triplet state.
This meta-stable state allows the photo-fragments to recombine or deactivate through non-radiative processes.
In this tutorial, a machine learning interatomic potential was trained based on the first-principles data to demonstrate how to perform multiple walkers metadynamics.

Two CVs can be chosen for the metadynamics simulation (see the figure below).
