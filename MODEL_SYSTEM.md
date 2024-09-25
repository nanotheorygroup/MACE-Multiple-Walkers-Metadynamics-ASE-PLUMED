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

![The schematic representation of the CVs used in the current study. The two CVs, _d_ are _n_ are colored in cyan and orange respectively](/CVs_Picture.png)

Two CVs are chosen for the current metadynamics simulations. The first CV, _d_, is the distance between atom N and O1, corresponding to the dissociation coordinate. However, it is not enough to differentiate the two dissociation channels. Since the two channels only differ in the net negative charge localization and are comparable in an aqueous solution, one can hypothesize that the solvation structure of $\mathrm{NO_3^-}$ is crucial to determine which channel it will eventually proceed. Based on this hypothesis, a second CV, _n_, is used to describe the number of H atoms surrounding atoms O2 and O3. This number, _n_, is defined as follow:

```math
n = \sum_{i\in\mathrm{O2, O3}}\sum_{j\in\mathrm{H_{w}}}s_{ij}
```
```math
s_{ij} = \frac{1}{2}\left[\frac{1-(\frac{r_{ij}}{r_{0}})^n}{1-(\frac{r_{ij}}{r_{0}})^m}\right]
```
where $r_{ij}$ is the distance between atoms _i_ and _j_. $r_{0}$ is set to be 2.55 Å, and n and m are set to 8 and 16, respectively.
