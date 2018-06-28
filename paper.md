---
title: 'effmass: An effective mass package'
tags:
  - Python
  - materials science
  - DFT
authors:
  - name: Lucy D. Whalley
    orcid: 0000-0002-2992-9871
    affiliation: 1
affiliations:
 - name: Department of Materials, Imperial College London, London, United Kingdom
   index: 1
date: 25 June 2018
bibliography: paper.bib
colorlinks: true
---

# Summary

Many semiconductor properties depend on the response of electrons to an external pertubation.
This perturbation could take the form of an electric field, change in temperature or an applied lattice stress. 
In a crystal, this response depends on the interaction of the electrons with a periodic potential. 
The effective mass approximation assumes that the response of an electron in a periodic potential is equivalent to that of a free electron with a renormalised mass (called the ``effective mass'').
This makes the effective mass a critical parameter in models for the optical and transport properties of a semiconductor.

The effective mass has a number of definitions, depending on the perturbation under consideration. 
The conventional definition of effective mass is inversely proportional to the second derivative of electron energy with respect to electron momentum [@Ashcroft1976, p. 227]. 
This allows the effective mass to be easily calculated from *ab-initio* band structures, and there are existing codes which have implemented this [@Fornari2012 ; @Morgan2018].

We must approximate the band structure with a parabola for the previous definition to be valid [@Ariel2012]. 
However, this approximation breaks down when there is a high concentration of electrons in the material - when, for example, the material is doped or excited under a laser. 
Instead, we can then approximate the band structure with the Kane quasi-linear dispersion [@Kane1957], and the definition of effective mass is adapted accordingly.

[``effmass``](https://github.com/lucydot/effmass) [@Whalley2018] is a Python 3 package for calculating various definitions of effective mass  from the electronic bandstructure of a semiconducting material. 
It contains a core class that calculates the effective mass and other associated properties of selected band structure segments.
[``effmass``](https://github.com/lucydot/effmass) also contains functions for locating band structure extrema, calculating the Kane quasi-linear dispersion parameters and plotting approximations to the true dispersion.
Parsing of electronic structure data is faciliated by the [``vasppy``](https://github.com/bjmorgan/vasppy) [@Morgan2018] package.

# Acknowledgements

LW would like to thank [Aron Walsh](https://github.com/aronwalsh), [Benjamin Morgan](https://github.com/bjmorgan) and [Jarvist Moore Frost](https://github.com/jarvist) for their guidance during this project.
This package was written during a PhD funded by the EPSRC through the Centre for Doctoral Training in New and Sustainable Photovoltaics (grant no. EP/L01551X/1).
The input data used for developing and testing this package was generated using the [ARCHER UK National Supercomputing Service](http://www.archer.ac.uk). We have access to Archer via our membership of the UK's HEC Materials Chemistry Consortium, which is funded by EPSRC (EP/L000202).

# References






