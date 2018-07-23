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

The `effmass` package is aimed towards theoretical solid state physicists and chemists who have a basic familiarity with Python. Depending on the functionality and level of approximation you are looking for, 
it may be that one of the packages listed below will suit your needs better.

# Related packages

Effective mass calculations are implemented in a number of other packages:

[vasppy](https://github.com/bjmorgan/vasppy/) [@Morgan2018]: This is installed as a dependancy of `effmass`. Calculates the effective mass using a least-squares quadratic fit for parabolic dispersions. 

[sumo](https://github.com/SMTG-UCL/sumo): Calculates the effective mass using a least-squares fit for parabolic and non-parabolic dispersions. 

[emc](https://github.com/afonari/emc): Calculates the effective mass *tensor* using a finite-difference method for parabolic dispersions.

[pymatgen](http://pymatgen.org/) [@Ong2013]: This is installed as a dependancy of `effmass`. Calculates an average effective mass *tensor* for non-parabolic dispersions with multiple bands and extrema. Also calculates the Seebeck effective mass as defined [here](https://perso.uclouvain.be/geoffroy.hautier/wp-content/papercite-data/pdf/gibbs2017.pdf).

If you have an update to the information above then please use the Github [issue tracker](https://github.com/lucydot/effmass/issues/). 

## Which features are unique to the effmass package? 

To our knowledge, the following features are unique to this package:

- easily compare the values of curvature effective mass calculated using multiple numerical techniques (least-squares and polynomial fitting)
- tailor the polynomial fitting used to approximate the DFT calculated dispersion: by choosing the order of the polynomial and the energy range to fit over.
- visualise the dispersions used to approximate the DFT calculated dispersion
- quantify non-parabolicity through the Kane dispersion parameters: effective mass at band-edge and alpha
- calculate the optical effective mass assuming a Kane dispersion.

# Acknowledgements

LW would like to thank [Aron Walsh](https://github.com/aronwalsh), [Benjamin Morgan](https://github.com/bjmorgan) and [Jarvist Moore Frost](https://github.com/jarvist) for their guidance during this project.
This package was written during a PhD funded by the EPSRC through the Centre for Doctoral Training in New and Sustainable Photovoltaics (grant no. EP/L01551X/1).
The input data used for developing and testing this package was generated using the [ARCHER UK National Supercomputing Service](http://www.archer.ac.uk). We have access to Archer via our membership of the UK's HEC Materials Chemistry Consortium, which is funded by EPSRC (EP/L000202).

# References






