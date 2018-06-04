# effmass

[![Build Status](https://travis-ci.org/lucydot/effmass.svg?branch=master)](https://travis-ci.org/lucydot/effmass)
[![Test Coverage](https://codeclimate.com/github/lucydot/effmass/badges/coverage.svg)](https://codeclimate.com/github/lucydot/effmass/coverage)

Python package for calculating various definitions of effective mass from the electronic bandstructure of a semiconducting material. 

The package can:

**1. Read in a bandstructure**
This requires the `VASP` output files `PROCAR` and `OUTCAR`. It is assumed you have walked through a 1D slice of the Brillouin Zone, capturing the maxima and minima of interest.

**2. Locate the extrema**
These correspond to the valence band maxima and conduction band minima.

**3. Fit a polynomial**
For each extremum, a polynomial over a specified energy range is fitted.

**4. Calculate effective mass**
The curvature (aka inertial), transport, and optical effective mass are calculated using the derivatives of each polynomial function. 

**5. Assess the extent of non-parabolicity**
The transport effective mass is used to calculate the parameters of the Kane quasi-linear dispersion.

Dependancies:

- [vasppy](https://github.com/lucydot/vasppy) 
- [numpy](http://www.numpy.org/)
- [matplotlib](https://matplotlib.org/)
