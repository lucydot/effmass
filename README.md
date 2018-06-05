# effmass

[![Build Status](https://travis-ci.com/lucydot/effmass.svg?branch=master)](https://travis-ci.com/lucydot/effmass)
[![Test Coverage](https://codeclimate.com/github/lucydot/effmass/badges/coverage.svg)](https://codeclimate.com/github/lucydot/effmass/coverage)

`effmass` is a Python package for calculating various definitions of effective mass from the electronic bandstructure of a semiconducting material. It consists of a core class that calculates properties of selected bandstructure segments. The module also contains functions for locating bandstructure extrema and plotting fits to the dispersion.

Examples are provided in a Jupyter notebook [here](nbviewer.jupyter.org/github/lucydot/effmass/blob/master/paper/notebook.ipynb).
API documentation is [here](effmass.readthedocs.io/en/latest/).
Source code is available as a git repository at [https://github.com/lucydot/effmass](https://github.com/lucydot/effmass).

## Features

`effmass` can:

**Read in a bandstructure**
This requires the `VASP` output files `PROCAR` and `OUTCAR`. It is assumed you have walked through a 1D slice of the Brillouin Zone, capturing the maxima and minima of interest. `effmass` uses the Python package [vasppy](https://github.com/bjmorgan/vasppy) for parsing `VASP` output.

**Locate the extrema**
These correspond to the valence band maxima and conduction band minima. Maxima and minima within a certain energy range can also be located.

**Calculate curvature, transport and optical effective mass**
The curvature (aka inertial) and transport effective masses are calculated using the derivatives of a fitted polynomial function. The optical effective mass

**Assess the extent of non-parabolicity**
Parameters of the Kane quasi-linear dispersion are calculated to quantify the extent of non-parabolicity over a given energy range. add

**Calculate electron fill level from DOS data**

**Plot fits to the dispersion**

## Installation

TODO.

## Tests

Automated testing of the latest commit happens [here](https://travis-ci.com/lucydot/effmass).

Manual tests can be run using 
```
python -m pytest
```

This code has been tested with Python versions 3.6.

## Documentation

An overview of the features of effmass along with example code is contained in a [Jupyter notebook](nbviewer.jupyter.org/github/lucydot/effmass/blob/master/paper/notebook.ipynb), which is available in the `paper` directory.

API documentation is available [here](effmass.readthedocs.io/en/latest/).
