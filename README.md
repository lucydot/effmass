# effmass

[![PyPI version](https://badge.fury.io/py/effmass.svg)](https://badge.fury.io/py/effmass)
[![Documentation Status](https://readthedocs.org/projects/effmass/badge/?version=latest)](https://effmass.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.com/lucydot/effmass.svg?branch=master)](https://travis-ci.com/lucydot/effmass)
[![Test Coverage](https://codeclimate.com/github/lucydot/effmass/badges/coverage.svg)](https://codeclimate.com/github/lucydot/effmass/coverage)
[![DOI](https://zenodo.org/badge/136037407.svg)](https://zenodo.org/badge/latestdoi/136037407)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![JOSS status](http://joss.theoj.org/papers/389754561f0710b756514b8cb9ac0e6a/status.svg)](http://joss.theoj.org/papers/389754561f0710b756514b8cb9ac0e6a)

`effmass` is a Python 3 package for calculating various definitions of effective mass from the electronic bandstructure of a semiconducting material. It consists of a core class that calculates the effective mass and other associated properties of selected bandstructure segments. The module also contains functions for locating bandstructure extrema and plotting approximations to the dispersion.

Examples are provided in a Jupyter notebook [here](https://nbviewer.jupyter.org/github/lucydot/effmass/blob/master/paper/notebook.ipynb).
API documentation is [here](https://effmass.readthedocs.io/en/latest/).
Source code is available as a git repository at [https://github.com/lucydot/effmass](https://github.com/lucydot/effmass).

## Features

`effmass` can:

**Read in a bandstructure:**
This requires the `VASP` output files `PROCAR` and `OUTCAR`. It is assumed you have walked through a 1D slice of the Brillouin Zone, capturing the maxima and minima of interest. `effmass` uses the Python package [vasppy](https://github.com/bjmorgan/vasppy) for parsing `VASP` output.

**Locate extrema:**
These correspond to the valence band maxima and conduction band minima. Maxima and minima within a certain energy range can also be located.

**Calculate curvature, transport and optical effective masses:**
The curvature (aka inertial) and transport masses are calculated using the derivatives of a fitted polynomial function. The optical effective mass can also be calculated assuming a Kane dispersion.

**Assess the extent of non-parabolicity:**
Parameters of the Kane quasi-linear dispersion are calculated to quantify the extent of non-parabolicity over a given energy range. 

**Calculate the quasi-fermi level for a given carrier concentration:**
This requires the `VASP` output file `DOSCAR`. Using density-of-states data and assuming no thermal smearing, `effmass` can calculate the energy to which states are occupied. This is a useful approximation to the quasi-Fermi level.

**Plot fits to the dispersion:**
Selected bandstructure segments and approximations to the dispersion (assuming a Kane, quadratic, or higher order fit) can be visualised.

## Development

Please use the Github [issue tracker](https://github.com/lucydot/effmass/issues/) for feature requests and bug reports. 

If you would like to contribute please do so via a pull request. All contributors must read and respect the [code of conduct](https://github.com/lucydot/effmass/blob/master/CODE_OF_CONDUCT.md). In particular, we welcome contributions which would extend `effmass` so that it is able to parse output from other electronic structure codes. 

## Installation

```
pip install effmass
```

Or download the latest release from [GitHub](https://github.com/lucydot/effmass/releases), and install
```
cd effmass
python setup.py install
```

Or clone the latest development version
```
git clone git@github.com:lucydot/effmass.git
```
and install the same way.
```
cd effmass
python setup.py install 
```

## Tests

Automated testing of the latest commit happens [here](https://travis-ci.com/lucydot/effmass).

Manual tests can be run using 
```
python -m pytest
```

This code has been tested with Python versions 3.6.

## Documentation

An overview of the features of effmass along with example code is contained in a [Jupyter notebook](https://nbviewer.jupyter.org/github/lucydot/effmass/blob/master/paper/notebook.ipynb), which is available in the `paper` directory.

API documentation is available [here](https://effmass.readthedocs.io/en/latest/).

## Citing `effmass`

This code can be cited as:

Whalley, Lucy D. (2018). *effmass - an effective mass package*. http://doi.org/10.5281/zenodo.1285137

### Bibtex

```
@misc{Whalley_Zenodo2018a,
  author       = {Lucy D. Whalley},
  title        = {effmass: 0.1.0},
  month        = jun,
  year         = 2018,
  doi          = {10.5281/zenodo.1285137},
  url          = {http://dx.doi.org/10.5281/zenodo.1285137}
}
```
