# effmass

[![PyPI version](https://badge.fury.io/py/effmass.svg)](https://badge.fury.io/py/effmass)
[![Documentation Status](https://readthedocs.org/projects/effmass/badge/?version=latest)](https://effmass.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.com/lucydot/effmass.svg?branch=master)](https://travis-ci.com/lucydot/effmass)
[![Test Coverage](https://codeclimate.com/github/lucydot/effmass/badges/coverage.svg)](https://codeclimate.com/github/lucydot/effmass/coverage)
[![DOI](https://zenodo.org/badge/136037407.svg)](https://zenodo.org/badge/latestdoi/136037407)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![JOSS status](http://joss.theoj.org/papers/389754561f0710b756514b8cb9ac0e6a/status.svg)](http://joss.theoj.org/papers/389754561f0710b756514b8cb9ac0e6a)

💃 Effmass now has a command line interface  
💃 Effmass now supports FHI-Aims, Castep and ASE  
⚠️ The `Data` class has now been renamed `DataVasp`. You may need to update your scripts!

`effmass` is a [peer-reviewed](https://joss.theoj.org/papers/10.21105/joss.00797) Python (3.6+) package for calculating various definitions of effective mass from the electronic bandstructure of a semiconducting material. It consists of a core `Segment` class that calculates the effective mass and other associated properties of selected bandstructure segments. The programme also contains functions for locating bandstructure extrema, constructing segments and plotting approximations to the dispersion. There is a command line interface for calculating parabolic effective mass values, and an API for the more complex, non-parabolic definitions of effective mass.

If you use `effmass` for your published research [please cite accordingly](#citing-effmass).

## What's new?

💃 `effmass` now interfaces with more codes:
- it can read in Castep output data (in addition to Vasp and FHI-aims)
- it can work with ASE bandstructure objects

💃 `effmass` now includes a command line interface

![](./cli2.gif)  

As a result of these changes, and with view to supporting more DFT codes in the future, the `Data` class has been renamed to `DataVasp` ⚠️ **On updating to the latest version of effmass you may need to update your scripts / Jupyter Notebook to reflect this change** ⚠️

## Features

`effmass` can:

**Read in a bandstructure:**
It is assumed you have used a DFT calculator to walk through a 1D slice of the Brillouin Zone, capturing the maxima and minima of interest. `effmass` uses the Python package [vasppy](https://github.com/bjmorgan/vasppy) for parsing `VASP` output.

**Locate extrema:**
These correspond to the valence band maxima and conduction band minima. Maxima and minima within a certain energy range can also be located.

**Calculate curvature, transport and optical effective masses:**
The curvature (aka inertial) and transport masses are calculated using the derivatives of a fitted polynomial function. The optical effective mass can also be calculated assuming a Kane dispersion.

**Assess the extent of non-parabolicity:**
Parameters of the Kane quasi-linear dispersion are calculated to quantify the extent of non-parabolicity over a given energy range. 

**Calculate the quasi-fermi level for a given carrier concentration:**
Using density-of-states data and assuming no thermal smearing, `effmass` can calculate the energy to which states are occupied. This is a useful approximation to the quasi-Fermi level. *Note: this is only supported for VASP and requires the output file `DOSCAR`.* 

**Plot fits to the dispersion:**
Selected bandstructure segments and approximations to the dispersion (assuming a Kane, quadratic, or higher order fit) can be visualised.

Depending on the functionality and level of approximation you are looking for, 
it may be that one of the packages listed [here](https://effmass.readthedocs.io/en/latest/Related%20packages.html) will suit your needs better.

## Supported Codes

`effmass` currently supports `VASP`, `FHI-Aims`, `Castep` and `ASE`. In the near future we hope to play nicely with other codes that interface with the ASE bandstructure class, and pymatgen. We especially welcome contributions that will help make `effmass` available to more researchers.

## Installation

`effmass` can then be installed using the Python package manager `pip`:

```
pip install effmass
```

If you use conda/anaconda the safest thing to do is to create a new environment and then install effmass:

```
conda create --name effmass
conda activate effmass
conda install pip
pip install effmass
```

If you do not use `pip` you can download and install the latest release from [GitHub](https://github.com/lucydot/effmass/releases) and install:
```
cd effmass
python setup.py install
```

## Command Line Interface

The [command line interface](#Installation) provides basic functionality for calculating parabolic effective masses.
For those who have a basic familiarity with Python there is an API which provides access to all features, including non-parabolic effective mass definitions. 

To start the command line interface simply type

```
effmass
```

and follow the prompts. You are asked if you would like to print a plot of the segments found - we recommend that you do this, to check that the segments are "sensible". You are also asked if you would like to print a summary file - again, we recommend that you do this, so that you have a record of the CLI options chosen.

## Documentation

- An overview of the features of effmass, along with example code for Vasp and FHI-aims output data, is contained in a Jupyter notebook [here](https://nbviewer.jupyter.org/github/lucydot/effmass/blob/master/Tutorial.ipynb).
- Additional Jupyter notebook examples for the Castep and ASE interfaces are [here](https://nbviewer.jupyter.org/github/lucydot/effmass/blob/master/tests/Castep_ASE_interface.ipynb).
- The API documentation is [here](https://effmass.readthedocs.io/en/latest/).
- Further details about the various effective mass definitions implemented in `effmass` can be found in Phys. Rev. B **99** (8), 085207, which is also [available on arXiv](https://arxiv.org/pdf/1811.02281.pdf).
- The source code is available as a git repository at [https://github.com/lucydot/effmass](https://github.com/lucydot/effmass).

### Running notebook examples

If you want to run the jupyter notebook examples/tutorials you will also need to install `notebook`:

```
pip install notebook
```

To run the notebook, run the following command at the Terminal (Mac/Linux) or Command Prompt (Windows):

``` 
jupyter notebook
```

This will open a web browser tab, which you can use to navigate to the notebook examples.

## Publications using `effmass`

A [number of publications](https://scholar.google.co.uk/scholar?oi=bibs&hl=en&cites=12032412581356217625) have used `effmass`.

`effmass` was initially developed for a project that has been published as *Impact of nonparabolic electronic band structure on the optical and transport properties of photovoltaic materials*  Phys. Rev. B **99** (8), 085207. This paper is also [avaiable on arXiv](https://arxiv.org/pdf/1811.02281.pdf). The [paper directory](https://github.com/lucydot/effmass/paper) contains the Vasp input data (POSCAR), Vasp output data (OUTCAR/PROCAR) and band structures generated for this study.

## Questions, bug reports, feature requests

Please use the Github [issue tracker](https://github.com/lucydot/effmass/issues/) for any questions, feature requests or bug reports. Please do not contact the developers via email unless there is a specific reason you do not want the conversation to be public.

## Development

If you would like to contribute please do so via a pull request. All contributors must read and respect the [code of conduct](https://github.com/lucydot/effmass/blob/master/CODE_OF_CONDUCT.md). In particular, we welcome contributions which would extend `effmass` so that it is able to parse output from other electronic structure codes. 

## Tests

Automated testing of the latest commit happens [here](https://travis-ci.com/lucydot/effmass).

You can also run tests locally:
```
pip install effmass[tests]
cd effmass
python -m pytest
```

## Citing `effmass`

If you use this code in your research, please cite the following paper:

Whalley, Lucy D. (2018). *effmass - an effective mass package*. The Journal of Open Source Software, 3(28) 797.
Link to paper [here](https://joss.theoj.org/papers/10.21105/joss.00797).

### Bibtex

```
@misc{Whalley_JOSS2018,
  author       = {Lucy D. Whalley},
  title        = {effmass: An effective mass package},
  volume       = {3},
  issue        = {28},
  pages        = {797},
  month        = {Aug},
  year         = {2018},
  doi          = {10.21105/joss.00797},
  url          = {http://joss.theoj.org/papers/10.21105/joss.00797}
}
```

## Contributors

Lead developer: 
[Lucy Whalley](https://lucydot.github.io), a.k.a [lucydot](https://github.com/lucydot)

Contributors: 
Matthias Goloumb (Support for FHI-Aims), a.k.a [MatthiasGolomb](https://github.com/MatthiasGolomb) //
Sean Kavanagh (Documentation), a.k.a [kavanase](https://github.com/kavanase) //
Benjamin Morgan (Vasppy compatability), a.k.a [bjmorgan](bjmorgan) //


