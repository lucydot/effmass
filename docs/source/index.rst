.. effmass documentation master file, created by
   sphinx-quickstart on Thu Mar 22 12:48:42 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Effmass
=======

Python package for calculating the effective mass from a bandstructure calculations. Rahter than fitting a quadratic and relying upont he parabolic approximation, we fit a higher order polynomial and use this calculate a non-parabolic optical effective mass. This documentation is best understood by reading this paper which makes use of the package.

.. figure:: .static/electronworries.png
    :align: right
    :figwidth: 300px
    :alt: "Do I look fat in this potential?" -- electron

The package can:

1. **Read in a bandstructure**: Currently this requires the VASP output files PROCAR and OUTCAR. It is assumed you have walked through a 1D slice of the Brillouin Zone, capturing the maxima and minima of interest.

2. **Locate the extrema**: These correspond to the valence band maxima and conduction band minima.

4. **Quantify non-parabolicity**: By calculating the alpha parmeter as defined in the Kane-quasi linear dispersion.

s. **Calculate bandedge (parabolic) effective mass and optical (non-parabolic) effective mass**: The optical effective mass can weighted according to the density of states (with option to read in from DOSCAR) and probability of occupation (fermi-dirac distribution).

Effmass is free and open source. If you have any issues, please raise them in the `Github issues tracker <https://github.com/lucydot/effmass/issues>`_.
If you use this package please cite.

Details
-------

:Authors: Lucy Whalley
:Contact: lucywhalley@gmail.com
:GitHub: https://github.com/lucydot
:License: MIT
:Citation: :download:`bibtex<.static/CITATION>`
:Build: Travis CI badge here

Documentation
-------------

.. toctree::
   :maxdepth: 1
   
   Installation
   Usage
   API docs




