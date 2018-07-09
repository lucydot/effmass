========
Features
========

`effmass` can:

**Read in a bandstructure:**
This requires the ``VASP`` output files ``PROCAR`` and ``OUTCAR``. It is assumed you have walked through a 1D slice of the Brillouin Zone, capturing the maxima and minima of interest. `effmass` uses the Python package `vasppy <https://github.com/bjmorgan/vasppy>`_ for parsing ``VASP`` output. 

**Locate extrema:**
These correspond to the valence band maxima and conduction band minima. Maxima and minima within a certain energy range can also be located.

**Calculate curvature, transport and optical effective masses:**
The curvature (aka inertial) and transport effective masses are calculated using the derivatives of a fitted polynomial function. The optical effective mass can also be calculated assuming a Kane dispersion.

**Assess the extent of non-parabolicity**
Parameters of the Kane quasi-linear dispersion are calculated to quantify the extent of non-parabolicity over a given energy range.

**Calculate the quasi-fermi level for a given carrier concentration**
This requires the ``VASP`` output file ``DOSCAR``. Using density-of-states data and assuming no thermal smearing, `effmass` can calculate the energy to which states are occupied. This is a useful proxy to the quasi-Fermi level.

**Plot fits to the dispersion**
Selected bandstructure segments and approximations to the dispersion (assuming a Kane, quadratic, or higher order fit) can be visualised.


.. figure:: .static/screenshot.png
    :figwidth: 400px
    :alt: "effmass screenshot"

*Which features are unique to the effmass package?*

Effective mass calculations are implemented in a number of other packages, including  `vasppy <https://github.com/bjmorgan/vasppy>`_ and `pymatgen <http://pymatgen.org/>`_, both of which are dependancies of `effmass`. However, to our knowledge, the following is unique to our implementation:

- easily compare the values of curvature effective mass calculated using multiple numerical techniques (least-squares and polynomial fitting)
- tailor the polynomial fitting used to approximate the DFT calculated dispersion: by choosing the order of the polynomial and the energy range to fit over.
- visualise the dispersions used to approximate the DFT calculated dispersion
- quantify non-parabolicity through the Kane dispersion parameters: effective mass at band-edge and alpha
- calculate the optical effective mass assuming a Kane dispersion
