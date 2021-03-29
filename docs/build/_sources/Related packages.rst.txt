Related packages
================

The `effmass` package is aimed towards theoretical solid state physicists and chemists who have a basic familiarity with Python. 
Depending on the functionality and level of approximation you are looking for, 
it may be that one of the packages listed below will suit your needs better.

`vasppy <https://github.com/bjmorgan/vasppy>`_: This is installed as a dependancy of `effmass`. Calculates the effective mass using a least-squares quadratic fit for parabolic dispersions. 

`sumo <https://github.com/SMTG-UCL/sumo>`_: Calculates the effective mass using a least-squares fit for parabolic and non-parabolic dispersions. 

`emc <https://github.com/afonari/emc>`_: Calculates the effective mass *tensor* using a finite-difference method for parabolic dispersions.

`pymatgen <http://pymatgen.org>`_: This is installed as a dependancy of `effmass`. Calculates an average effective mass *tensor* for non-parabolic dispersions with multiple bands and extrema. Also calculates the Seebeck effective mass as defined `here <https://perso.uclouvain.be/geoffroy.hautier/wp-content/papercite-data/pdf/gibbs2017.pdf>`_.

If you have an update to the information above then please use the Github `issue tracker <https://github.com/lucydot/effmass/issues/>`_. 

===================================================
Which features are unique to the `effmass` package?
===================================================

To our knowledge, the following features are unique to `effmass`:

- easily compare the values of curvature effective mass calculated using multiple numerical techniques (least-squares and polynomial fitting)
- tailor the polynomial fitting used to approximate the DFT calculated dispersion: by choosing the order of the polynomial and the energy range to fit over.
- visualise the dispersions used to approximate the DFT calculated dispersion
- quantify non-parabolicity through the Kane dispersion parameters: effective mass at band-edge and alpha
- calculate the optical effective mass assuming a Kane dispersion.
