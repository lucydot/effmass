#! /usr/bin/env python3

from importlib.metadata import version, PackageNotFoundError

angstrom_to_bohr = 1.88973
ev_to_hartree = 0.0367493
kt_to_ev = 0.02588716 # kT/q at T=300K 
AA = 1E-10
q= 1.60217662E-19
boltzmann= 1.38064852E-23

__all__ = ['inputs','outputs','extrema','analysis','dos']

try:
    __version__ = version("effmass")
except PackageNotFoundError:
    # package is not installed
    __version__ = None
