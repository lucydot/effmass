#! /usr/bin/env python3

"""
A module for storing electronic structure data and user settings.

The module contains a :class:`Data` class which parses OUTCAR and PROCAR files using the `vasppy <https://github.com/bjmorgan/vasppy>`_ package. 
A function for parsing DOSCAR files is also provided. 
A :class:`Settings` class stores analysis parameters set by the user.

"""
from vasppy import procar, outcar
from effmass import extrema
import math
import warnings
import numpy as np


class Settings:
    """Class for setting analysis parameters.

    Attributes:     energy_range (float): energy in kT over which the
    segment extends.     extrema_search_depth (float): energy in kT from
    bandedge over which to search for extrema.     degree_bandfit (int):
    the degree of the polynomial which is used to fit to dispersion data
    when calculating the transport mass.
    """

    def __init__(self,
                 energy_range=0.25,
                 extrema_search_depth=0.025,
                 bandfit=6):
        """Initialises an instance of the Settings class and checks input using
        :meth:`check_settings()`.

        Args:
            energy_range (float): energy in eV over which the segment extends. Defaults to 0.25 eV.
            extrema_search_depth (float): energy in eV from bandedge over which to search for extrema. Defaults to 0.025 eV.
            degree_bandfit (int): the degree of the polynomial which is used to fit to dispersion data when calculating the transport mass.

        Returns:
            None.
        """
        self.energy_range = energy_range
        self.extrema_search_depth = extrema_search_depth
        self.degree_bandfit = bandfit
        self.check_settings()

    def check_settings(self):
        """Check that Settings class attributes are sane.

        Args:
            None.

        Returns:
            None.
        """
        assert (self.energy_range >
                0), "The energy range must be a positive number"
        assert (self.extrema_search_depth >
                0), "The energy depth must be a positive number"
        assert (
            type(self.degree_bandfit) == int and self.degree_bandfit > 1
        ), "The bandfit degree must be a positive integer greater than 1"


class Data():
    r"""
    Class for parsing and storing data from a vasp calculation.

    Attributes:
        spin_channels (int): 1 (non-spin-polarised), 2 (spin-polarised), 4 (spin-orbit coupling).
        number_of_kpoints (int): the number of k-points per band.
        number_of_bands (int): the number of bands.
        number_of_ions (int): the number of ions.
        kpoints (array(float)): 2-dimensional array with shape (number_of_kpoints, 3). Each row contains the fractional coordinates of a kpoint [kx,ky,kz].
        energies (array(float)): 2-dimensional array with shape (number_of_bands,number_of_kpoints). Each row contains energies of eigenstates in eV for a particular band. 
        occupancy (array(float)): 2-dimensional array with shape (number_of_bands,number_of_kpoints). Each row contains occupation number of the eigenstates for a particular band. Values range from 0-1 (spin-polarised) or 0-2 (non-spin-polarised).     
        reciprocal_lattice (list(float)): the reciprocal lattice vectors in format [[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]], units Angstrom :math:`^{-1}`.
        CBM (float): the conduction band minimum energy in eV.
        VBM (float): the valence band maximum in eV.
        fermi_energy (float): the fermi energy in eV. Automatically set to the mean of Data.CBM and Data.VBM. 
        dos (array): 2-dimensional array. Each row contains density of states data (units "number of states / unit cell")  at a given energy: [energy(float),dos(float)]. 
        integrated_dos: 2-dimensional array. Each row contains integrated density of states data at a given energy: [energy(float),integrated_dos(float)].
    """

    def __init__(self, outcar_path, procar_path, ignore=0, **kwargs):
        r"""
        Initialises an instance of the :class:`~effmass.inputs.Data` class and checks data using :meth:`check_data`.

        Args:
            outcar_path (str): The path to the OUTCAR file
            procar_path (:obj:`str` or :obj:`list`): The path(s) to one or more PROCAR files.
            
            ignore (int): The number of kpoints to ignore at the beginning of the bandstructure slice through kspace (useful for hybrid calculations where zero weightings are appended to a previous self-consistent calculation).
            **kwargs: Additional keyword arguments for reading the PROCAR file(s). 
 
        Returns: 
            None.
        """
        assert (type(outcar_path) == str), "The OUTCAR path must be a string"
        assert (type(ignore) == int and ignore >= 0
                ), "The number of kpoints to ignore must be a positive integer"

        reciprocal_lattice = outcar.reciprocal_lattice_from_outcar(outcar_path)
        if isinstance(procar_path, list):
            vasp_data = procar.Procar.from_files(procar_path, **kwargs)
        elif isinstance(procar_path, str):
            vasp_data = procar.Procar.from_file(procar_path, **kwargs)
        else:
            raise TypeError('procar_path must be a string or list of strings')

        self.spin_channels = vasp_data.spin_channels
        self.number_of_bands = vasp_data.number_of_bands
        self.number_of_ions = vasp_data.number_of_ions

        number_of_kpoints = vasp_data.number_of_k_points
        vasp_data_energies = np.array( [ band.energy for band in np.ravel( vasp_data.bands ) ] )
        vasp_data_occupancies = np.array( [ band.occupancy for band in np.ravel( vasp_data.bands ) ] )
        if vasp_data.calculation['spin_polarised']: # to account for the change in PROCAR format for calculations with 2 spin channels (1 k-point block ---> 2 k-point blocks)
            energies = np.zeros([self.number_of_bands*2,number_of_kpoints]) # This is a very ugly way to slice 'n' dice. Should avoid creating new array and use array methods instead. But it does the job so will keep for now.
            for i in range(self.number_of_bands):
                energies[i] = vasp_data_energies.reshape(
                                            number_of_kpoints*2,  # factor of 2 for each kpoint block
                                            self.number_of_bands).T[i][:number_of_kpoints]
                energies[self.number_of_bands+i] = vasp_data_energies.reshape(
                                            number_of_kpoints*2,
                                            self.number_of_bands).T[i][number_of_kpoints:]
            occupancy = np.zeros([self.number_of_bands*2,number_of_kpoints])
            for i in range(self.number_of_bands):
                occupancy[i] = vasp_data_occupancies.reshape(
                                                 number_of_kpoints*2,
                                                 self.number_of_bands).T[i][:number_of_kpoints]
                occupancy[self.number_of_bands+i] = vasp_data_occupancies.reshape(
                                                 number_of_kpoints*2,
                                                 self.number_of_bands).T[i][number_of_kpoints:]
        else:
            energies = vasp_data_energies.reshape(
                                            number_of_kpoints,
                                            self.number_of_bands).T
            occupancy = vasp_data_occupancies.reshape(
                                                 number_of_kpoints,
                                                 self.number_of_bands).T
        
        # remove values which are from the self-consistent calculation prior to the bandstructure calculation (workflow for hybrid functionals)
        self.energies = np.delete(energies,list(range(ignore)),1)
        self.occupancy = np.delete(occupancy,list(range(ignore)),1)
        self.number_of_kpoints = number_of_kpoints - ignore
        
        # handle negative occupancy values
        if np.any(self.occupancy < 0):
            warnings.warn("One or more occupancies in your PROCAR file are negative. All negative occupancies will be set to zero.")
            self.occupancy[ self.occupancy < 0 ] = 0.0

        self.kpoints = np.array( [ kp.frac_coords 
            for kp in vasp_data.k_points[ignore:vasp_data.number_of_k_points] ] )
        self.reciprocal_lattice = reciprocal_lattice * 2 * math.pi
        self.CBM = extrema._calc_CBM(self.occupancy, self.energies)
        self.VBM = extrema._calc_VBM(self.occupancy, self.energies)
        self.fermi_energy = (self.CBM + self.VBM) / 2
        self.dos = []
        self.integrated_dos = []
        self.check_data()

    def check_data(self):
        """Check that data class attributes are sane.

        Args:
            None.

        Returns:
            None.

        Notes:
            There is a `sanity_check` method which runs automatically when reading data in using the `vasppy.procar <http://vasppy.readthedocs.io/en/latest/vasppy.html#module-vasppy.procar>`_ module.
        """
        assert (
            ((self.spin_channels == 1) | (self.spin_channels == 2) |
             (self.spin_channels == 4)) is True
        ), "Spin channels must have value 1 (non spin-polarised) or 2 (spin-polarised)"
        assert (type(self.number_of_kpoints) == int
                and self.number_of_kpoints > 0
                ), "The number of kpoints is not a positive integer"
        assert (type(self.number_of_bands) == int and self.number_of_bands > 0
                ), "The number of bands is not a positive integer"
        assert (type(self.number_of_ions) == int and self.number_of_ions > 0
                ), "The number of ions is not a positive integer"
        assert (self.CBM >
                self.VBM), "The CBM energy is lower than than the VBM energy"
        if self.fermi_energy < self.VBM:
            warnings.warn("The fermi energy is lower than the VBM")
        if self.fermi_energy > self.CBM:
            warnings.warn("The fermi energy is higher than the CBM")
        if ((self.occupancy == 0) | (self.occupancy == 1) |
            (self.occupancy == 2)).all() is False:
            warnings.warn("You have partial occupancy of bands")

    def parse_DOSCAR(self, filename='./DOSCAR'):
        """Parses the DOS and integrated DOS from a vasp DOSCAR file.

        Args:
            filename (str, optional): The location and filename of the DOSCAR to read in. Defaults to `'./DOSCAR'`.

        Returns:
            None.

        Notes:
            If the DOS has been sampled at more than 10000 points then this function will break at the expression for `num_data_points`.
            In this case, edit your DOSCAR file so that in the header there is a space preceding the number of points.
        """
        with open(filename, 'r') as f:
            lines = f.readlines()

        num_data_points = int(lines[5].split()[2])
        if len(lines[6].split()) == 5:
            self.dos = np.array([[
                float(x.split()[0]),
                float(x.split()[1]) + float(x.split()[2])
            ] for x in lines[6:num_data_points + 6]])
            self.integrated_dos = np.array([[
                float(x.split()[0]),
                float(x.split()[3]) + float(x.split()[4])
            ] for x in lines[6:num_data_points + 6]])
        elif len(lines[6].split()) == 3:
            self.dos = np.array([[float(x.split()[0]),
                                  float(x.split()[1])]
                                 for x in lines[6:num_data_points + 6]])
            self.integrated_dos = np.array(
                [[float(x.split()[0]),
                  float(x.split()[2])] for x in lines[6:num_data_points + 6]])
        else:
            print("problem parsing DOSCAR")
        return
