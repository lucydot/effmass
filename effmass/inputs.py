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


class Data_Vasp():
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

    def __init__(self, outcar_path, procar_path, ignore=0):
        r"""
        Initialises an instance of the :class:`~effmass.inputs.Data` class and checks data using :meth:`check_data`.

        Args:
            outcar_path (str): The path to the OUTCAR file
            procar_path (str): The path to the PROCAR file
            ignore (int): The number of kpoints to ignore at the beginning of the bandstructure slice through kspace (useful for hybrid calculations where zero weightings are appended to a previous self-consistent calculation).

        Returns:
            None.
        """
        assert (type(outcar_path) == str), "The OUTCAR path must be a string"
        assert (type(procar_path) == str), "The PROCAR path must be a string"
        assert (type(ignore) == int and ignore >= 0
                ), "The number of kpoints to ignore must be a positive integer"

        reciprocal_lattice = outcar.reciprocal_lattice_from_outcar(outcar_path)
        vasp_data = procar.Procar()
        vasp_data.read_from_file(procar_path)

        self.spin_channels = vasp_data.spin_channels
        self.number_of_bands = vasp_data.number_of_bands
        self.number_of_ions = vasp_data.number_of_ions

        number_of_kpoints = vasp_data.number_of_k_points
        if vasp_data.calculation['spin_polarised']: # to account for the change in PROCAR format for calculations with 2 spin channels (1 k-point block ---> 2 k-point blocks)
            energies = np.zeros([self.number_of_bands*2,number_of_kpoints]) # This is a very ugly way to slice 'n' dice. Should avoid creating new array and use array methods instead. But it does the job so will keep for now.
            for i in range(self.number_of_bands):
                energies[i] = vasp_data.bands[:, 1:].reshape(
                                            number_of_kpoints*2,  # factor or 2 for each kpoint block
                                            self.number_of_bands).T[i][:number_of_kpoints]
                energies[self.number_of_bands+i] = vasp_data.bands[:, 1:].reshape(
                                            number_of_kpoints*2,
                                            self.number_of_bands).T[i][number_of_kpoints:]
            occupancy = np.zeros([self.number_of_bands*2,number_of_kpoints])
            for i in range(self.number_of_bands):
                occupancy[i] = vasp_data.occupancy[:, 1:].reshape(
                                                 number_of_kpoints*2,
                                                 self.number_of_bands).T[i][:number_of_kpoints]
                occupancy[self.number_of_bands+i] = vasp_data.occupancy[:, 1:].reshape(
                                                 number_of_kpoints*2,
                                                 self.number_of_bands).T[i][number_of_kpoints:]
        else:
            energies = vasp_data.bands[:, 1:].reshape(
                                            number_of_kpoints,
                                            self.number_of_bands).T
            occupancy = vasp_data.occupancy[:, 1:].reshape(
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

        self.kpoints = vasp_data.k_points[ignore:vasp_data.number_of_k_points]
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


class Data_Aims():
    r"""
    Class for parsing and storing data from a FHI-AIMS calculation.

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
    """

    def __init__(self, file_path):
        r"""
        Initialises an instance of the :class:`~effmass.inputs.Data_Aims` class and checks data using :meth:`check_data`.

        Args:
            file_path (str): The path to the directory containing output, geometry.in, control.in and bandstructure files
            ignore (int): The number of kpoints to ignore at the beginning of the bandstructure slice through kspace (useful for hybrid calculations where zero weightings are appended to a previous self-consistent calculation).

        Returns:
            None.
        """
        assert (type(file_path) == str), "The file path must be a string"

        "Finding reciprocal lattice vectors"

        latvec = []

        for line in open("{}/geometry.in".format(file_path)):
            line = line.split("\t")[0]
            words = line.split()
            if len(words) == 0:
                continue
            if words[0] == "lattice_vector":
                if len(words) != 4:
                    raise Exception("geometry.in: Syntax error in line '"+line+"'")
                latvec.append(np.array(words[1:4]))

        if len(latvec) != 3:
            raise Exception("geometry.in: Must contain exactly 3 lattice vectors")

        latvec = np.asarray(latvec)
        latvec = latvec.astype(np.float)

        #Calculate reciprocal lattice vectors
        rlatvec = []
        volume = (np.dot(latvec[0,:],np.cross(latvec[1,:],latvec[2,:])))
        rlatvec.append(np.array(2*math.pi*np.cross(latvec[1,:],latvec[2,:])/volume))
        rlatvec.append(np.array(2*math.pi*np.cross(latvec[2,:],latvec[0,:])/volume))
        rlatvec.append(np.array(2*math.pi*np.cross(latvec[0,:],latvec[1,:])/volume))

        reciprocal_lattice = np.asarray(rlatvec)
        self.reciprocal_lattice = reciprocal_lattice


        "Finding spin channels"

        spin_channels = 0

        for line in open("{}/control.in".format(file_path)):
            line = line.split("\t")[0]
            words = line.split()
            if len(words) == 0:
                continue
            if words[0] == "spin" and words[1] == "none":
                spin_channels = 1
            elif words[0] == "spin" and words[1] == "collinear":
                spin_channels = 2
            elif words[0] == "include_spin_orbit":
                spin_channels = 4

        self.spin_channels = spin_channels

        "Finding number of bands"

        number_of_bands = 0

        for line in open("{}/calculation.out".format(file_path)):
            line = line.split("\t")[0]
            if "Number of states" in line:
                words = line.split()
                number_of_bands = int(words[-1])

                break

        if spin_channels == 2: #Doubling for spin-polarised calculation
            number_of_bands = 2*number_of_bands

        self.number_of_bands = number_of_bands

        "Finding number of ions"

        number_of_ions = 0

        for line in open("{}/calculation.out".format(file_path)):
            line = line.split("\t")[0]
            if "Number of atoms" in line:
                words = line.split()
                number_of_ions = int(words[-1])

                break

        self.number_of_ions = number_of_ions

        "Finding number of kpoints and determining number of BZ paths"

        number_of_kpoints = 0
        number_of_BZ_paths = 0
        path_list = []

        for line in open("{}/control.in".format(file_path)):
            line = line.split("\t")[0]
            if not line.startswith("#") and "output band" in line:
                words = line.split()
                path_list.append(int(words[8]))
                number_of_BZ_paths += 1

        number_of_kpoints = sum(path_list)

        "Reading out bandstructure files to determine kpoint, energy and occupation matrices"

        kpoints = np.zeros([number_of_kpoints,3])
        energies = np.zeros([number_of_bands,number_of_kpoints])
        occupancy = np.zeros([number_of_bands,number_of_kpoints])
        path_counter = 0

        if spin_channels == 1 or spin_channels == 4:
            kpoint_counter = 0
            while path_counter<number_of_BZ_paths:
                kpoint_counter = sum(path_list[:path_counter])
                for line in open("{}/band1{:03d}.out".format(file_path, path_counter+1)):
                    line = line.split("\t")[0]
                    words = line.split()
                    kpoints[kpoint_counter,0] = float(words[1])
                    kpoints[kpoint_counter,1] = float(words[2])
                    kpoints[kpoint_counter,2] = float(words[3])
                    for i in range(number_of_bands):
                        energies[i,kpoint_counter] = float(words[5+2*i])
                        occupancy[i,kpoint_counter] = float(words[4+2*i])
                    kpoint_counter += 1
                path_counter +=1

        if spin_channels == 2:
            while path_counter<number_of_BZ_paths:
                kpoint_counter = sum(path_list[:path_counter])
                for line in open("{}/band1{:03d}.out".format(file_path, path_counter+1)):
                    line = line.split("\t")[0]
                    words = line.split()
                    kpoints[kpoint_counter,0] = float(words[1])
                    kpoints[kpoint_counter,1] = float(words[2])
                    kpoints[kpoint_counter,2] = float(words[3])
                    for i in range(number_of_bands//2):
                        energies[i,kpoint_counter] = float(words[5+2*i])
                        occupancy[i,kpoint_counter] = float(words[4+2*i])
                    kpoint_counter += 1
                kpoint_counter = sum(path_list[:path_counter])
                for line in open("{}/band2{:03d}.out".format(file_path, path_counter+1)):
                    line = line.split("\t")[0]
                    words = line.split()
                    for i in range(number_of_bands//2):
                        energies[number_of_bands//2+i,kpoint_counter] = float(words[5+2*i])
                        occupancy[number_of_bands//2+i,kpoint_counter] = float(words[4+2*i])
                    kpoint_counter += 1
                path_counter += 1

        "Delete double kpoints at path edges"

        index_count = len(kpoints)
        index = 0
        while index < index_count-1:
            if np.array_equal(kpoints[index],kpoints[index+1]):
                kpoints = np.delete(kpoints,index+1,axis=0)
                energies = np.delete(energies,index+1,axis=1)
                occupancy = np.delete(occupancy,index+1,axis=1)
                index_count = len(kpoints)
            index += 1

        self.number_of_kpoints = len(kpoints)


        self.CBM = extrema._calc_CBM(occupancy, energies)
        self.VBM = extrema._calc_VBM(occupancy, energies)
        self.fermi_energy = (self.CBM + self.VBM) / 2

        "Cutting energy values in a range of 30 eV above and below the Fermi level. FHI AIMS is all electron, but not all states are needed for a meaningful effmass calculation"

        index_count = len(occupancy)
        index = 0
        while index < index_count-1:
            if all(item < self.fermi_energy - 30 for item in energies[index]):
                energies = np.delete(energies, index, axis = 0)
                occupancy = np.delete(occupancy, index, axis = 0)
                index_count = len(occupancy)
            elif all(item > self.fermi_energy + 30 for item in energies[index]):
                energies = np.delete(energies, index, axis = 0)
                occupancy = np.delete(occupancy, index, axis = 0)
                index_count = len(occupancy)
            else:
                index += 1

        self.energies = energies
        self.occupancy = occupancy
        self.kpoints = kpoints
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

