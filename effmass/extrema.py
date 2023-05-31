#! /usr/bin/env python3

"""
A module for finding the band structure extrema and instantiating 
a :class:`Segment` object for each extrema point.

The extrema are found within an energy range given by the :class:`Settings` class.
Each `Segment` object contains data for kpoints within an energy range given by the :class:`Settings` class. 
"""

import math
import warnings
import numpy as np
import numpy.ma as ma
from effmass import analysis
from effmass import inputs


def _check_partial_occupancy(occupancy):
    """Raises warning if there is partial occupancy of bands.

    Args:
        occupancy (array): Array with shape (number_of_bands,number_of_kpoints). Each row contains occupation number of the eigenstates for a particular band. Values range from 0-1 (spin-polarised) or 0-2 (non-spin-polarised). See :attr:`effmass.inputs.Data.occupancy`.

    Returns:
        None.
    """
    if not np.all(np.in1d(occupancy,([0, 1, 2]))):
        warnings.warn(
            "You have partial occupation numbers in Data.occupancy. You should check that the attributes Data.VBM, Data.CBM and Data.fermi_energy are correct, and if not, set them manually."
        )


def _calc_CBM(occupancy, energy):
    """Finds the minimum unoccupied energy eigenstate (conduction band minimum)
    .

    Args:
        occupancy (array): occupancy of eigenstates with shape (number_of_kpoints, number_of_bands)
        energy (array): energy (eV) of eigenstates with shape (number_of_kpoints, number_of_bands)

    Returns:
        float: the minimum unoccupied energy (eV)
    """
    _check_partial_occupancy(occupancy)
    energy_unoccupied = ma.masked_where(occupancy > 0.5, energy)
    return np.amin(energy_unoccupied)


def _calc_VBM(occupancy, energy):
    """Finds the minimum unoccupied energy eigenstate (valence band maximum).

    Args:
        occupancy (array): occupancy of eigenstates with shape (number_of_kpoints, number_of_bands)
        energy (array): energy (eV) of eigenstates with shape (number_of_kpoints, number_of_bands)

    Returns:
        float: the minimum unoccupied energy (eV)
    """
    _check_partial_occupancy(occupancy)
    energy_occupied = ma.masked_where(occupancy < 0.5, energy)
    return np.amax(energy_occupied)


def calculate_direction(a, b):
    """Calculates the direction vector between two points.

    Args:
        a (list): the position vector of point a.
        b (list): the position vector of point b.

    Returns:
        array: The (unnormalised) direction vector between points a and b. The smallest magnitude of an element is 1 (eg: [1,1,2]).
    """
    difference = np.subtract(a, b)
    if np.count_nonzero(difference) < 1:
        print("The two k-points are equal")
        return np.array([0, 0, 0])

    # we need to find the smallest non-zero value within a-b
    a = np.array(a)
    b = np.array(b)
    direction_masked = ma.masked_equal(
        a - b, 0)  # return array with invalid entries where values are equal
    direction_filled = ma.filled(
        direction_masked, 10
        **6)  # fill invalid elements of array with a large number s
    direction_absolute = np.absolute(
        direction_filled)  # return absolute values of each element
    smallest = np.amin(direction_absolute)
    direction = (
        b - a) / smallest  # use the minimum absolute value as a divisor a-b
    if -1 in direction:
        direction = np.multiply(direction, -1)
    return direction


def change_direction(kpoints, kpoint_indices):
    """Finds the index of the kpoint (if any) where there is a change of
    direction in reciprocal space.

    Args:
        kpoints (array): array of kpoints with shape (number_of_kpoints, 3). Each row contains the fractional coordinates of a kpoint [kx,ky,kz]. See :attr:`effmass.inputs.Data.kpoints`.
        kpoint_indices (list (int)): the kpoint indices over which to check for change in direction

    Returns:
        int: the index of the kpoint where there is a change of direction. If there is no change of direction, returns None.
    """
    new_direction = calculate_direction(kpoints[kpoint_indices[0]],
                                        kpoints[kpoint_indices[1]])
    change_index = None
    for i, value in enumerate(kpoint_indices[:-1]):
        old_direction = new_direction
        new_direction = calculate_direction(kpoints[kpoint_indices[i]],
                                            kpoints[kpoint_indices[i + 1]])
        if np.linalg.norm(new_direction - old_direction) > 0.005:
            change_index = i + 1
            break
    return change_index

def calc_CBM_VBM_from_Fermi(Data, CBMVBM_search_depth=4.0):
    """ This function is used to find the CBM and VBM when there is no occupancy data. It relies upon the Fermi level being in the middle of the band gap.
    The CBMVBM_search_depth is refereced from the fermi energy.
    
    Args:
        DataASE (DataASE): instance of the :class:`DataASE` class.

    Returns:
        (float, float): A tuple containing the conduction band minimum and valence band maximum in eV.
"""
    Data.CBM = Data.fermi_energy
    Data.VBM = Data.fermi_energy

    Settings = inputs.Settings(extrema_search_depth=CBMVBM_search_depth)

    CB_indices = find_CB_indices(Data, Settings)
    VB_indices = find_VB_indices(Data, Settings)

    CBM = min([Data.energies[i][j] for i,j in CB_indices])
    VBM = max([Data.energies[i][j] for i,j in VB_indices])

    return CBM, VBM

def get_minimum_indices(Data,extrema_search_depth):
    """Finds the kpoint indices and band indices of all minimum turning points in CB within `extrema_search_depth`.

    Args:
        Data (Data): instance of the :class:`Data` class.
        extrema_search_depth (float): energy in kT from bandedge over which to search for minima.

    Returns:
        array: A 2-dimensional array. Each row contains [:attr:`efmmas.inputs.Data.bands` index, :attr:`effmass.inputs.Data.kpoints` index] for each minimum in the CB band.
    """
    energy_CB = ma.masked_where(Data.energies < Data.CBM, Data.energies)
    # returns array of energies within energy_range
    electrons_in_range = ma.masked_where(
        np.absolute(energy_CB - Data.CBM) > extrema_search_depth, energy_CB)
    CB_minima = ma.masked_where(
        _mark_minima(electrons_in_range) == 0, electrons_in_range)
    CB_min_indices = np.argwhere(CB_minima.mask == 0)
    return CB_min_indices

def get_maximum_indices(Data,extrema_search_depth):
    """Finds the kpoint indices and band indices of all maximum turning points in VB within `extrema_search_depth`.

    Args:
        Data (Data): instance of the :class:`Data` class.
        extrema_search_depth (float): energy in kT from bandedge over which to search for maxima.

    Returns:
        array: A 2-dimensional array. Each row contains [:attr:`efmmas.inputs.Data.bands` index, :attr:`effmass.inputs.Data.kpoints` index] for each maximum in the VB band.
    """
    energy_VB = ma.masked_where(Data.energies > Data.VBM, Data.energies)
    # returns array of energies within energy_range
    holes_in_range = ma.masked_where(
        np.absolute(energy_VB - Data.VBM) >
        extrema_search_depth, energy_VB)
    VB_maxima = ma.masked_where(
        _mark_maxima(holes_in_range) == 0, holes_in_range)
    VB_max_indices = np.argwhere(VB_maxima.mask == 0)
    return VB_max_indices

def get_frontier_CB_indices(Data,CB_min_indices, degeneracy_condition):
    """Returns the indices of the lowest energy minima across the Brillouin Zone

    Args:
        Data (Data): instance of the :class:`Data` class.
        CB_min_indices (array(int)): A 2-dimensional array. Each row contains [:attr:`efmmas.inputs.Data.bands` index, :attr:`effmass.inputs.Data.kpoints` index] for each minimum in the CB band.

    Returns:
        frontier_indices (array(int)): A 2-dimensional array. Each row contains [:attr:`efmmas.inputs.Data.bands` index, :attr:`effmass.inputs.Data.kpoints` index] for each minimum in the frontier conduction band(s).
    """
    frontier_indices = np.where(CB_min_indices[:, 0] == CB_min_indices[:, 0].min())[0]
    frontier_indices = CB_min_indices[frontier_indices]
    for band, kpoint in frontier_indices:
        i = 1
        while math.isclose(Data.energies[band+i, kpoint],Data.energies[band, kpoint], abs_tol=degeneracy_condition):
            frontier_indices = np.append(frontier_indices, np.array([[band+i, kpoint]]), axis=0)
            i += 1
    return frontier_indices

def get_frontier_VB_indices(Data,VB_max_indices, degeneracy_condition):
    """Returns the indices of the highest energy maxima across the Brillouin Zone

    Args:
        Data (Data): instance of the :class:`Data` class.
        VB_max_indices (array(int)): A 2-dimensional array. Each row contains [:attr:`efmmas.inputs.Data.bands` index, :attr:`effmass.inputs.Data.kpoints` index] for each maximum in the VB band.

    Returns:
        frontier_indices (array(int)): A 2-dimensional array. Each row contains [:attr:`efmmas.inputs.Data.bands` index, :attr:`effmass.inputs.Data.kpoints` index] for each maximum in the frontier valence band(s).
    """
    frontier_indices = np.where(VB_max_indices[:, 0] == VB_max_indices[:, 0].max())[0]
    frontier_indices = VB_max_indices[frontier_indices]
    for band, kpoint in frontier_indices:
        i = 1
        while math.isclose(Data.energies[band-i, kpoint],Data.energies[band, kpoint], abs_tol=degeneracy_condition):
            frontier_indices = np.append(frontier_indices, np.array([[band-i, kpoint]]), axis=0)
            i += 1
    return frontier_indices

def find_CB_indices(Data, Settings):
    """Finds the kpoint index and band index of the minimum energy
    turning points within :attr:`effmass.inputs.Settings.energy_range` of the
    conduction band minimum (:attr:`effmass.inputs.Data.CBM`). Return indices for the 
    lowest energy CB only if `frontier_bands_only` is True.

    Args:
        Data (Data): instance of the :class:`Data` class.
        Settings (Settings): instance of the :class:`Settings` class.

    Returns:
        array: A 2-dimensional array. Contains [:attr:`efmmas.inputs.Data.bands` index, :attr:`effmass.inputs.Data.kpoints` index] for each minima.
    """
    if Settings.conduction_band is True:
        CB_min_indices = get_minimum_indices(Data, Settings.extrema_search_depth)

        if Settings.frontier_bands_only is True:
            # At a given k-point there may be multiple bands with highest energy
            CB_min_indices = get_frontier_CB_indices(Data, CB_min_indices, Settings.degeneracy_condition)

    else:
        CB_min_indices = None

    return CB_min_indices

def find_VB_indices(Data, Settings):
    """Finds the kpoint index and band index of the maximum energy
    turning points within :attr:`effmass.inputs.Settings.energy_range` of the
    valence band maximum (:attr:`effmass.inputs.Data.VBM`). Return indices for the 
    highest energy VB only if `frontier_bands_only` is True.

    Args:
        Data (Data): instance of the :class:`Data` class.
        Settings (Settings): instance of the :class:`Settings` class.

    Returns:
        array: A 2-dimensional array. Contains [:attr:`efmmas.inputs.Data.bands` index, :attr:`effmass.inputs.Data.kpoints` index] for each maxima.
    """
    
    VB_max_indices = get_maximum_indices(Data, Settings.extrema_search_depth)

    if Settings.frontier_bands_only is True:
        # At a given k-point there may be multiple bands with highest energy
        VB_max_indices = get_frontier_VB_indices(Data, VB_max_indices, Settings.degeneracy_condition)


    return VB_max_indices

    # Returns an array of band numbers and k-points for extrema, selected according to user settings.
    # The first index differentiates between the valence band and conduction band.
    


def _mark_maxima(holes_array):
    """Helper function which takes an array of hole energies and returns an
    array for masking those which are not maxima."""
    # create array to append to
    not_maxima = []
    width = holes_array.shape[1]
    height = holes_array.shape[0]

    # handle edge cases
    for band in range(0, height):
        if (holes_array[band, 1] >= holes_array[band, 0]):
            not_maxima.append([band, 0])

    for band in range(0, height):
        if (holes_array[band, -2] >= holes_array[band, -1]):
            not_maxima.append([band, -1])

    # find the bands where numbers either side aren't bigger (the maxima)
    # will also include inflexion but needed to get HSP's
    for band in range(0, height):
        for k in range(1, width - 1):
            if (holes_array[band, k - 1] >= holes_array[band, k]
                    or holes_array[band, k + 1] >= holes_array[band, k]):
                not_maxima.append([band, k])

    # Need to assign false after inspection
    # or it screws with the algorithm
    for row in not_maxima:
        holes_array[row[0], row[1]] = False

    # returns matrix with all non-maxima marked false
    return holes_array


def _mark_minima(electrons_array):
    """Helper function which takes an array of electron energies and returns an
    array for masking those which are not minima."""
    # create array to append to
    not_minima = []
    width = electrons_array.shape[1]
    height = electrons_array.shape[0]

    # handle edge cases
    for band in range(0, height):
        if (electrons_array[band, 1] <= electrons_array[band, 0]):
            not_minima.append([band, 0])

    for band in range(0, height):
        if (electrons_array[band, -2] <= electrons_array[band, -1]):
            not_minima.append([band, -1])

    # find the bands where numbers either side aren't bigger (the maxima)
    for band in range(0, height):
        for k in range(1, width - 1):
            if (electrons_array[band, k - 1] <= electrons_array[band, k] or
                    electrons_array[band, k + 1] <= electrons_array[band, k]):
                not_minima.append([band, k])

    # Need to assign false after inspection
    # or it screws with the algorithm
    for row in not_minima:
        electrons_array[row[0], row[1]] = False

    # returns matrix with all non-maxima marked false
    return electrons_array


def _within(band_index, kpoint_index, trial_kpoint_index, Settings, Data):
    """Helper function which checks whether trial_kpoint_index is less than
    :attr:`effmass.inputs.Settings.energy_range`.

    away from kpoint_index for a given band
    """
    if (np.absolute(Data.energies[band_index, trial_kpoint_index] -
                    Data.energies[band_index, kpoint_index])
        ) < Settings.energy_range:
        return True
    else:
        return False


def _kpoints_after(band_index, kpoint_index, Settings, Data):
    """Helper function that returns eigenstates which are 1) on the same band
    2) before the given kpoint in the bandstructure route and 3) within
    :attr:`effmass.inputs.Settings.energy_range`."""
    trial_kpoint_index = kpoint_index + 1
    kpoint_after = [kpoint_index]
    while (trial_kpoint_index < Data.number_of_kpoints
           and _within(band_index, kpoint_index, trial_kpoint_index, Settings,
                       Data) is True):
        kpoint_after.append(trial_kpoint_index)
        trial_kpoint_index = trial_kpoint_index + 1

    return kpoint_after


def _kpoints_before(band_index, kpoint_index, Settings, Data):
    """Helper function that returns eigenstates which are 1) on the same band
    2) before the given eigenstate in the route through reciprocal space 3)
    within :attr:`effmass.inputs.Settings.energy_range`."""
    trial_kpoint_index = kpoint_index - 1
    kpoint_before = [kpoint_index]
    while trial_kpoint_index >= 0 and _within(band_index, kpoint_index,
                                              trial_kpoint_index, Settings,
                                              Data) is True:
        kpoint_before.append(trial_kpoint_index)
        trial_kpoint_index = trial_kpoint_index - 1

    return kpoint_before


def get_kpoints_before(band_index,
                       kpoint_index,
                       Settings,
                       Data,
                       truncate_dir_change=True):
    """For a given eigenstate, finds eigenstates which 1) belong to the same
    band 2) come before the given eigenstate in the route through reciprocal
    space 3) are within :attr:`effmass.inputs.Settings.energy_range`.

    Args:
        band_index (int): index of :attr:`effmass.inputs.Data.bands`.
        kpoint_index (int): index of :attr:`effmass.inputs.Data.kpoints`.
        Settings (Settings): instance of the :class:`Settings` class.
        Data (Data): instance of the :class:`Data` class.
        truncate_dir_change (bool): If True, truncates eigenstates when there is a change in direction. If False, there is no truncation. Defaults to True.

    Returns:
        list(int): indices of :attr:`effmass.inputs.Data.kpoints`.
    """
    kpoints_before = _kpoints_before(band_index, kpoint_index, Settings, Data)

    if len(kpoints_before) > 2:
        change_index = change_direction(Data.kpoints, kpoints_before)
        if change_index:
            if truncate_dir_change is True:
                kpoints_before = kpoints_before[:change_index]
                if len(kpoints_before) < 3:
                    kpoints_before = None
    else:
        kpoints_before = None

    return kpoints_before


def get_kpoints_after(band_index,
                      kpoint_index,
                      Settings,
                      Data,
                      truncate_dir_change=True):
    """For a given eigenstate, finds eigenstates which 1) belong to the same
    band 2) come after the given eigenstate in the route through reciprocal
    space 3) are within :attr:`effmass.inputs.Settings.energy_range`.

    Args:
        band_index (int): index of :attr:`effmass.inputs.Data.bands`.
        kpoint_index (int): index of :attr:`effmass.inputs.Data.kpoints`.
        Settings (Settings): instance of the :class:`Settings` class.
        Data (Data): instance of the :class:`Data` class.
        truncate_dir_change (bool): If True, truncates eigenstates when there is a change in direction. If False, there is no truncation. Defaults to True.

    Returns:
        list(int): indices of :attr:`effmass.inputs.Data.kpoints`.
    """
    kpoints_after = _kpoints_after(band_index, kpoint_index, Settings, Data)

    # now make sure that there is not a change in direction for this set of kpoints and that long enough
    if len(kpoints_after) > 2:
        change_index = change_direction(Data.kpoints, kpoints_after)
        if change_index:
            if truncate_dir_change is True:
                kpoints_after = kpoints_after[:change_index]
                if len(kpoints_after) < 3:
                    kpoints_after = None
    else:
        kpoints_after = None

    return kpoints_after


def _dot_product(vector1, vector2):
    return np.dot(vector1, vector2) / (np.linalg.norm(vector1)*np.linalg.norm(vector2))


def filter_segments_by_direction(segment_list, direction):
    """Filter a list of Segments so that only those in a particular direction remain.
    
    Args:
        segment_list (list(Segment)):  A list of :class:`Segment` objects.
        direction (array(float)): The direction array, length 3.
         
    Returns:
        segment_list (list(Segment)):  A list of :class:`Segment` objects.
    """
    
    return [segment for segment in segment_list if math.isclose(_dot_product(segment.direction, direction), 1)]    

    # need to test for equality of direction not magnitude of vector


def generate_segments(Settings, Data, bk=None, truncate_dir_change=True):
    """Generates a list of Segment objects.

    Args:
        Settings (Settings): instance of the :class:`Settings` class.
        Data (Data): instance of the :class:`Data` class.
        truncate_dir_change (bool): If True, truncates eigenstates when there is a change in direction. If False, there is no truncation. Defaults to True.
        bk (list(int)): To manually set an extrema point, in format [:attr:`effmass.inputs.Data.energies` row index, :attr:`effmass.inputs.Data.kpoints` row index]. Defaults to None.
   
   Returns:
        list(Segment): A list of :class:`Segment` objects.
    """
    if bk:
        extrema_array = bk

    else:
        if (Settings.valence_band is True)and (Settings.conduction_band is True):
            extrema_array = np.concatenate((find_VB_indices(Data, Settings),find_CB_indices(Data, Settings)))
        elif Settings.valence_band is True:
            extrema_array = find_VB_indices(Data, Settings)
        elif Settings.conduction_band is True:
            extrema_array = find_CB_indices(Data, Settings)
    
    kpoints_list = []
    band_list = []
    for band, kpoint in extrema_array: # flattened CB and VB arrays to a single array
        kpoints_before = get_kpoints_before(
            band,
            kpoint,
            Settings,
            Data,
            truncate_dir_change=truncate_dir_change)
        kpoints_after = get_kpoints_after(
            band,
            kpoint,
            Settings,
            Data,
            truncate_dir_change=truncate_dir_change)
        if kpoints_before:
            kpoints_list.append(kpoints_before)
            band_list.append(band)
        if kpoints_after:
            kpoints_list.append(kpoints_after)
            band_list.append(band)
    segment_list = [
        analysis.Segment(Data, band, kpoints)
        for band, kpoints in zip(band_list, kpoints_list)
    ]

    if Settings.direction:
        segment_list = filter_segments_by_direction(segment_list,np.array(Settings.direction))

    if len(segment_list) == 0:
        print("No segments found for current criteria in Settings")

    return segment_list
