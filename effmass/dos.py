#! /usr/bin/env python3
"""
A module for analysing DOSCAR data.
"""

import numpy as np


def _check_integrated_dos_loaded(Data):
    """Helper function to check if :attr:`~effmass.inputs.Data.integrated_dos`
    is loaded.

    Args:
        Data (Data): Instance of the :class:`Data` class.

    Returns:
        None.
    """
    assert Data.integrated_dos != [], "Data.integrated_dos is empty. Please set attribute, perhaps using Data.parse_DOSCAR, and try again."
    return


def _check_dos_loaded(Data):
    """Helper function to check if :attr:`~effmass.inputs.Data.dos` is loaded.

    Args:
        Data (Data): Instance of the :class:`Data` class.

    Returns:
        None.
    """
    assert Data.integrated_dos != [], "Data.dos is empty. Please set attribute, perhaps using Data.parse_DOSCAR, and try again."
    return


def find_dos_VBM_index(Data):
    """Finds the lowest index of the
    :attr:`~effmass.inputs.Data.integrated_dos` array where the energy exceeds
    :attr:`~effmass.inputs.Data.VBM`.

    Args:
        Data (Data): Instance of the Data class.

    Returns:
        int: the lowest index of the :attr:`~effmass.inputs.Data.integrated_dos` array where the energy exceeds :attr:`~effmass.inputs.Data.VBM`.
    """
    _check_integrated_dos_loaded(Data)

    for i in range(len(Data.integrated_dos)):
        if Data.VBM < Data.integrated_dos[i][0]:
            if Data.dos[i][1] == 0.0:
                return i


def find_dos_CBM_index(Data):
    """Finds the highest index of the
    :attr:`~effmass.inputs.Data.integrated_dos` array where the energy is less
    than :attr:`~effmass.inputs.Data.CBM`.

    Args:
        Data (Data): Instance of the Data class.

    Returns:
        int: the highest index of the :attr:`~effmass.inputs.Data.integrated_dos` array where the energy is less than :attr:`~effmass.inputs.Data.CBM`.
    """
    _check_integrated_dos_loaded(Data)

    for i in range(len(Data.integrated_dos))[::-1]:
        if Data.CBM > Data.integrated_dos[i][0]:
            if Data.dos[i][1] == 0.0:
                return i


def electron_fill_level(Data, volume, concentration, CBM_index):
    r"""
    Finds the energy to which a given electron concentration will fill the density of states in :attr:`~effmass.inputs.Data.integrated_dos`.
    
    Uses linear interpolation to estimate the energy between two points given in the DOSCAR.

    Args:
        Data (Data): Instance of the :class:`Data` class. 
        volume (float): volume of the unit cell in angstrom :math:`^3`.
        concentration (float): electron concentration in cm :math:`^{-3}`.
        CBM_index (int): highest index of the :attr:`~effmass.inputs.Data.integrated_dos` array where the energy is less than :attr:`~effmass.inputs.Data.CBM`.

    Returns:
        float: the energy (eV, referenced from the CBM) to which the electrons will fill. For the case where the concentration specified would fill all states specified by :attr:`~effmass.inputs.Data.integrated_dos`, None is returned.

    Notes:
        The precision of the result will depend upon the energy resolution in the DOSCAR.
    """
    _check_integrated_dos_loaded(Data)
    CBM_index = CBM_index
    states_per_unit_cell = volume * 1E-30 * concentration * 1E6
    assert (
        states_per_unit_cell < np.absolute(Data.integrated_dos[-1][1] -
                                           Data.integrated_dos[CBM_index][1])
    ), "the concentration specified would fill all available energy states"

    upper_index = len(Data.integrated_dos) - 1
    lower_index = CBM_index

    # this function is made a little more complicated because the dos can be a step function
    # therefore, we cannot interpolate between consecutive indices, but need to find the range of the step.

    for i in range(CBM_index + 1, len(Data.integrated_dos)):
        if states_per_unit_cell < np.absolute(
                Data.integrated_dos[i][1] - Data.integrated_dos[CBM_index][1]):
            upper_index = i  # marks the maximum energy for this concentration
            break

    for i in range(1, upper_index - CBM_index):
        if Data.integrated_dos[upper_index -
                               1][1] - Data.integrated_dos[upper_index - 1 -
                                                           i][1] != 0:
            lower_index = upper_index - i  # marks the minimum energy for this concentration
            break

    # linear interpolation
    proportion = (states_per_unit_cell - (Data.integrated_dos[lower_index][1] -
                                          Data.integrated_dos[CBM_index][1])
                  ) / np.absolute(Data.integrated_dos[upper_index][1] -
                                  (Data.integrated_dos[lower_index][1]))
    energy = ((Data.integrated_dos[lower_index][0] -
               Data.integrated_dos[CBM_index][0]) +
              (Data.integrated_dos[upper_index][0] -
               Data.integrated_dos[lower_index][0]) * proportion)

    # where lower index has been calculated to be below the CBM, as the concentration is smaller than the first step so set lower energy equal to CBM (0eV)
    if Data.integrated_dos[lower_index][0] - Data.integrated_dos[CBM_index][0] < 0:
        energy = (Data.integrated_dos[upper_index][0] -
                  Data.integrated_dos[CBM_index][0]) * proportion

    return energy


def hole_fill_level(Data, volume, concentration, VBM_index):
    r"""
    Finds the energy to which a given hole concentration will fill the density of states in :attr:`~effmass.inputs.Data.integrated_dos`.
    
    Uses linear interpolation to estimate the energy between two points given in the DOSCAR.

    Args:
        Data (Data): Instance of the :class:`Data` class. 
        volume (float): volume of the unit cell in angstrom :math:`^3`.
        concentration: hole concentration in cm :math:`^{-3}`.
        VBM_index (int): lowest index of the :attr:`~effmass.inputs.Data.integrated_dos` array where the energy is more than than :attr:`~effmass.inputs.Data.VBM`.

    Returns:
        float: the energy (eV, referenced from the VBM) to which the holes will fill. For the case where the concentration specified would fill all states specified by :attr:`~effmass.inputs.Data.integrated_dos`, None is returned.

    Notes:
        The precision of the result will depend upon the energy resolution in the DOSCAR. 
    """
    _check_integrated_dos_loaded(Data)
    VBM_index = VBM_index
    states_per_unit_cell = volume * 1E-30 * concentration * 1E6
    assert (
        states_per_unit_cell < np.absolute(Data.integrated_dos[0][1] -
                                           Data.integrated_dos[VBM_index][1])
    ), "the concentration specified would fill all available energy states"

    upper_index = 0
    lower_index = VBM_index

    # this function is made a little more complicated because the dos can be a step function
    # therefore, we cannot interpolate between consecutive indices, but need to find the range of the step.

    for i in range(0, VBM_index)[::-1]:
        if states_per_unit_cell < np.absolute(
                Data.integrated_dos[i][1] -
                Data.integrated_dos[VBM_index - 1][1]):
            upper_index = i  # marks the maximum hole energy for this concentration
            break

    for i in range(1, VBM_index):
        if Data.integrated_dos[upper_index +
                               1][1] - Data.integrated_dos[upper_index + 1 +
                                                           i][1] != 0:
            lower_index = upper_index + i  # marks the minimum hole energy for this concentration
            break

    # linear interpolation
    proportion = (states_per_unit_cell -
                  np.absolute(Data.integrated_dos[lower_index][1] -
                              Data.integrated_dos[VBM_index][1])
                  ) / np.absolute(Data.integrated_dos[upper_index][1] -
                                  (Data.integrated_dos[lower_index][1]))
    energy = ((Data.integrated_dos[lower_index][0] -
               Data.integrated_dos[VBM_index][0]) -
              (Data.integrated_dos[lower_index][0] -
               Data.integrated_dos[upper_index][0]) * proportion)

    # where lower index has been calculated to be above the CBM, as the concentration is smaller than the first step and so set lower energy equal to VBM (0eV)
    if Data.integrated_dos[lower_index][0] - Data.integrated_dos[VBM_index][0] > 0:
        energy = (Data.integrated_dos[upper_index][0] -
                  Data.integrated_dos[VBM_index][0]) * proportion

    return energy
