#! /usr/bin/env python3
import pytest
from effmass import extrema
import numpy as np

def test_calc_CBM(toy_data_object):
    
    assert extrema._calc_CBM(toy_data_object.occupancy,toy_data_object.energies) == 1

def test_calc_VBM(toy_data_object):

    assert extrema._calc_VBM(toy_data_object.occupancy,toy_data_object.energies) == 0.5

def test_calculate_direction(toy_data_object):
   
    a = toy_data_object.kpoints[0]
    b = toy_data_object.kpoints[1]

    assert np.allclose(extrema.calculate_direction(a,b),np.array([1,0,0]))

def test_calculate_direction_opposite(toy_data_object):
    a = toy_data_object.kpoints[1]
    b = toy_data_object.kpoints[0]

    assert np.allclose(extrema.calculate_direction(a,b),np.array([1,0,0]))

def test_change_direction(toy_data_object):
    
    assert extrema.change_direction(toy_data_object.kpoints, [0,1,2,3,4,5]) == 3

def test_find_extrema_indices(toy_data_object,toy_settings_object):

    assert np.array_equal(extrema.find_extrema_indices(toy_data_object, toy_settings_object),np.array([[0,0],[0,3],[1,0],[1,4]]))

def test_get_kpoints_before(toy_data_object,toy_settings_object):

    assert extrema.get_kpoints_before(1, 4, toy_settings_object, toy_data_object, truncate_dir_change=True) == [4,3,2]

def test_get_kpoints_before_smaller_range(toy_data_object,toy_settings_object_smaller_range):

    assert extrema.get_kpoints_before(1, 4, toy_settings_object_smaller_range, toy_data_object, truncate_dir_change=True) == None

def test_get_kpoints_before_no_truncate(toy_data_object,toy_settings_object):
 
    assert extrema.get_kpoints_before(1, 4, toy_settings_object, toy_data_object, truncate_dir_change=False) ==[4,3,2,1,0]

def test_generate_segments(toy_data_object,toy_settings_object):  

    assert len(extrema.generate_segments(toy_settings_object, toy_data_object, truncate_dir_change=True)) == 3
    assert extrema.generate_segments(toy_settings_object, toy_data_object, truncate_dir_change=True)[-1].kpoint_indices == [4,3,2]
    assert extrema.generate_segments(toy_settings_object, toy_data_object, truncate_dir_change=True)[1].kpoint_indices == [0,1,2]
    assert extrema.generate_segments(toy_settings_object, toy_data_object, truncate_dir_change=True)[0].kpoint_indices == [0,1,2]
