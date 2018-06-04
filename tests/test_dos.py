#! /usr/bin/env python3
import pytest
from effmass import dos

@pytest.mark.parametrize("data_object,VBM_index", [
	(pytest.lazy_fixture('MAPI_soc_data_object_with_DOSCAR'), 185),
	(pytest.lazy_fixture('toy_data_object'),6)
])
def test_find_VBM_dos_index(data_object, VBM_index):
	assert dos.find_dos_VBM_index(data_object) == VBM_index


@pytest.mark.parametrize("data_object,hole_fill_level,volume,concentration,VBM_index", [
	(pytest.lazy_fixture('MAPI_soc_data_object_with_DOSCAR'), -0.325845, 248.50, 1E20, 185),
	(pytest.lazy_fixture('toy_data_object'), -0.0275, 200, 1E20, 6)
])
def test_hole_fill_level(data_object, hole_fill_level, volume, concentration, VBM_index):
	assert dos.hole_fill_level(data_object, volume, concentration, VBM_index) == pytest.approx(hole_fill_level)