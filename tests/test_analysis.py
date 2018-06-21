import pytest
from effmass import analysis, extrema, angstrom_to_bohr
import numpy as np
import math

def test_check_poly_order():
    with pytest.raises(AssertionError):
        analysis._check_poly_order(1)

def test_weighted_mean():
    values = [1,1,1,3,1.5,1.5,0.5]
    weighting = [2,2,2,1,1,1,1]
    assert math.isclose(analysis.weighted_mean(values,weighting),1.25)

def test_solve_quadratic():
	a = 2
	b = 11
	c = [-21,-40]
	assert math.isclose(analysis.solve_quadratic(a,b,c),[1.5,2.5])

def test_fermi_function(MAPI_soc_segment_object_hole):
	assert math.isclose(MAPI_soc_segment_object_hole.fermi_function(0,0),0.5)
	assert MAPI_soc_segment_object_hole.fermi_function(0,1) < 0.5
	assert MAPI_soc_segment_object_hole.fermi_function(1,0) > 0.5

def test_fermi_function_toy(toy_segments):
    online_calc = 1 - 0.999970880324388
    assert math.isclose(toy_segments[0].fermi_function(0.48,0.75),online_calc,rel_tol=1E-5)
    online_calc = 1 - 0.999986565914777
    assert math.isclose(toy_segments[0].fermi_function(0.46,0.75),online_calc,rel_tol=1E-5)

def test_fermi_function(MAPI_soc_segment_object_electron):
	assert math.isclose(MAPI_soc_segment_object_electron.fermi_function(0,0),0.5)
	assert MAPI_soc_segment_object_electron.fermi_function(0,1) > 0.5
	assert MAPI_soc_segment_object_electron.fermi_function(1,0) < 0.5

def test_return_ptype(MAPI_soc_segment_object_hole):
	assert MAPI_soc_segment_object_hole.return_ptype() == "hole"

def test_return_ptype(MAPI_soc_segment_object_electron):
	assert MAPI_soc_segment_object_electron.return_ptype() == "electron"

def test_return_bandedge_energy(toy_segments, toy_data_object):
	assert math.isclose(toy_segments[0].return_bandedge_energy(toy_data_object),extrema.calc_VBM(toy_data_object.occupancy,toy_data_object.energies))
	assert math.isclose(toy_segments[1].return_bandedge_energy(toy_data_object),extrema.calc_CBM(toy_data_object.occupancy,toy_data_object.energies))
	assert math.isclose(toy_segments[2].return_bandedge_energy(toy_data_object),extrema.calc_CBM(toy_data_object.occupancy,toy_data_object.energies))

def test_return_weighting(toy_segments):
    weighting_1 = (1-0.999955005754945) # by hand
    weighting_2 = (1-0.999983701573141)
    effmass_result = toy_segments[0].return_weighting()
    ratio_1 = weighting_1/weighting_2
    ratio_2 = effmass_result[1]/effmass_result[2]
    assert math.isclose(ratio_1,ratio_2,rel_tol=1E-5)

def test_calc_poly_first_derivative(toy_segments):
    assert len(set(toy_segments[0].calc_poly_derivatives(2, False, dk=toy_segments[0].dk_bohr)[1])) <= 1 # for quadratic all second derivatives should be equal
    toy_segments[0].dk_bohr=np.array([0,0.1,0.2]) # change to simple values
    toy_segments[0].dE_hartree=np.array([0,0.005,0.02])
    assert math.isclose(toy_segments[0].calc_poly_derivatives(polyfit_order=2, polyfit_weighting=False)[0][-1],0.2) # first derivative

def test_calc_poly_second_derivative(toy_segments):
    toy_segments[0].dk_bohr=np.array([0,0.1,0.2]) # change to simple values
    toy_segments[0].dE_hartree=np.array([0,0.005,0.02])
    assert math.isclose(toy_segments[0].calc_poly_derivatives(polyfit_order=2, polyfit_weighting=False)[1][-1],2*0.5) # second derivative (factor 2 as to the power 2)
  
def test_dot_product(toy_segments):
	assert np.isclose(toy_segments[0].cartesian_kpoints,toy_segments[0].kpoints)
	dk_bohr_by_hand = np.array([0,0.25/1.88973,0.5/1.88973])
	assert np.isclose(toy_segments[0].dk_bohr, dk_bohr_by_hand)

def test_return_polyfit(toy_segments):
    toy_segments[0].dk_bohr=np.array([0,0.1,0.2]) # change to simple values
    toy_segments[0].dE_hartree=np.array([0,0.005,0.02])
    assert math.isclose(toy_segments[0].return_polyfit(2, False)[0],0.0,abs_tol=1E-9)
    assert math.isclose(toy_segments[0].return_polyfit(2, False)[-1],0.02)

def test_calc_finite_difference_effmass(toy_segments):
 	toy_segments[0].dk_bohr=np.array([0,0.1,0.2]) # change to simple values
 	toy_segments[0].dE_hartree=np.array([0,0.005,0.02])
 	assert math.isclose(toy_segments[0].calc_finite_difference_effmass(),1.0)

def test_explosion_index(toy_segments):
 	toy_segments[0].dk_bohr = np.array([0.01,0.02,0.03,0.04,0.05,0.06,0.07]) # toy data with characteristic curve
 	toy_segments[0].dE_hartree = np.array([0,0.04,0.08,0.18,0.28,0.30,0.32])
 	assert toy_segments[0].explosion_index() == 4

def test_check_kanefit_points(toy_segments):
	with pytest.raises(AssertionError):
		toy_segments[0].check_kanefit_points(2)
	
def test_calc_alpha(toy_segments):
	toy_segments[0].dk_bohr = np.array([0.1,0.2,0.3,0.4,0.5]) # toy data with characteristic curve
	toy_segments[0].dE_hartree = np.array([0.00124534733,0.004927169016,0.01089396461,0.01892547876,0.02876732333])
	assert math.isclose(toy_segments[0].calc_alpha(polyfit_order=6,truncate=False),3,rel_tol=1E-1)

def test_calc_kane_mass_bandedge(toy_segments):
	toy_segments[0].dk_bohr = np.array([0.1,0.2,0.3,0.4,0.5]) # toy data with characteristic curve
	toy_segments[0].dE_hartree = np.array([0.00124534733,0.004927169016,0.01089396461,0.01892547876,0.02876732333])
	assert math.isclose(toy_segments[0].calc_kane_mass_bandedge(polyfit_order=6,truncate=False),4,1E-3)

# TODO: average_allband_optical_effmass, calc_allband_optical_effmass, functions associated with analytic optical Kane effmass, optical mass and use of weighting





