#! /usr/bin/env python3
import pytest
from effmass import outputs
import matplotlib

@pytest.mark.parametrize("data_object,settings_object,segment_object", [
    (pytest.lazy_fixture('MAPI_soc_data_object_with_DOSCAR'), pytest.lazy_fixture("MAPI_settings_object"),pytest.lazy_fixture("MAPI_soc_segment_object_hole")),
])
def test_plot_segments(data_object,settings_object,segment_object):
    fig, ax = outputs.plot_segments(data_object,settings_object,[segment_object])
    assert type(fig) == matplotlib.figure.Figure

@pytest.mark.parametrize("data_object", [
    (pytest.lazy_fixture('MAPI_soc_data_object_with_DOSCAR')),
])
def test_plot_integrated_dos(data_object):
    fig, ax = outputs.plot_integrated_dos(data_object)
    assert type(fig) == matplotlib.figure.Figure

@pytest.mark.parametrize("data_object", [
    (pytest.lazy_fixture('MAPI_soc_data_object_with_DOSCAR')),
])
def test_plot_dos(data_object):
    fig, ax = outputs.plot_dos(data_object)
    assert type(fig) == matplotlib.figure.Figure

@pytest.mark.parametrize("segment_object,poly_degree", [
    (pytest.lazy_fixture('MAPI_soc_segment_object_hole'),0),
])
def test_print_results_invalid_degree(segment_object,poly_degree):   
    with pytest.raises(AssertionError):
        outputs.print_results(segment_object,poly_degree)

