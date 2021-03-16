#! /usr/bin/env python3
import pytest
from effmass import inputs, extrema
from vasppy import procar, outcar
import numpy as np
import math
import os

@pytest.fixture()
def MAPI_settings_object():
	return inputs.Settings(extrema_search_depth=0.025, energy_range=0.25)

@pytest.fixture()
def Si_nsp_data_object(): # Si nsp calculation data
        return inputs.DataCastep(os.path.join(os.path.dirname(__file__), 'data_castep/Si_nsp/'), 'Si')

@pytest.fixture()
def GaAs_sp_data_object(): # GaAs sp calculation data
        return inputs.DataCastep(os.path.join(os.path.dirname(__file__), 'data_castep/GaAs_sp/'), 'GaAs')

@pytest.fixture()
def MAPI_soc_data_object(): # MAPI SoC calculation data
	return inputs.DataVasp(os.path.join(os.path.dirname(__file__), 'data_vasp/MAPI_soc_OUTCAR'), os.path.join(os.path.dirname(__file__), 'data_vasp/MAPI_soc_PROCAR'), ignore=216)

@pytest.fixture()
def Ge_SP_data_object(): # Ge spin-polarised calculation from Adam Jackson
 	return inputs.DataVasp(os.path.join(os.path.dirname(__file__), 'data_vasp/Ge_SP_OUTCAR'), os.path.join(os.path.dirname(__file__), 'data_vasp/Ge_SP_PROCAR'))

@pytest.fixture()
def MAPI_cl_data_object(): # MAPI cl calculation data
	return inputs.DataVasp(os.path.join(os.path.dirname(__file__), 'data_vasp/MAPI_cl_OUTCAR'), os.path.join(os.path.dirname(__file__), 'data_vasp/MAPI_cl_PROCAR'), ignore=216)

@pytest.fixture()
def Ge_nsp_aims_data_object(): # Ge non spin-polarised calculation done with FHI-Aims
        return inputs.DataAims(os.path.join(os.path.dirname(__file__), 'data_aims/Ge_nsp_aims'))

@pytest.fixture()
def Ge_sp_aims_data_object(): # Ge spin-polarised calculation done with FHI-Aims
        return inputs.DataAims(os.path.join(os.path.dirname(__file__), 'data_aims/Ge_sp_aims'))

@pytest.fixture()
def Ge_soc_aims_data_object(): # Ge spin-orbit coupling calculation done with FHI-Aims
        return inputs.DataAims(os.path.join(os.path.dirname(__file__), 'data_aims/Ge_soc_aims'))

@pytest.fixture()
def MAPI_soc_data_object_with_DOSCAR(MAPI_settings_object): # MAPI spin orbit coupling calculation data with doscar loaded
    data = inputs.DataVasp(os.path.join(os.path.dirname(__file__), 'data_vasp/MAPI_soc_OUTCAR'), os.path.join(os.path.dirname(__file__), 'data_vasp/MAPI_soc_PROCAR'), ignore=216)
    data.parse_DOSCAR(os.path.join(os.path.dirname(__file__), 'data_vasp/MAPI_soc_DOSCAR'))
    return data


@pytest.fixture()
def MAPI_soc_segment_object_valence_band(MAPI_soc_data_object_with_DOSCAR,MAPI_settings_object):
    segments = extrema.generate_segments(MAPI_settings_object, MAPI_soc_data_object_with_DOSCAR, truncate_dir_change=False)
    return segments[3]

@pytest.fixture()
def MAPI_soc_segment_object_conduction_band(MAPI_soc_data_object_with_DOSCAR,MAPI_settings_object):
    segments = extrema.generate_segments(MAPI_settings_object, MAPI_soc_data_object_with_DOSCAR, truncate_dir_change=False)
    return segments[-1]

@pytest.fixture()
def toy_settings_object():
	return inputs.Settings(extrema_search_depth=0.02588716 , energy_range=0.02588716*3)

@pytest.fixture()
def toy_settings_object_smaller_range():
	return inputs.Settings(extrema_search_depth=0.02588716 , energy_range=0.02588716*2 )

@pytest.fixture()
def toy_data_object():
    data_object = inputs.DataVasp(os.path.join(os.path.dirname(__file__), 'data_vasp/MAPI_soc_OUTCAR'), os.path.join(os.path.dirname(__file__), 'data_vasp/MAPI_soc_PROCAR'),ignore=0)
    
    data_object.number_of_kpoints = 5
    data_object.number_of_bands = 2
    data_object.number_of_ions = 1
    data_object.kpoints = np.array([[0,0,0],[0.25,0,0],[0.5,0,0],[0.5,0.25,0],[0.5,0.5,0]])
    data_object.energies = np.array([[1,0.5],[2,1],[1,0.49124914462],[2,1.04],[1,0.4649965785],[2,1.07],[1,0.48],[2,1.04],[1,0.46],[2,1.0005]])[:,1:].reshape(data_object.number_of_kpoints,data_object.number_of_bands).T
    data_object.occupancy = np.array([[1,1],[2,0],[1,1],[2,0],[1,1],[2,0],[1,1],[2,0],[1,1],[2,0]])[:,1:].reshape(data_object.number_of_kpoints,data_object.number_of_bands).T   
    data_object.reciprocal_lattice = np.array([[1,0,0],[0,1,0],[0,0,1]])
    data_object.CBM = 1
    data_object.VBM = 0.5
    data_object.fermi_energy = 0.75
    data_object.dos = np.array([[0.42, 0],[0.46,7],[0.47,3],[0.48,5],[0.49,3],[0.5,24],[0.6,0],[0.7,0],[0.8,0],[0.9,0],[1.0,3.1],[1.02,0.2],[1.04, 3.20],[1.06,1.215]])
    data_object.integrated_dos = np.array([[0.42,3],[0.46,3.02],[0.47,3.04],[0.48,3.08],[0.49,3.08],[0.5,3.16],[0.6,3.16],[0.7,3.16],[0.8,3.16],[0.9,3.16],[1.0,3.17],[1.02,3.19],[1.04, 3.20],[1.06,3.215]])


    return data_object

@pytest.fixture()
def toy_segments(toy_data_object,toy_settings_object):
	segments = extrema.generate_segments(toy_settings_object, toy_data_object, truncate_dir_change=True)
	return segments
