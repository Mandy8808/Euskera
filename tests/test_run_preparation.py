"""Fail before output, and persist trustworthy run coordinates and inputs."""
import json
import h5py
import numpy as np
import pytest
import euskera
from euskera.io import OutputConfig
from euskera.io.provenance import write_run_metadata

MODEL={'gaussian_function':[{'positions_gaussiana':[0.,0.,0.],'sigma':[1.,1.,1.],'amplitude':1.}]}

@pytest.mark.parametrize('typed',[False,True])
@pytest.mark.parametrize('values',[{'gridlength':float('nan')},{'t0':1.},{'methodEnerg':3},{'resol':0},{'cmass':float('inf')}])
def test_invalid_numerics_do_not_create_output(tmp_path,typed,values):
    target=tmp_path/'absent'
    kwargs={'evolution_config' if typed else 'simulation_parameters_update':values}
    with pytest.raises((ValueError,TypeError)):
        euskera.evolve(MODEL,salva_data_update={'address':str(target)},**kwargs)
    assert not target.exists()

@pytest.mark.parametrize('output',[{'format':'csv'},{'data_save':{'save_typo':True}}])
def test_invalid_output_has_no_files(tmp_path,output):
    target=tmp_path/'absent'
    with pytest.raises(ValueError):
        euskera.evolve(MODEL,salva_data_update={'address':str(target),**output})
    assert not target.exists()

@pytest.mark.parametrize('sigma',[[1.,0.,1.],[1.,-1.,1.],[1.,float('nan'),1.]])
def test_invalid_model_before_output(tmp_path,sigma):
    target=tmp_path/'absent';model={'gaussian_function':[{**MODEL['gaussian_function'][0],'sigma':sigma}]}
    with pytest.raises(ValueError,match='sigma'):
        euskera.evolve(model,salva_data_update={'address':str(target)})
    assert not target.exists()

@pytest.mark.parametrize('fmt',['npz','hdf5'])
def test_saved_physical_times_and_metadata(tmp_path,fmt):
    euskera.evolve(MODEL,evolution_config={'resol':4,'gridlength':4.,'tmax':.2},
                  output_config=OutputConfig(address=str(tmp_path),format=fmt,save_number=2))
    path=tmp_path/f'end_save_line_rho.{fmt}'
    archive=np.load(path) if fmt=='npz' else h5py.File(path)
    with archive as data:
        np.testing.assert_allclose(data['time'][()], [0.,.1,.2])
        np.testing.assert_array_equal(data['snapshot_index'][()], [0,1,2])
        np.testing.assert_array_equal(data['t'][()], [0,1,2])
    metadata=json.loads((tmp_path/'run_metadata.json').read_text())
    assert metadata['run_id'] and metadata['python']
    assert metadata['fft_backend'] in ('numpy','pyfftw')
    assert metadata['steps']*metadata['timestep']==pytest.approx(.2)
    assert metadata['evolution']['resol']==4


def test_profile_hashes_and_optional_copies(tmp_path):
    model={'soliton':[{'profiles':[np.array([1.,.5,0.])],'alphas':[2.]}]}
    for copy in [False,True]:
        write_run_metadata(str(tmp_path),model,{'lambda_value':1.},
                           {'copy_profiles':copy},{},1,'numpy',.1,2)
        metadata=json.loads((tmp_path/'run_metadata.json').read_text())
        profile=metadata['models']['soliton'][0]['profiles'][0]
        assert len(profile['sha256'])==64
        assert metadata['models']['soliton'][0]['alphas']==[1.]
        if copy:
            with np.load(tmp_path/'input_profiles.npz') as saved:
                np.testing.assert_array_equal(saved[profile['copy_key']],model['soliton'][0]['profiles'][0])
        else:
            assert not (tmp_path/'input_profiles.npz').exists()
