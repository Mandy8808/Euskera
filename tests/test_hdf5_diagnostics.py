"""Named diagnostics round trips and full evolution output."""
from itertools import product
import h5py
import numpy as np
import pytest
import euskera
from euskera.core.grids import RealGrid, KGrid
from euskera.io import read_hdf5_diagnostics
from euskera.io.save_data import data_Objgenerator, StoreSolution
from euskera.observables.simulation_conserv_quant import Conserv

@pytest.mark.parametrize('enabled', list(product([False, True], repeat=4)))
def test_diagnostic_subsets_and_complex_samples_roundtrip(tmp_path, enabled):
    # Non-default order makes mislabeling positional diagnostics visible.
    options = dict(zip(['Frequency', 'Pi', 'Energ', 'Numb_Part'], enabled))
    xyz, radius = RealGrid(2*np.pi, 4)
    kv, k2 = KGrid(2*np.pi, 4)
    psi = np.stack([np.broadcast_to(a*np.exp(1j*xyz[0]), (4,4,4)) for a in (1., .5)])
    rho = np.sum(abs(psi)**2, axis=0)
    params = {'resol':4, 'gridlength':2*np.pi, 'num_threads':1, 'cmass':0., 'lambda_value':0.}
    values = Conserv([psi,rho,np.zeros_like(rho),radius,k2,kv],options,params,max_pos=[(0,0,0)]*2)
    obj = data_Objgenerator({'save_energies':True},str(tmp_path),'hdf5',options)['save_energies']
    for index in [0, 2, 10]:
        obj.save_file(values,index)
    obj.close_file('end_save_energies')
    path = tmp_path/'end_save_energies.hdf5'
    result = read_hdf5_diagnostics(path)
    assert list(result) == [0,2,10]
    expected = {}
    for key,value in zip([key for key,on in options.items() if on],values):
        if key == 'Numb_Part':
            expected.update(mass_total=value[0],mass_components=np.asarray(value[1],float))
        else:
            expected[{'Frequency':'frequency_samples','Pi':'momentum','Energ':'energy'}[key]] = np.asarray(value,dtype=complex if key=='Frequency' else float)
    for snapshot in result.values():
        assert snapshot.keys() == expected.keys()
        for key,value in expected.items():
            np.testing.assert_allclose(snapshot[key],value)
    with h5py.File(path) as archive:
        assert list(archive['save_energies_0'].attrs['diagnostic_names']) == [key for key,on in options.items() if on]
        for value in archive['save_energies_0'].values():
            assert value.dtype.kind in 'fc'


def test_hdf5_full_evolution_default_diagnostics(tmp_path):
    euskera.evolve(
        {'gaussian_function':[{'amplitude':1., 'sigma':[1.,1.,1.], 'positions_gaussiana':[0.,0.,0.]}]},
        evolution_config=euskera.EvolutionConfig(gridlength=4.,resol=4,tmax=.01),
        output_config=euskera.OutputConfig(address=str(tmp_path),format='hdf5',save_number=1),
    )
    data = read_hdf5_diagnostics(tmp_path/'end_save_energies.hdf5')
    assert list(data) == [0,1]
    for snapshot in data.values():
        assert set(snapshot) == {'energy','mass_total','mass_components'}
        assert np.isfinite(snapshot['energy'])
    np.testing.assert_allclose(data[0]['mass_total'],data[1]['mass_total'],rtol=1e-12)


def test_unlabelled_diagnostics_rejected_before_file_creation(tmp_path):
    obj = StoreSolution(str(tmp_path),'save_energies','hdf5')
    with pytest.raises(ValueError, match='diagnostic_names'):
        obj.save_file(np.array([1.],dtype=object),0)
    assert not list(tmp_path.iterdir())
