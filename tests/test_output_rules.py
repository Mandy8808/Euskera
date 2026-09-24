"""Spatial and temporal selection must preserve physical samples."""
import json
import numpy as np
import pytest
import euskera
from euskera import All, Last, Final, TimeRange, SaveRule, OutputConfig, read_output
from euskera.io.selected_output import OutputSchedule


@pytest.mark.parametrize('fmt', ['npz', 'hdf5'])
def test_spatial_cuts_and_independent_times(tmp_path, fmt):
    axes = [np.array([-2., -.1, 1.]), np.array([-1., .2, 3.]), np.array([-3., -.2, 2.])]
    fields = np.arange(27).reshape(3, 3, 3)
    psi = np.stack([fields + 1j, 2*fields - 1j])
    rules = [SaveRule('plane', ['rho', 'psi', 'phi'], ['xy', 'xz', 'yz']),
             SaveRule('line', ['rho', 'psi', 'phi'], ['x', 'y', 'z']),
             SaveRule('volume', ['rho', 'psi', 'phi'], selection=Last(2))]
    schedule = OutputSchedule(OutputConfig(address=str(tmp_path), format=fmt, save_number=4, rules=rules), 1., {})
    schedule.prepare(axes)
    for j in range(5):
        schedule.save(schedule.outputs_at(j, j/4), j, j/4, fields+j, psi+j, -fields+j, {})
    schedule.close()
    for variable, field in [('rho', fields), ('psi', psi), ('phi', -fields)]:
        for geometry, orientations in [('plane', ['xy','xz','yz']), ('line', ['x','y','z']), ('volume', [None])]:
            for orientation in orientations:
                name = '_'.join(v for v in [geometry, orientation, variable] if v)
                result = read_output(tmp_path, name)
                expected_indices = [3,4] if geometry == 'volume' else list(range(5))
                np.testing.assert_array_equal(result['snapshot_index'], expected_indices)
                np.testing.assert_allclose(result['time'], np.array(expected_indices)/4)
                varying = orientation or 'xyz'
                slices = (Ellipsis,) + tuple(slice(None) if a in varying else 1 for a in 'xyz')
                for index, value in zip(expected_indices, result['data']):
                    np.testing.assert_array_equal(value, (field+index)[slices])
                assert result['axes'] == list(varying)


def test_selection_roundtrip_and_validation():
    for selection, expected in [(All(), [0,1,2,3,4]), (Last(2), [3,4]), (Final(), [4]), (TimeRange(.25,.75), [1,2,3])]:
        config = OutputConfig(rules=[SaveRule('volume', ['psi'], selection=selection)])
        rule = OutputConfig.from_legacy(config.to_dict()).rules[0]
        assert [i for i in range(5) if rule.matches(i,4,i/4,1)] == expected
    rule = SaveRule('volume', ['rho'], selection=Last(2), every=3, save_initial=True)
    assert [i for i in range(11) if rule.matches(i,10,i/10,1)] == [0,9,10]
    for kwargs in [dict(geometry='plane', variables=['rho']), dict(geometry='line',variables=['phi'],orientations=['xy']), dict(geometry='volume',variables=['bad'])]:
        with pytest.raises(ValueError): SaveRule(**kwargs)
    with pytest.raises(ValueError): Last(0)
    with pytest.raises(ValueError): TimeRange(2,1)
    with pytest.raises(ValueError): OutputSchedule(OutputConfig(rules=[SaveRule('diagnostics',['Energ'])]),1,{'Energ':False})


def run_model(path, fmt, rules):
    euskera.evolve({'gaussian_function':[{'positions_gaussiana':[.2,0.,0.], 'amplitude':1., 'sigma':[.5,.7,.9]}]},
        evolution_config=euskera.EvolutionConfig(gridlength=4.,resol=4,tmax=.04,lambda_value=.2),
        output_config=OutputConfig(address=str(path),format=fmt,save_number=4,rules=rules),
        diagnostics_config=euskera.DiagnosticsConfig(Numb_Part=True,Energ=True,Pi=True,Frequency=True))


@pytest.mark.parametrize('fmt', ['npz','hdf5'])
def test_sparse_evolution_matches_full_and_diagnostics(tmp_path, fmt):
    all_rules = [SaveRule('volume',['rho','psi','phi']), SaveRule('diagnostics',['Numb_Part','Energ','Pi','Frequency'])]
    run_model(tmp_path/'all',fmt,all_rules)
    sparse_rules = [SaveRule('volume',['rho','psi','phi'],selection=Last(2)),
                    SaveRule('diagnostics',['Numb_Part','Energ','Pi','Frequency'],selection=Last(2))]
    run_model(tmp_path/'sparse',fmt,sparse_rules)
    for variable in ['rho','psi','phi']:
        full = read_output(tmp_path/'all',f'volume_{variable}')
        sparse = read_output(tmp_path/'sparse',f'volume_{variable}')
        np.testing.assert_array_equal(sparse['snapshot_index'], [3,4])
        np.testing.assert_allclose(sparse['data'],full['data'][-2:],rtol=1e-12,atol=1e-12)
    for diagnostic in ['Numb_Part','Energ','Pi','Frequency']:
        full = read_output(tmp_path/'all',f'save_energies_{diagnostic}')
        sparse = read_output(tmp_path/'sparse',f'save_energies_{diagnostic}')
        for name in full['data']:
            np.testing.assert_allclose(sparse['data'][name],full['data'][name][-2:],rtol=1e-11,atol=1e-11)
    # Diagnostic-only sampling between sparse field outputs also keeps the same field.
    run_model(tmp_path/'mixed',fmt,[SaveRule('volume',['psi'],selection=Last(2)), all_rules[1]])
    np.testing.assert_allclose(read_output(tmp_path/'mixed','volume_psi')['data'],read_output(tmp_path/'all','volume_psi')['data'][-2:],atol=1e-12)


def test_empty_window_and_overlapping_rules(tmp_path):
    rules=[SaveRule('volume',['rho'],selection=TimeRange(2,3)),
           SaveRule('line',['rho'],['z']), SaveRule('line',['rho'],['z'],selection=Final())]
    run_model(tmp_path,'npz',rules)
    assert read_output(tmp_path,'volume_rho')['data'].size == 0
    assert len(read_output(tmp_path,'line_z_rho')['time']) == 5


def test_plot_uses_orientation_and_physical_time(tmp_path):
    import matplotlib.pyplot as plt
    from euskera.visualization import plot_output
    run_model(tmp_path, 'npz', [SaveRule('plane',['phi'],['yz'],selection=Final()),
                              SaveRule('line',['psi'],['z'],selection=Final())])
    ax = plot_output(tmp_path, 'plane_yz_phi')
    assert ax.get_xlabel() == 'y' and ax.get_ylabel() == 'z'
    assert 'sample=4' in ax.get_title()
    ax2 = plot_output(tmp_path, 'line_z_psi')
    assert ax2.get_xlabel() == 'z'
    plt.close('all')


def test_animation_aligns_sparse_streams(monkeypatch):
    from euskera.visualization.video_make import Visualization
    from euskera.visualization import video_make
    view = object.__new__(Visualization)
    view.fig = object()
    captured = {}
    def draw(dataS, dataL, *args):
        captured['line'] = dataL
        captured['plane'] = dataS
        return None, None, None
    view.imagshow02D = draw
    monkeypatch.setattr(video_make.animation, 'FuncAnimation', lambda *a, **kw: kw)
    line = {'snapshot_index':np.array([0,2,4]), 'time':np.array([0.,.2,.4]),
            'line0':np.zeros(3), 'line2':np.ones(3), 'line4':np.ones(3)*2}
    plane = {'snapshot_index':np.array([3,4]), 'time':np.array([.3,.4]),
             'plane3':np.zeros((3,3)), 'plane4':np.ones((3,3))}
    result = view.fplot2D(0,[np.arange(3),np.arange(3)],line,'line','density',(-1,1),(-1,1),None,50,plane,'plane',999)
    assert list(result['frames']) == [0]
    assert captured['line'][2] == captured['plane'][1] == .4
    np.testing.assert_array_equal(captured['line'][0], line['line4'])
