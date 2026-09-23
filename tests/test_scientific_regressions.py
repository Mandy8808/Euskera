"""Independent analytic checks for repaired numerical diagnostics."""
import numpy as np
import pytest
from scipy.linalg import block_diag, eigvals
from scipy.optimize import linear_sum_assignment
from euskera.core.grids import KGrid
from euskera.observables.simulation_conserv_quant import kintE, Pi, centpotetE, selfinterCondensateE
from euskera.spectral import espectro, backgroundOper
from euskera.spectral.blocks import linBlock

@pytest.fixture(params=['numpy', 'pyfftw'])
def fft_objects(request):
    if request.param == 'numpy':
        return (lambda a: np.fft.fftn(a, axes=(1,2,3)),
                lambda a: np.fft.ifftn(a, axes=(1,2,3)))
    fftw = pytest.importorskip('pyfftw')
    template = np.zeros((2,8,8,8), complex)
    return (fftw.builders.fftn(template, axes=(1,2,3)),
            fftw.builders.ifftn(template, axes=(1,2,3)))

@pytest.mark.parametrize('mode', [(1,0,0),(0,-2,0),(0,0,3),(1,-2,3),(4,0,0),(0,4,0),(0,0,4)])
def test_kinetic_energy_and_momentum_analytic(mode, fft_objects):
    x = np.arange(8) * (2*np.pi/8)
    coords = np.meshgrid(x,x,x,indexing='ij',sparse=True)
    wave = np.broadcast_to(np.exp(1j*sum(k*a for k,a in zip(mode,coords))), (8,8,8))
    psi = np.stack([wave, .5*wave])
    kv, k2 = KGrid(2*np.pi,8)
    cell = (2*np.pi/8)**3
    mass = (2*np.pi)**3 * 1.25
    for method in [1,2]:
        np.testing.assert_allclose(cell*kintE(psi,k2,fft_objects,kv,method),
                                   .5*mass*sum(k*k for k in mode),rtol=1e-12)
    if 4 not in mode:  # Nyquist samples do not determine the sign of momentum.
        np.testing.assert_allclose(Pi(psi,kv,cell,fft_objects),mass*np.array(mode),atol=1e-11)

@pytest.mark.parametrize('central_mass', [0., 2.])
def test_external_potential_has_full_weight(central_mass):
    r = np.array([1.,2.,3.]);rho=np.array([.2,.5,.8]);own=np.array([-.4,-.6,-.7])
    external = -central_mass/r
    params = {'cmass':central_mass}
    self_energy = selfinterCondensateE(rho,own+external,r,params)
    np.testing.assert_allclose(self_energy,.5*np.sum(rho*own))
    np.testing.assert_allclose(self_energy+centpotetE(rho,r,params),
                               .5*np.sum(rho*own)+np.sum(rho*external))

def test_linear_spectrum_contains_all_three_sectors():
    functions=([lambda r: np.exp(-r*r)],[lambda r: -.3*np.ones_like(r)])
    util=[1,6,4.]; couplings=(1.,.2)
    sigma,u,inv,_,d2,_,scale=backgroundOper(functions,util,couplings,0)
    m11,_,_,m22,m33=linBlock(6,sigma,u,inv,d2,couplings,scale)
    expected=1j*eigvals(block_diag(m11,m22,m33))
    actual=espectro(functions,util,couplings,0,'linear')[0]
    costs=abs(actual[:,None]-expected[None,:])
    rows,cols=linear_sum_assignment(costs)
    np.testing.assert_allclose(actual[rows],expected[cols],atol=1e-9,rtol=1e-10)
