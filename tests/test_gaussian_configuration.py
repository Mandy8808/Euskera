"""Gaussian configuration order and completeness regressions."""
from itertools import permutations
import numpy as np
import pytest
from euskera.models.models import Models
from euskera.models.gaussiana_model import Gaussiana_Model

CONFIG = {'positions_gaussiana': [0.2, -0.1, 0.3], 'amplitude': 2., 'sigma': [0.5, 0.8, 1.1]}
PARAMS = {'resol': 4, 'gridlength': 4., 'num_threads': 1}

@pytest.mark.parametrize('keys', list(permutations(CONFIG)))
def test_gaussian_order_independent(keys):
    model = Models(gaussian_function=[{key: CONFIG[key] for key in keys}])
    grid, (psi, rho) = model.call_model(1, PARAMS)
    expected = 2. * np.exp(-sum((x-c)**2/(2*s*s) for x,c,s in zip(grid[:3], CONFIG['positions_gaussiana'], CONFIG['sigma'])))
    np.testing.assert_allclose(psi[0], expected)
    np.testing.assert_allclose(rho[0], expected**2)

@pytest.mark.parametrize('key', list(CONFIG))
def test_missing_parameter_in_one_configuration_is_rejected(key):
    with pytest.raises(ValueError, match='missing'):
        Models(gaussian_function=[CONFIG, {k:v for k,v in CONFIG.items() if k != key}])

def test_direct_gaussian_lists_cannot_truncate():
    model = Gaussiana_Model({'positions_gaussiana': [[0,0,0]]*2, 'amplitude': [1.], 'sigma': [[1,1,1]]*2})
    with pytest.raises(ValueError, match='equal lengths'):
        model.apply(1, PARAMS)
