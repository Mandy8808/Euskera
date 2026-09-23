"""Validate initial-condition dictionaries before allocation or output."""
import numpy as np


def validate_models(models, components):
    if isinstance(components, bool) or not isinstance(components, int) or components < 1:
        raise ValueError('field_components must be a positive integer')
    if not isinstance(models, dict) or not models:
        raise ValueError('model_parameters must be a nonempty dictionary')
    for name, configurations in models.items():
        if not isinstance(configurations, (list, tuple)) or not configurations:
            raise ValueError(f'{name} requires a nonempty configuration list')
        for index, config in enumerate(configurations):
            prefix = f'{name}[{index}]'
            if not isinstance(config, dict):
                raise ValueError(f'{prefix} must be a dictionary')
            sizes = ({'positions_gaussiana': 3, 'sigma': 3, 'amplitude': None}
                     if name == 'gaussian_function' else
                     {'positions': 3, 'velocities': 3, 'phases': 1, 'alphas': 1,
                      'dr': 1, 'betas': components})
            if name == 'ell_boson':
                sizes['ell'] = 1
            for key, size in sizes.items():
                if key not in config:
                    raise ValueError(f'{prefix}.{key} is required')
                value = np.asarray(config[key])
                if value.dtype.kind not in 'ifu' or not np.isfinite(value).all():
                    raise ValueError(f'{prefix}.{key} must contain finite real numbers')
                if value.shape != (() if size is None else (size,)):
                    raise ValueError(f'{prefix}.{key} has an invalid shape')
                if key in ('sigma', 'alphas', 'dr') and np.any(value <= 0):
                    raise ValueError(f'{prefix}.{key} must be positive')
            if name != 'gaussian_function':
                profiles = config.get('profiles', [])
                if len(profiles) != components:
                    raise ValueError(f'{prefix}.profiles needs one profile per component')
                for j, profile in enumerate(profiles):
                    array = np.asarray(profile)
                    if (array.ndim != 1 or array.size < 2 or array.dtype.kind not in 'ifu'
                            or not np.isfinite(array).all()):
                        raise ValueError(f'{prefix}.profiles[{j}] must be a finite real 1D array with at least two values')
            if name == 'proca':
                if components != 3 or config.get('polarization') not in (
                    ['radial'], ['circular'], ['linear_x'], ['linear_y'], ['linear_z'],
                    ('radial',), ('circular',), ('linear_x',), ('linear_y',), ('linear_z',)):
                    raise ValueError(f'{prefix}.polarization requires a supported polarization and three components')
            if name == 'ell_boson':
                ell = config['ell'][0]
                if ell < 0 or int(ell) != ell or components != 2 * ell + 1:
                    raise ValueError(f'{prefix}.ell requires a nonnegative integer and 2*ell+1 components')
