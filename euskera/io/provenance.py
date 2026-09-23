"""Portable run metadata and optional copies of input radial profiles."""
from datetime import datetime, timezone
from hashlib import sha256
from importlib.metadata import version, PackageNotFoundError
import json
from pathlib import Path
import platform
import subprocess
from uuid import uuid4
import numpy as np


def write_run_metadata(address, models, evolution, output, diagnostics, components,
                       backend, timestep, steps):
    """Record effective settings; replace profile arrays by content descriptors."""
    versions = {}
    for package in ('euskera', 'numpy', 'scipy', 'numba', 'numexpr', 'h5py', 'pyfftw'):
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = None
    root = Path(__file__).resolve().parents[2]
    git = {'revision': None, 'dirty': None}
    try:
        git['revision'] = subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=root, stderr=subprocess.DEVNULL,
            timeout=5, text=True).strip()
        status = subprocess.check_output(
            ['git', '--no-optional-locks', 'status', '--porcelain'], cwd=root,
            stderr=subprocess.DEVNULL, timeout=5, text=True)
        git['dirty'] = bool(status.strip())
    except (OSError, subprocess.SubprocessError):
        pass
    copies = {}
    descriptions = {}
    for name, configs in models.items():
        descriptions[name] = []
        for i, config in enumerate(configs):
            item = dict(config)
            if 'profiles' in item:
                profiles = []
                for j, values in enumerate(item['profiles']):
                    array = np.ascontiguousarray(values)
                    digest = sha256(array.dtype.str.encode() + str(array.shape).encode() + array.tobytes()).hexdigest()
                    key = f'{name}_{i}_{j}'
                    profiles.append({'sha256': digest, 'dtype': array.dtype.str,
                                     'shape': list(array.shape),
                                     'copy_key': key if output['copy_profiles'] else None})
                    if output['copy_profiles']:
                        copies[key] = array
                item['profiles'] = profiles
            # Model builders enforce alpha=1 when self-interaction is enabled.
            if evolution['lambda_value'] != 0 and 'alphas' in item:
                item['alphas'] = [1.0]
            descriptions[name].append(item)
    path = Path(address)
    if copies:
        np.savez(path / 'input_profiles.npz', **copies)
    metadata = {
        'schema_version': 1, 'run_id': str(uuid4()),
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'python': platform.python_version(), 'versions': versions, 'git': git,
        'fft_backend': backend, 'field_components': components,
        'evolution': evolution, 'output': output, 'diagnostics': diagnostics,
        'models': descriptions, 'timestep': timestep, 'steps': steps,
        'time_convention': 't0=0; tmax is duration; time is measured in simulation code units',
    }
    def convert(value):
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, np.ndarray):
            return value.tolist()
        raise TypeError(f'Cannot serialize {type(value)}')
    (path / 'run_metadata.json').write_text(json.dumps(metadata, indent=2, default=convert, allow_nan=False)+'\n')
