"""Execute each workflow notebook in a fresh kernel and temporary workspace."""
import argparse
import json
import os
from pathlib import Path
import sys
import tempfile


def main():
    from nbclient import NotebookClient
    import nbformat
    from jupyter_client import KernelManager
    from jupyter_client.kernelspec import KernelSpecManager
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, help='Retain executed notebooks and generated data here')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix='euskera-notebooks-') as temporary:
        work = args.output_dir.resolve() if args.output_dir else Path(temporary)/'results'
        work.mkdir(parents=True, exist_ok=True)
        kernels = Path(temporary)/'kernels'
        spec = kernels/'euskera'
        spec.mkdir(parents=True)
        (spec/'kernel.json').write_text(json.dumps({
            'argv':[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'],
            'display_name':'Euskera smoke', 'language':'python'}))
        for path in sorted((root/'examples'/'workflows').glob('*.ipynb')):
            destination = work/path.stem
            destination.mkdir(exist_ok=True)
            environment = dict(os.environ, EUSKERA_SMOKE_TEST='1',
                EUSKERA_EXAMPLE_OUTPUT=str(destination), MPLBACKEND='Agg',
                MPLCONFIGDIR=str(Path(temporary)/'matplotlib'),
                IPYTHONDIR=str(Path(temporary)/'ipython'),
                PYTHONPATH=str(root)+os.pathsep+os.environ.get('PYTHONPATH',''))
            manager = KernelManager(kernel_name='euskera',
                kernel_spec_manager=KernelSpecManager(kernel_dirs=[str(kernels)]),
                connection_file=str(Path(temporary)/(path.stem+'.json')))
            notebook = nbformat.read(path, as_version=4)
            client = NotebookClient(notebook, km=manager, timeout=120, allow_errors=False,
                                    resources={'metadata':{'path':str(destination)}})
            try:
                client.execute(env=environment)
            finally:
                if manager.has_kernel:
                    manager.shutdown_kernel(now=True)
            nbformat.write(notebook, destination/path.name)
            print(f'PASS {path.name}', flush=True)


if __name__ == '__main__':
    main()
