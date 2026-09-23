# Executable workflow notebooks

Install Euskera with its notebook dependencies:

```bash
python -m pip install -e ".[notebooks]"
```

Open a notebook using that environment and run every cell:

1. [Background profile](background_profile.ipynb): short radial integration,
   observables, saving and plotting. This is an initial-value demonstration,
   not a converged bound-state profile.
2. [Gaussian simulation](gaussian_simulation.ipynb): simulation, physical times,
   named HDF5 diagnostics, metadata and mass-conservation check.
3. [Generated-data visualization](visualize_generated_data.ipynb): generate
   NPZ density planes, read the last snapshot and save a plot.

Each notebook runs independently and writes under `example_output` by default.
Set `EUSKERA_EXAMPLE_OUTPUT` to choose another output directory. The examples
need no checked-in simulation data. Research notebooks remain in the other
example directories.

For fast automated execution of all three, from the repository root:

```bash
python scripts/run_notebook_smoke.py --output-dir notebook-results
```

The runner sets `EUSKERA_SMOKE_TEST=1`, executes all cells in fresh kernels,
and isolates each output directory. Only numerical sizes change in smoke mode.
Omit `--output-dir` to discard generated files after verification. The runner
uses the current Python interpreter and does not install a global Jupyter kernel.
