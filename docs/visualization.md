# Visualisation

Saved simulation data can be inspected with the plotting and saving helpers
in `euskera.visualization` and with the notebook
[`example_visualizing_simulation_data.ipynb`](../examples/example_visualizing_simulation_data.ipynb).
The visualisation workflow can plot slices, profiles, densities, potentials
and diagnostic histories; FFmpeg is only required for video export.

Run simulations into a dedicated directory, then pass that directory to the
existing notebook or plotting helper. Representative figures and videos are
available under [`galleries/`](../galleries/), while sample outputs are in
[`simulation_data/`](../simulation_data/).

For frequency and oscillation analysis, see
[`example_computing_oscilation_soliton_frequency.ipynb`](../examples/example_computing_oscilation_soliton_frequency.ipynb).
The [API guide](api.md) lists diagnostic and frequency helpers.
