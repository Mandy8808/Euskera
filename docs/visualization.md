# Visualisation

Saved simulation data can be inspected with the plotting and saving helpers
in `euskera.visualization` and with the notebook
[`example_visualizing_simulation_data.ipynb`](../examples/simulations/example_visualizing_simulation_data.ipynb).
The visualisation workflow can plot slices, profiles, densities, potentials
and diagnostic histories; FFmpeg is only required for video export.

Simulation-specific figures are implemented in
`euskera.visualization.simulation_plots`. Reusable helpers such as
`colorBar_and_normaliz` and `colored_line` live in
`euskera.visualization.plot_tools`.

Background profile plots are part of the same canonical namespace:
`plotUsingPerf`, `plotUsingDiscSol`, `plotUsingSol` and `plotPerf` are
available from `euskera.visualization` (and from the package root). Shared
Matplotlib defaults and the `get_colors` palette live in
`euskera.visualization.configuration`.

Run simulations into a dedicated directory, then pass that directory to the
existing notebook or plotting helper. Representative figures and videos are
available under [`galleries/`](../galleries/), while sample outputs are in
[`simulation_data/`](../simulation_data/).

For frequency and oscillation analysis, see
[`example_computing_oscilation_soliton_frequency.ipynb`](../examples/simulations/example_computing_oscilation_soliton_frequency.ipynb).
The [API guide](api.md) lists diagnostic and frequency helpers.

## Output by variable, geometry, and time

See [selective output](selective-output.md) for `SaveRule`, `All`, `Last`, `TimeRange`, `Final`, independent diagnostics, all three planes and axes, and reading actual coordinates.
