# Selective simulation output

Output rules select **variables, geometry, orientation, and time** before evolution begins. The simulation remains three-dimensional: saving a plane does not change the equations or the integration grid.

## Example: complete plane histories and the last 100 volumes

```python
from euskera import OutputConfig, SaveRule, All, Last, Final, TimeRange

output = OutputConfig(
    address="Data/selective_run",
    format="hdf5",  # "npz" is also supported.
    save_number=1000,
    rules=[
        SaveRule("plane", ["rho", "psi", "phi"], ["xy", "xz", "yz"]),
        SaveRule("line", ["rho", "psi", "phi"], ["x", "y", "z"]),
        SaveRule("volume", ["rho", "psi", "phi"], selection=Last(100)),
        SaveRule("diagnostics", ["Numb_Part", "Energ"], every=10),
    ],
)
# evolve(model_parameters, output_config=output, ...)
```

`rules=None` preserves the existing `data_save` behavior. An explicit rule list replaces `data_save`, including its defaults; `rules=[]` saves only the grid and metadata. Variables in one rule share its selection. Use separate rules for different time selections. Overlapping rules for the same output are combined: each sample is written only once.

Requested diagnostics must also be enabled in `DiagnosticsConfig`. Particle number and energy are enabled by default; `Pi` and `Frequency` are not. `Ji` is not implemented yet.

## Calendar and selections

`save_number=N` defines N intervals and N+1 candidate samples, numbered 0 through N. A sample does not necessarily correspond to one integration step. Global indices and physical times are preserved even when an output starts late.

| Selection | Result |
| --- | --- |
| `All()` (default) | All eligible samples, including the initial sample |
| `Last(100)` | The last 100 eligible samples, including the final sample |
| `Final()` | Only the final sample |
| `TimeRange(start, stop)` | Samples in the closed physical-time interval [start, stop] |

`every=k` first selects indices 0, k, 2k, ... and adds N if it is not a multiple of k. The temporal selection is then applied. Thus, `Last(100), every=2` keeps the last 100 samples of that reduced calendar. If more samples are requested than are available, all eligible samples are saved and a warning is emitted. An interval containing no samples produces an empty output readable with `read_output`.

`save_initial=True` belongs to **each SaveRule**: it adds sample 0 even if it falls outside the temporal selection. It never removes samples: combined with `Last(100)`, it can produce 101 samples. Its default value, `False`, means that zero is not added outside the selection; `All()` still includes zero. The grid and metadata are saved even when no rule selects zero.

With `save_number=1000`, `Last(100)`, and `every=1`, the indices are 901–1000. To save only this interval for every field output, apply `Last(100)` to the line, plane, and volume rules. Diagnostics can still use `All()`.

The time step remains fixed. The existing `dtime` adjustment is preserved: `save_number` contributes to the total integration-step count. Changing only the rules does not change the time step; changing `save_number` can. Adaptive termination with rolling retention of the last frames is not supported here.

## Orientations, variables, and coordinates

| Geometry | Orientation | Slice |
| --- | --- | --- |
| `plane` | `xy` | z nearest to 0 |
| `plane` | `xz` | y nearest to 0 |
| `plane` | `yz` | x nearest to 0 |
| `line` | `x` | y and z nearest to 0 |
| `line` | `y` | x and z nearest to 0 |
| `line` | `z` | x and y nearest to 0 |
| `volume` | No orientations | Full volume |

`rho`, `psi`, and `phi` support every geometry. All complex components of `psi` are retained. Each fixed coordinate uses the grid node nearest to zero; ties use the first index. `output_manifest.json` records the actual indices and positions of the slices. Spatial axes retain x, y, z order: for example, a `yz` plane has y, z dimensions, preceded by the component dimension for `psi`.

## Diagnostics and the potential half-step

Keeping particle number and energy throughout evolution is recommended. Their series require little storage, although energy and momentum can be expensive to compute. `every` reduces that cost. `Frequency` saves complex field samples for subsequent analysis; it is not a conserved quantity. For spectral analysis, choose a uniform calendar: if `every` does not divide N, the last interval can be shorter.

Each diagnostic has its own output: `save_energies_Numb_Part`, `save_energies_Energ`, `save_energies_Pi`, and `save_energies_Frequency`. They can share a rule or use independent calendars. `Conserv` is called once per sampling instant for the union of requested diagnostics.

When no sample is needed, the method combines the final potential half-step with the initial half-step of the next iteration. If any rule requests fields or diagnostics, it completes the pending half-step **once**, obtains all requested outputs, and starts the next iteration with a half-step. This uniform closure also applies when only `rho` or `phi` is requested. Any pending half-step is completed at termination, even when the selection excludes the endpoint. Numerical comparisons between selections allow for rounding differences.

Enabled initial diagnostics are still computed to validate the initial state, even when they are not written. All evolution FFTs and potential updates remain necessary before a final output window; the savings come from omitted output extraction, writes, and diagnostic calculations.

## Reading and visualization

```python
from euskera import read_output
from euskera.visualization import plot_output

volume = read_output("Data/selective_run", "volume_psi")
print(volume["snapshot_index"], volume["time"])
# volume["data"]: (samples, components, nx, ny, nz)

plane = read_output("Data/selective_run", "plane_yz_phi")
print(plane["coordinates"], plane["fixed_coordinates"])
plot_output("Data/selective_run", "plane_yz_phi", sample=-1)
plot_output("Data/selective_run", "line_z_psi", sample=-1, component=0)

energy = read_output("Data/selective_run", "save_energies_Energ")
print(energy["time"], energy["data"]["energy"])
```

`read_output` supports NPZ and HDF5 and loads one complete output into memory. For large volumes, access individual samples directly in HDF5. Fields are stacked along the sample axis; diagnostics are returned as a dictionary of numeric arrays. NPZ retains the existing object encoding for diagnostics: only open files from trusted sources. HDF5 uses numeric datasets.

`plot_output` draws lines and planes using their actual coordinates and times. `sample` is a position within the output (−1 means the last sample), not a global index. For `psi`, it displays the modulus of the selected component. Existing files continue to use the existing readers. The existing animation now intersects output indices and uses physical times when available; for new slices, `plot_output` avoids manual axis specification.

The formats still use per-sample files during execution and `end_...` archives at closure. Use a new directory for each run; resuming into existing results is not implemented.

## Backward compatibility

The existing `salva_data_update` argument and `OutputConfig(data_save=...)` remain supported. `rules` defaults to `None`, so old calls use `data_Objgenerator` and `fdata_save`, retaining existing filenames, the initial frame, and the common saving cadence. Legacy planes remain central xy slices and legacy lines remain central x slices.

When `rules` is an explicit list, `OutputSchedule` handles output instead. `data_save` remains a validated configuration field but does not choose outputs in that mode. The two output mechanisms do not both write results. Passing `rules=[]` is intentionally different from omitting rules: it requests no field or diagnostic outputs.

If both `salva_data_update` and `output_config` are passed, the latter takes precedence, as before. The notebook switch `USE_SELECTIVE_OUTPUT` is only an example convenience; it is not an argument to `evolve`. Setting it to `False` passes `output_config=None` and preserves the existing notebook configuration.

New selected outputs use explicit names such as `plane_yz_phi` and require `read_output` or direct access to their files; backward compatibility does not mean that existing analysis code automatically understands the new filenames.

## Examples and verification

The [selective_output.ipynb](../examples/workflows/selective_output.ipynb) notebook executes a small case, compares a plane with its volume slice, reads energy, and plots outputs. Simulation notebooks include an optional section before evolution; they initially retain the old mode so their existing analysis cells continue to work.

Tests cover both formats, all orientations and variables, complex components, independent times, empty and overlapping selections, and numerical equivalence between complete and partial output.

## Recoverable consolidation

`cleanup_policy="after_success"` and `consolidation_batch_size=16` apply to every writer, including the grid, selected fields, and diagnostics. They are independent of `SaveRule` temporal selections. See [consolidation and recovery](simulations.md#consolidation-and-disk-space) for disk requirements, incremental cleanup, and retry instructions. For a selected output, use its actual stream name when retrying, for example `StoreSolution(address, "volume_psi", format="hdf5").close_file("end_volume_psi")`. Match the original cleanup policy. Do not rerun evolution into an interrupted output directory.
