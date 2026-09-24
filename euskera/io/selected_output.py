"""Write selected spatial outputs and read them with their own coordinates."""

import json
from pathlib import Path
import warnings
import numpy as np
import h5py

from .schedule import Last
from .save_data import StoreSolution, _diagnostic_datasets

__all__ = ["read_output"]


class OutputSchedule:
    """Resolve rules to unique streams; overlapping rules form a union."""

    def __init__(self, config, duration, diagnostics):
        self.config = config
        self.duration = duration
        self.total = config.save_number
        self.streams = {}
        for rule in config.rules:
            if isinstance(rule.selection, Last):
                available = self.total // rule.every + 1 + bool(self.total % rule.every)
                if rule.selection.count > available:
                    warnings.warn(f"Last({rule.selection.count}) exceeds {available} available samples; saving all eligible samples", UserWarning)
            for variable in rule.variables:
                if rule.geometry == "diagnostics" and not diagnostics.get(variable):
                    raise ValueError(f"Diagnostic {variable} must also be enabled in DiagnosticsConfig")
                for orientation in rule.orientations or [None]:
                    name = (f"save_energies_{variable}" if rule.geometry == "diagnostics" else
                            "_".join(part for part in (rule.geometry, orientation, variable) if part))
                    stream = self.streams.setdefault(name, {
                        "geometry": rule.geometry, "variable": variable,
                        "orientation": orientation, "rules": [],
                    })
                    stream["rules"].append(rule)

    def prepare(self, grid):
        self.axes = {axis: np.asarray(values).reshape(-1) for axis, values in zip("xyz", grid)}
        self.center = {axis: int(np.argmin(abs(values))) for axis, values in self.axes.items()}
        path = Path(self.config.address)
        path.mkdir(parents=True, exist_ok=True)
        options = dict(cleanup_policy=self.config.cleanup_policy,
                       consolidation_batch_size=self.config.consolidation_batch_size)
        grid_writer = StoreSolution(str(path), "grid", self.config.format, **options)
        grid_writer.save_file(list(self.axes.values()), ti=None)
        grid_writer.close_file("end_grid")
        metadata = {"schema_version": 1, "streams": {}}
        for name, stream in self.streams.items():
            diag = stream["geometry"] == "diagnostics"
            stream["writer"] = StoreSolution(str(path), name, self.config.format,
                diagnostic_names=[stream["variable"]] if diag else None, **options)
            axes = "" if diag else stream["orientation"] or "xyz"
            metadata["streams"][name] = {
                "geometry": stream["geometry"], "variable": stream["variable"],
                "orientation": stream["orientation"], "axes": list(axes),
                "fixed_coordinates": {} if diag else {
                    a: float(self.axes[a][self.center[a]]) for a in "xyz" if a not in axes},
                "fixed_indices": {} if diag else {a: self.center[a] for a in "xyz" if a not in axes},
                "format": self.config.format,
            }
        (path / "output_manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")

    def outputs_at(self, index, time):
        return {name: stream for name, stream in self.streams.items()
                if any(rule.matches(index, self.total, time, self.duration) for rule in stream["rules"])}

    @staticmethod
    def diagnostic_flags(selected):
        names = {s["variable"] for s in selected.values() if s["geometry"] == "diagnostics"}
        return {name: name in names for name in ("Numb_Part", "Energ", "Pi", "Ji", "Frequency")}

    def save(self, selected, index, time, rho, psi, phi, diagnostics):
        fields = {"rho": rho, "psi": psi, "phi": phi}
        for stream in selected.values():
            if stream["geometry"] == "diagnostics":
                value = np.empty(1, dtype=object)
                value[0] = diagnostics[stream["variable"]]
            else:
                varying = stream["orientation"] or "xyz"
                slices = tuple(slice(None) if a in varying else self.center[a] for a in "xyz")
                value = fields[stream["variable"]][(Ellipsis,) + slices]
            stream["writer"].save_file(value, index, physical_time=time)

    def close(self):
        for name, stream in self.streams.items():
            stream["writer"].close_file("end_" + name)


def read_output(address, name):
    """Read one selected stream (NPZ/HDF5) with original indices and times.

    Returns metadata, coordinate vectors, time, snapshot_index, and data.
    Fields are stacked along a leading sample axis. Diagnostics are a mapping
    of numeric arrays. Only read NPZ diagnostics from trusted simulations:
    the legacy object-array encoding requires pickle.
    """
    path = Path(address)
    metadata = json.loads((path / "output_manifest.json").read_text())["streams"][name]
    diag = metadata["geometry"] == "diagnostics"
    fmt = metadata["format"]
    filename = path / f"end_{name}.{fmt}"
    if fmt == "npz":
        with np.load(filename, allow_pickle=diag) as archive:
            indices = archive["snapshot_index"].astype(int)
            times = archive["time"] if "time" in archive else np.empty(0)
            samples = [archive[f"{name}{i}"] for i in indices]
        if diag:
            samples = [_diagnostic_datasets(s, [metadata["variable"]]) for s in samples]
        with np.load(path / "end_grid.npz") as grid:
            coordinates = {a: grid[a] for a in metadata["axes"]}
    else:
        with h5py.File(filename) as archive:
            indices = archive["snapshot_index"][:].astype(int)
            times = archive["time"][:] if "time" in archive else np.empty(0)
            samples = [({k: v[()] for k, v in archive[f"{name}_{i}"].items()} if diag
                        else archive[f"{name}_{i}"][()]) for i in indices]
        with h5py.File(path / "end_grid.hdf5") as grid:
            coordinates = {a: grid[a][:] for a in metadata["axes"]}
    data = ({k: np.stack([s[k] for s in samples]) for k in samples[0]} if samples else {}) if diag else (
        np.stack(samples) if samples else np.empty((0,)))
    return {**metadata, "coordinates": coordinates, "snapshot_index": indices,
            "time": times, "data": data}
