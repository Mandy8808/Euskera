"""Regression tests for HDF5/npz output consolidation in euskera.io.save_data."""

import json
import zipfile
from pathlib import Path

import h5py
import numpy as np
import pytest

from euskera.io import OutputConfig
from euskera.io import consolidation as cs
from euskera.io.save_data import StoreSolution, data_Objgenerator


def test_hdf5_close_file_consolidates_single_result(tmp_path):
    """Saving a single HDF5 timestep must not raise BadZipFile on close.

    HDF5 files are not zip archives, so close_file must not try to merge
    them with the npz-only zip-merging logic.
    """
    obj = StoreSolution(address=str(tmp_path), filename="save_rho", format="hdf5")
    obj.save_file(np.zeros((2, 2, 2)), ti=0)
    obj.close_file("end_save_rho")

    archive = tmp_path / "end_save_rho.hdf5"
    assert archive.exists()
    with h5py.File(archive, "r") as f:
        assert "save_rho_0" in f
        assert "t" in f

    # Intermediate per-step files should have been cleaned up.
    assert not list(tmp_path.glob("save_rho_u_*.h5"))


def test_hdf5_close_file_consolidates_multiple_results(tmp_path):
    """Multiple HDF5 timesteps should be merged into a single archive."""
    obj = StoreSolution(address=str(tmp_path), filename="save_rho", format="hdf5")
    for ti, value in enumerate([0.0, 1.0, 2.0]):
        obj.save_file(np.full((2, 2, 2), value), ti=ti)
    obj.close_file("end_save_rho")

    archive = tmp_path / "end_save_rho.hdf5"
    with h5py.File(archive, "r") as f:
        for ti, value in enumerate([0.0, 1.0, 2.0]):
            assert np.allclose(f[f"save_rho_{ti}"][()], value)
        assert np.allclose(f["t"][()], [0.0, 1.0, 2.0])

    assert not list(tmp_path.glob("save_rho_u_*.h5"))


def test_npz_close_file_still_consolidates(tmp_path):
    """Existing npz consolidation behavior must be unaffected."""
    obj = StoreSolution(address=str(tmp_path), filename="save_rho", format="npz")
    obj.save_file(np.zeros((2, 2, 2)), ti=0)
    obj.save_file(np.ones((2, 2, 2)), ti=1)
    obj.close_file("end_save_rho")

    archive = tmp_path / "end_save_rho.npz"
    assert archive.exists()
    with np.load(archive) as data:
        assert "save_rho0" in data.files
        assert "save_rho1" in data.files

    assert not list(tmp_path.glob("save_rho_u_*.dat.npz"))


@pytest.fixture(params=["npz", "hdf5"])
def snapshots(request, tmp_path):
    fmt = request.param
    store = StoreSolution(str(tmp_path), "save_rho", format=fmt,
                          consolidation_batch_size=1)
    for step in range(3):
        store.save_file(np.full((2, 3), step, dtype=float), step, physical_time=step * .25)
    extension = "npz" if fmt == "npz" else "h5"
    sources = sorted(tmp_path.glob(f"save_rho_u_*.{extension}"))
    output = tmp_path / f"end_save_rho.{fmt}"
    return store, sources, output


def read_values(output, fmt):
    if fmt == "npz":
        with np.load(output) as archive:
            return {name: archive[name] for name in archive.files}
    with h5py.File(output) as archive:
        return {name: archive[name][()] for name in archive}


@pytest.mark.parametrize("policy", ["after_success", "incremental"])
def test_cleanup_timing_and_content(snapshots, monkeypatch, policy):
    store, sources, output = snapshots
    store.cleanup_policy = policy
    original = cs._copy_batch
    calls = []

    def observe(path, batch, fmt):
        calls.append(batch)
        if len(calls) == 2:
            assert Path(calls[0][0]["path"]).exists() == (policy == "after_success")
            assert not output.exists()
        return original(path, batch, fmt)

    monkeypatch.setattr(cs, "_copy_batch", observe)
    store.close_file("end_save_rho")
    values = read_values(output, store.format)
    for step in range(3):
        key = f"save_rho{step}" if store.format == "npz" else f"save_rho_{step}"
        np.testing.assert_array_equal(values[key], np.full((2, 3), step))
    np.testing.assert_array_equal(values["time"], [0, .25, .5])
    assert all(not source.exists() for source in sources)
    assert not list(output.parent.glob("*.consolidation.json"))
    assert not list(output.parent.glob(".consolidation-*"))


@pytest.mark.parametrize("policy", ["after_success", "incremental"])
def test_copy_failure_preserves_old_output_and_resumes(snapshots, monkeypatch, policy):
    store, sources, output = snapshots
    output.write_bytes(b"previous result")
    copy = cs._copy_batch
    calls = 0

    def fail_second(*args):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("simulated disk full")
        return copy(*args)

    monkeypatch.setattr(cs, "_copy_batch", fail_second)
    with pytest.raises(OSError, match="disk full"):
        cs.consolidate(sources, output, store.format, policy, 1)
    assert output.read_bytes() == b"previous result"
    assert sources[0].exists() == (policy == "after_success")
    assert sources[1].exists() and sources[2].exists()
    monkeypatch.setattr(cs, "_copy_batch", copy)
    cs.consolidate(sources, output, store.format, policy, 1)
    assert len(read_values(output, store.format)) == 3
    assert all(not source.exists() for source in sources)


@pytest.mark.parametrize("policy", ["after_success", "incremental"])
def test_verification_failure_does_not_delete_batch(snapshots, monkeypatch, policy):
    store, sources, output = snapshots
    verify = cs._verify

    def fail_new_batch(path, fmt, expected, exact=False):
        if expected:
            raise ValueError("simulated verification failure")
        return verify(path, fmt, expected, exact)

    monkeypatch.setattr(cs, "_verify", fail_new_batch)
    with pytest.raises(ValueError, match="verification failure"):
        cs.consolidate(sources, output, store.format, policy, 1)
    assert all(source.exists() for source in sources)
    assert not output.exists()
    monkeypatch.setattr(cs, "_verify", verify)
    cs.consolidate(sources, output, store.format, policy, 1)
    assert len(read_values(output, store.format)) == 3


@pytest.mark.parametrize("policy", ["after_success", "incremental"])
def test_publication_failure_can_retry_with_new_store(snapshots, monkeypatch, policy):
    store, sources, output = snapshots
    store.cleanup_policy = policy
    output.write_bytes(b"previous result")
    replace = cs.os.replace

    def fail_publish(source, destination):
        if str(destination) == str(output):
            raise OSError("publication interrupted")
        return replace(source, destination)

    monkeypatch.setattr(cs.os, "replace", fail_publish)
    with pytest.raises(OSError, match="publication interrupted"):
        store.close_file("end_save_rho")
    assert output.read_bytes() == b"previous result"
    assert all(source.exists() == (policy == "after_success") for source in sources)
    monkeypatch.setattr(cs.os, "replace", replace)
    fresh = StoreSolution(str(output.parent), "save_rho", format=store.format,
                          cleanup_policy=policy)
    fresh.close_file("end_save_rho")
    values = read_values(output, store.format)
    np.testing.assert_array_equal(values["time"], [0, .25, .5])


def test_cleanup_failure_after_publication_resumes(snapshots, monkeypatch):
    store, sources, output = snapshots
    remove = cs.os.remove

    def fail_remove(path):
        if str(path) == str(sources[1]):
            raise OSError("cleanup interrupted")
        return remove(path)

    monkeypatch.setattr(cs.os, "remove", fail_remove)
    with pytest.raises(OSError, match="cleanup interrupted"):
        cs.consolidate(sources, output, store.format)
    assert output.exists() and not sources[0].exists()
    monkeypatch.setattr(cs.os, "remove", remove)
    cs.consolidate([], output, store.format)
    assert all(not source.exists() for source in sources)


def test_corrupt_checkpoint_stops_without_more_deletions(snapshots, monkeypatch):
    store, sources, output = snapshots
    copy = cs._copy_batch
    count = 0

    def fail_second(*args):
        nonlocal count
        count += 1
        if count == 2:
            raise OSError("interrupted")
        return copy(*args)

    monkeypatch.setattr(cs, "_copy_batch", fail_second)
    with pytest.raises(OSError):
        cs.consolidate(sources, output, store.format, "incremental", 1)
    state = json.loads((output.parent / (output.name + ".consolidation.json")).read_text())
    with open(state["temporary"], "wb") as stream:
        stream.write(b"broken")
    monkeypatch.setattr(cs, "_copy_batch", copy)
    with pytest.raises((OSError, ValueError, zipfile.BadZipFile)):
        cs.consolidate([], output, store.format, "incremental", 1)
    assert sources[1].exists() and sources[2].exists()


def test_rejects_destination_duplicate_and_empty_sources(snapshots):
    store, sources, output = snapshots
    with pytest.raises(ValueError, match="Destination"):
        cs.consolidate(sources, sources[0], store.format)
    with pytest.raises(ValueError, match="Duplicate"):
        cs.consolidate([sources[0], sources[0]], output, store.format)
    with pytest.raises(ValueError, match="No snapshot"):
        cs.consolidate([], output, store.format)
    assert all(source.exists() for source in sources)


@pytest.mark.parametrize("options", [
    {"cleanup_policy": "unknown"}, {"consolidation_batch_size": 0},
    {"consolidation_batch_size": True}, {"consolidation_batch_size": 1.5},
])
def test_invalid_cleanup_configuration(tmp_path, options):
    with pytest.raises(ValueError):
        OutputConfig(**options)
    with pytest.raises(ValueError):
        StoreSolution(str(tmp_path), "rho", **options)


def test_cleanup_options_propagate_to_stores(tmp_path):
    options = OutputConfig(cleanup_policy="incremental", consolidation_batch_size=2)
    assert OutputConfig.from_legacy(options.to_dict()) == options
    stores = data_Objgenerator({"save_rho": True}, str(tmp_path), "npz",
                               cleanup_policy=options.cleanup_policy,
                               consolidation_batch_size=options.consolidation_batch_size)
    assert stores["save_rho"].cleanup_policy == "incremental"
    assert stores["save_rho"].consolidation_batch_size == 2


def test_changed_source_is_not_deleted_on_retry(snapshots, monkeypatch):
    store, sources, output = snapshots
    copy = cs._copy_batch
    calls = 0

    def fail_second(*args):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("interrupted")
        return copy(*args)

    monkeypatch.setattr(cs, "_copy_batch", fail_second)
    with pytest.raises(OSError):
        cs.consolidate(sources, output, store.format, "after_success", 1)
    sources[1].write_bytes(b"replacement must survive")
    monkeypatch.setattr(cs, "_copy_batch", copy)
    with pytest.raises(ValueError, match="Source changed"):
        cs.consolidate([], output, store.format)
    assert all(source.exists() for source in sources)
    assert sources[1].read_bytes() == b"replacement must survive"
    assert not output.exists()


def test_incremental_warning_is_explicit(snapshots):
    store, sources, output = snapshots
    with pytest.warns(RuntimeWarning, match="corruption can lose"):
        cs.consolidate(sources, output, store.format, "incremental")


def test_repeated_close_preserves_completed_archive(snapshots):
    store, sources, output = snapshots
    store.close_file("end_save_rho")
    original = output.read_bytes()
    store.close_file("end_save_rho")
    fresh = StoreSolution(str(output.parent), "save_rho", format=store.format)
    fresh.close_file("end_save_rho")
    assert output.read_bytes() == original
