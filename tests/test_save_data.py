"""Regression tests for HDF5/npz output consolidation in euskera.io.save_data."""

import h5py
import numpy as np

from euskera.io.save_data import StoreSolution


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
