"""Recoverable, single-writer consolidation of snapshot archives.

An incremental archive is still mutable: its journal cannot repair corruption.
Never run two consolidators or modify snapshots in the same output concurrently.
"""

import glob
import hashlib
import json
import os
import shutil
import tempfile
import warnings
import zipfile

import h5py
import numpy as np


def validate_cleanup(cleanup_policy, consolidation_batch_size):
    if cleanup_policy not in ("after_success", "incremental"):
        raise ValueError("cleanup_policy must be after_success or incremental")
    if (isinstance(consolidation_batch_size, bool)
            or not isinstance(consolidation_batch_size, int)
            or consolidation_batch_size <= 0):
        raise ValueError("consolidation_batch_size must be a positive integer")


def _sync_directory(path):
    if os.name == "posix":
        fd = os.open(os.path.dirname(path) or ".", os.O_RDONLY)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)


def _sync_file(path):
    with open(path, "rb") as stream:
        os.fsync(stream.fileno())


def _write_state(path, state):
    fd, temporary = tempfile.mkstemp(prefix=".consolidation-", dir=os.path.dirname(path))
    try:
        with os.fdopen(fd, "w") as stream:
            json.dump(state, stream)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        _sync_directory(path)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)


def _signature(path):
    stat = os.stat(path)
    return [stat.st_size, stat.st_mtime_ns, stat.st_dev, stat.st_ino]


def _open(path, fmt, mode):
    if fmt == "npz":
        return zipfile.ZipFile(path, mode, compression=zipfile.ZIP_DEFLATED, allowZip64=True)
    return h5py.File(path, mode)


def _members(archive, fmt):
    return archive.namelist() if fmt == "npz" else list(archive.keys())


def _hash_value(digest, value):
    array = np.asarray(value)
    digest.update(str((array.shape, array.dtype.str)).encode())
    if array.dtype.hasobject:
        # Variable-length HDF5 strings must be hashed by value, not pointers.
        digest.update(repr(array.tolist()).encode())
    else:
        digest.update(array.tobytes())


def _digest(archive, name, fmt):
    digest = hashlib.sha256()
    if fmt == "npz":
        with archive.open(name) as stream:
            while True:
                block = stream.read(1024 * 1024)
                if not block:
                    break
                digest.update(block)
    else:
        def visit(relative, obj):
            digest.update(relative.encode())
            digest.update(b"dataset" if isinstance(obj, h5py.Dataset) else b"group")
            for key in sorted(obj.attrs):
                digest.update(key.encode())
                _hash_value(digest, obj.attrs[key])
            if isinstance(obj, h5py.Dataset):
                digest.update(str((obj.shape, obj.dtype.str)).encode())
                if obj.shape is None:
                    return
                if obj.ndim == 0:
                    _hash_value(digest, obj[()])
                else:
                    row_bytes = max(1, int(np.prod(obj.shape[1:])) * obj.dtype.itemsize)
                    rows = max(1, (1024 * 1024) // row_bytes)
                    for start in range(0, obj.shape[0], rows):
                        _hash_value(digest, obj[start:start + rows])
        root = archive[name]
        visit("", root)
        if isinstance(root, h5py.Group):
            root.visititems(visit)
    return digest.hexdigest()


def _copy_batch(temporary, batch, fmt):
    expected = {}
    with _open(temporary, fmt, "a") as target:
        for entry in batch:
            if _signature(entry["path"]) != entry["signature"]:
                raise ValueError(f"Source changed: {entry['path']}")
            with _open(entry["path"], fmt, "r") as source:
                for original, name in entry["members"]:
                    expected[name] = _digest(source, original, fmt)
                    if name in _members(target, fmt):
                        # A prior attempt may have closed the batch before journaling it.
                        if _digest(target, name, fmt) != expected[name]:
                            raise ValueError(f"Incomplete or corrupt member: {name}; keep recovery files")
                        continue
                    if fmt == "npz":
                        with source.open(original) as reader, target.open(name, "w", force_zip64=True) as writer:
                            shutil.copyfileobj(reader, writer, length=1024 * 1024)
                    else:
                        source.copy(original, target, name=name)
    _sync_file(temporary)
    _verify(temporary, fmt, expected)
    return expected


def _verify(path, fmt, expected, exact=False):
    with _open(path, fmt, "r") as archive:
        names = _members(archive, fmt)
        if len(names) != len(set(names)) or (exact and set(names) != set(expected)):
            raise ValueError("Consolidated archive has unexpected or duplicate members")
        for name, checksum in expected.items():
            if name not in names or _digest(archive, name, fmt) != checksum:
                raise ValueError(f"Consolidated archive failed verification: {name}")


def _remove_sources(entries):
    for entry in entries:
        path = entry["path"]
        if os.path.exists(path):
            if _signature(path) != entry["signature"]:
                raise ValueError(f"Refusing to delete changed source: {path}")
            os.remove(path)
            _sync_directory(path)


def consolidate(file_list, output, fmt, cleanup_policy="after_success", consolidation_batch_size=16):
    """Publish a verified archive, retaining a journal on failure for retry.

    Retry with the same output and policy. A readable partial archive can be
    resumed; a corrupt one is retained and raises without deleting more sources.
    """
    validate_cleanup(cleanup_policy, consolidation_batch_size)
    output = os.path.realpath(output)
    journal = output + ".consolidation.json"
    if cleanup_policy == "incremental":
        warnings.warn("Incremental cleanup deletes snapshots before final publication; "
                      "archive corruption can lose already deleted data.", RuntimeWarning, stacklevel=2)
    filenames = glob.glob(file_list) if isinstance(file_list, str) else list(file_list)
    filenames = sorted(os.path.realpath(name) for name in filenames)
    if output in filenames or any(os.path.exists(output) and os.path.samefile(name, output)
                                  for name in filenames if os.path.exists(name)):
        raise ValueError("Destination must not be a source")
    if os.path.exists(journal):
        with open(journal) as stream:
            state = json.load(stream)
        if (state["output"], state["format"], state["policy"]) != (output, fmt, cleanup_policy):
            raise ValueError("Recovery requires the original output, format and cleanup policy")
        if set(filenames) - {entry["path"] for entry in state["sources"]}:
            raise ValueError("New sources found during recovery; use a separate output")
    else:
        if not filenames:
            raise ValueError("No snapshot files to consolidate")
        entries, seen = [], set()
        for path in filenames:
            with _open(path, fmt, "r") as source:
                names = _members(source, fmt)
            members = [(name, name[:-4] if fmt == "npz" and name.endswith(".npy") else name)
                       for name in names]
            for _, name in members:
                if name in seen:
                    raise ValueError(f"Duplicate archive member: {name}")
                seen.add(name)
            entries.append({"path": path, "signature": _signature(path), "members": members})
        fd, temporary = tempfile.mkstemp(prefix=".consolidation-", suffix=".partial", dir=os.path.dirname(output))
        os.close(fd)
        with _open(temporary, fmt, "w"):
            pass
        _sync_file(temporary)
        state = {"output": output, "format": fmt, "policy": cleanup_policy,
                 "temporary": temporary, "sources": entries, "completed": 0,
                 "checksums": {}, "ready": False}
        _write_state(journal, state)

    temporary = state["temporary"]
    # A crash can occur between atomic publication and journal removal.
    if not os.path.exists(temporary):
        if not state["ready"]:
            raise ValueError("Recovery temporary is missing; originals were not cleaned further")
        _verify(output, fmt, state["checksums"], exact=True)
        _sync_file(output)
        _sync_directory(output)
    else:
        _verify(temporary, fmt, state["checksums"])
        if cleanup_policy == "incremental":
            _remove_sources(state["sources"][:state["completed"]])
        while state["completed"] < len(state["sources"]):
            start = state["completed"]
            batch = state["sources"][start:start + consolidation_batch_size]
            checksums = _copy_batch(temporary, batch, fmt)
            state["checksums"].update(checksums)
            state["completed"] += len(batch)
            _write_state(journal, state)
            if cleanup_policy == "incremental":
                _remove_sources(batch)
        _verify(temporary, fmt, state["checksums"], exact=True)
        state["ready"] = True
        _write_state(journal, state)
        os.replace(temporary, output)
        _sync_directory(output)
    _remove_sources(state["sources"])
    os.remove(journal)
    _sync_directory(output)
