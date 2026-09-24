"""DATA SAVING UTILITIES FOR EUSKERA SIMULATIONS"""

import os
import h5py
import glob
import zipfile
import numpy as np

from .consolidation import consolidate, validate_cleanup

########### Save data
#############################################################################
def data_Objgenerator(data_save, address, format, comp_conserv=None,
                      cleanup_policy="after_success", consolidation_batch_size=16):
    """
    Generates a dictionary of StoreSolution objects for enabled data saving options.

    Parameters:
    - data_save (dict): Dictionary of boolean flags indicating which data to save.
    - address (str): Directory where data will be stored.
    - format (str): File format (e.g., "npz", "hdf5").
    - comp_conserv (dict, optional): Diagnostic flags in Conserv output order.
    - cleanup_policy: "after_success" (default) or riskier "incremental" cleanup.
    - consolidation_batch_size: Positive source-file count per verified batch.

    Returns:
    - dict: Dictionary with keys as data types and values as StoreSolution objects.
    """
    data_save_obj = {
                    name: StoreSolution(address=address, filename=name, format=format,
                        cleanup_policy=cleanup_policy, consolidation_batch_size=consolidation_batch_size,
                        diagnostic_names=None if comp_conserv is None else [
                            key for key, enabled in comp_conserv.items()
                            if enabled and key in ("Numb_Part", "Energ", "Pi", "Frequency")])
                          for name, save in data_save.items() if save}
    return data_save_obj

def fdata_save(ti, data, data_save_obj, resol, end=False, physical_time=None):
    """
    Saves simulation data.

    Parameters:
    - ti: Temporal index.
    - data: Tuple containing ([x, y, z], rho, psi, phi, cData).
    - data_save_obj: Dictionary specifying the names of the data to be saved.
    - resol: Spatial resolution used in the simulation.
    - end: If True, the files are closed with the subname 'end_'.
    """

    [xarray, yarray, zarray], rho, psi, phi, cData = data

    if end:
        for name, obj in data_save_obj.items():
            obj.close_file("end_" + name)
        print("\n All data was save")
    else:
        save_map = {
            "grid": "[xarray[:, 0, 0], yarray[0, :, 0], zarray[0, 0, :]] if xarray is not None else None",
            #
            "save_rho": "rho if rho is not None else None",
            "save_psi": "psi if psi is not None else None",
            "save_phi": "phi if phi is not None else None",
            #
            "save_plane_rho": "rho[:, :, resol // 2] if rho is not None else None",
            "save_plane_psi": "psi[:, :, :, resol // 2] if psi is not None else None",
            "save_plane_phi": "phi[:, :, resol // 2] if phi is not None else None",
            #
            "save_line_rho": "rho[:, resol // 2, resol // 2] if rho is not None else None",
            "save_line_psi": "psi[:, :, resol // 2, resol // 2] if psi is not None else None",
            "save_line_phi": "phi[:, resol // 2, resol // 2] if phi is not None else None",
            #
            "save_energies": "cData if cData is not None else None",
        }
        #fdata_save(ti=0, data=data, data_save_obj=data_save_obj, resol=resol, end=False)
        for name, obj in data_save_obj.items():
            if name in save_map:
                if eval(save_map[name]) is None and name != "grid":  # the second condition is that after save the grid, its value change to None
                    print(f"Warning: Skipping {name} save because data is missing.")
                    continue
                if name == "grid" and ti != 0:
                    continue
                obj.save_file(eval(save_map[name]), ti=None if name == "grid" else ti, physical_time=physical_time)
    return None

def nameData(address, file_format, info=True):
    """
    Retrieves the names of files with a specific format in the given directory.

    Parameters:
    - address (str): The directory path to search.
    - file_format (str): The file extension to filter (e.g., 'npz', 'hdf5').
    - info (bool): If True, prints the number of files found.

    Returns:
    - List of filenames matching the given format.
    """
    DataName = []
    count = 0  # Track the number of matching files

    for file in os.listdir(address):
        if file.endswith(f".{file_format}"):
            DataName.append(os.path.basename(file))
            count += 1

    if info:
        #print(f"Se encontraron {count} archivos con la extensión .{file_format}")
        print(f"Found {count} files with the .{file_format} extension.")
    return DataName

def JoinFilesInOneZip(file_list, output_zip, cleanup_policy="after_success", consolidation_batch_size=16):
    """Consolidate NPZ snapshots; retain sources until publication by default.

    Incremental cleanup deletes verified batches before publication and can lose
    data if the partial archive is corrupted. Retry with the same output/policy.
    """
    consolidate(file_list, output_zip, "npz", cleanup_policy, consolidation_batch_size)


def JoinFilesInOneHDF5(file_list, output_hdf5, cleanup_policy="after_success", consolidation_batch_size=16):
    """Consolidate HDF5 snapshots using the same cleanup policies as NPZ."""
    consolidate(file_list, output_hdf5, "hdf5", cleanup_policy, consolidation_batch_size)

class StoreSolution:
    """
    Class to save the results.
    """

    def __init__(self, address, filename, format="npz", info=False, diagnostic_names=None,
                 cleanup_policy="after_success", consolidation_batch_size=16):
        """
        Initialize the storage class.

        Parameters:
        - address (str): Directory to save files.
        - filename (str): Base filename for saved data.
        - file_format (str): Format to save files ("npz" or "hdf5").
        - info (bool): Print information about file saving.
        - diagnostic_names: Enabled Conserv names in output order. Required
          for HDF5 save_energies; ignored for NPZ and numeric field datasets.
        - cleanup_policy: Delete sources after publication ("after_success") or
          after each verified batch ("incremental", with a data-loss risk).
        - consolidation_batch_size: Positive source-file count per batch.
        """
        if info:
            print(f"Using address: {address}, Filename: {filename}, Format: {format}")

        validate_cleanup(cleanup_policy, consolidation_batch_size)
        self.cleanup_policy = cleanup_policy
        self.consolidation_batch_size = consolidation_batch_size

        # Instance attributes
        self.address = address
        self.name = filename
        self.format = format.lower()  # Normalize format to lowercase
        self.time = []
        self.physical_time = []
        self.diagnostic_names = diagnostic_names

        # Ensure directory exists
        os.makedirs(self.address, exist_ok=True)

    def save_file(self, data, ti, physical_time=None):
        """
        Save data to a file.

        Parameters:
        - data: Data to be saved.
        - ti: Time step.
        - n (int): Index number (-1 for grid data, other values for simulation steps).
        """
        if self.format not in ["npz", "hdf5"]:
            print("WARNING: Unsupported format. Defaulting to 'npz'.")
            self.format = "npz"

        if ti is None:
            # Save grid data
            if self.format == "npz":
                kwargs = {"x": data[0], "y": data[1], "z": data[2]}
                np.savez(os.path.join(self.address, f"{self.name}_grid.dat"), **kwargs)
            elif self.format == "hdf5":
                with h5py.File(os.path.join(self.address, f"{self.name}_grid.h5"), "w") as f:
                    f.create_dataset("x", data=data[0])
                    f.create_dataset("y", data=data[1])
                    f.create_dataset("z", data=data[2])
        else:
            # Save simulation step data
            if self.format == "npz":
                kwargs = {f"{self.name}{ti}": data}
                fname = os.path.join(self.address, f"{self.name}_u_{ti}.dat.npz")
                np.savez(fname, **kwargs)
            elif self.format == "hdf5":
                diagnostics = None
                if self.name == "save_energies" or self.name.startswith("save_energies_"):
                    diagnostics = _diagnostic_datasets(data, self.diagnostic_names)
                fname = os.path.join(self.address, f"{self.name}_u_{ti}.h5")
                with h5py.File(fname, "w") as f:
                    if diagnostics is None:
                        f.create_dataset(f"{self.name}_{ti}", data=data)
                    else:
                        group = f.create_group(f"{self.name}_{ti}")
                        group.attrs["schema_version"] = 1
                        group.attrs["diagnostic_names"] = np.asarray(
                            self.diagnostic_names, dtype=h5py.string_dtype())
                        for key, value in diagnostics.items():
                            group.create_dataset(key, data=value)

            # Store time step
            self.time.append(ti)
            self.physical_time.append(physical_time)

    def close_file(self, zip_name):
        """
        Combine snapshots into one verified NPZ or HDF5 archive.

        Parameters:
        - zip_name (str): Name of the final zip archive (without extension).
        """
        if self.address is not None:
            archive_name = os.path.join(self.address, f"{zip_name}.{self.format}")
            # Recovery must use the original metadata, including physical times.
            journal = os.path.realpath(archive_name) + ".consolidation.json"
            if os.path.exists(journal):
                consolidate([], archive_name, self.format, self.cleanup_policy,
                            self.consolidation_batch_size)
                return
            file_extension = "*.npz" if self.format == "npz" else "*.h5"
            pattern = os.path.join(self.address, f"{self.name}_*{file_extension}")
            # A repeated close after successful cleanup must not replace data
            # with an archive containing only newly generated time metadata.
            if os.path.exists(archive_name) and not glob.glob(pattern):
                return
            times = {"t": np.array(self.time, dtype=float),
                     "snapshot_index": np.array(self.time, dtype=int)}
            if self.physical_time and all(t is not None for t in self.physical_time):
                times["time"] = np.array(self.physical_time, dtype=float)
            # Save all the time points where solutions are saved
            if self.format == "npz" and self.name != "grid":
                np.savez(os.path.join(self.address, f"{self.name}_t.dat"), **times)
            elif self.format == "hdf5" and self.name != "grid":
                with h5py.File(os.path.join(self.address, f"{self.name}_t.h5"), "w") as f:
                    for key, value in times.items():
                        f.create_dataset(key, data=value)

            # Find all files of the selected format
            filenames = glob.glob(pattern)

            # Archive files
            if self.format == "npz":
                JoinFilesInOneZip(filenames, archive_name, self.cleanup_policy, self.consolidation_batch_size)
            else:
                JoinFilesInOneHDF5(filenames, archive_name, self.cleanup_policy, self.consolidation_batch_size)


def _diagnostic_datasets(data, names):
    """Convert Conserv's ordered object array into named numeric datasets."""
    if names is None or len(names) != len(data):
        raise ValueError("HDF5 diagnostics require diagnostic_names matching the Conserv output order.")
    result = {}
    for name, value in zip(names, data):
        if name == "Numb_Part":
            result["mass_total"] = np.asarray(value[0], dtype=float)
            result["mass_components"] = np.asarray(value[1], dtype=float)
        elif name == "Energ":
            result["energy"] = np.asarray(value, dtype=float)
        elif name == "Pi":
            result["momentum"] = np.asarray(value, dtype=float)
        elif name == "Frequency":
            result["frequency_samples"] = np.asarray(value, dtype=complex)
        else:
            raise ValueError(f"Unsupported diagnostic: {name}")
    return result


def read_hdf5_diagnostics(path):
    """Return {snapshot_index: {quantity: numeric value/array}} for schema 1.

    Indices are save counters, not physical times. frequency_samples contains
    complex wavefunction samples, not estimated frequencies. Disabled
    quantities are absent; a snapshot with all diagnostics disabled is empty.
    """
    snapshots = {}
    with h5py.File(path, "r") as archive:
        names = sorted((name for name in archive if name.startswith("save_energies_")),
                       key=lambda name: int(name.rsplit("_", 1)[1]))
        for name in names:
            group = archive[name]
            if not isinstance(group, h5py.Group) or group.attrs.get("schema_version") != 1:
                raise ValueError(f"Unsupported diagnostics schema in {name}")
            snapshots[int(name.rsplit("_", 1)[1])] = {
                key: dataset[()] for key, dataset in group.items()}
    return snapshots


################ Old versions

def JoinFilesInOneZip_old(file_list, output_zip):
    """
    Uniendo todos los archivos en un .zip

    IN:
    output_zip -> nombre que daremos al .zip
    file_list -> lista o tupla de datos, o dirección donde encontraremos los datos

    Out:
    Crea un archivo: archive_name.zip
    """

    # creando el .zip
    archive = zipfile.ZipFile(output_zip, mode='w', compression=zipfile.ZIP_DEFLATED, allowZip64=True)

    # identificando los archivos a comprimir
    if isinstance(file_list, (list, tuple)):
        filenames = file_list
    elif isinstance(file_list, str):
        filenames = glob.glob(file_list)

    # Open each archive and write to the common archive
    for filename in filenames:
        f = zipfile.ZipFile(filename, mode='r', compression=zipfile.ZIP_DEFLATED)
        for name in f.namelist():
            data = f.open(name, 'r')
            # Save under name without .npy
            archive.writestr(name[:-4], data.read())
        f.close()

        # Remove original files
        os.remove(filename)
    archive.close()

def extract_zip_to_memory(filename):
    """
    Extrae los archivos desde un zip a un diccionario en memoria
    """
    extracted = []
    with zipfile.ZipFile(filename, 'r') as f:
        for name in f.namelist():
            data = f.read(name)
            # Save under name without .npy
            extracted.append((name[:-4], data))
    return (filename, extracted)

def JoinFilesInOneZip_old(file_list, output_zip):
    """
    Une todos los archivos de múltiples .zip en un solo .zip, usando multiprocesamiento
    """
    from concurrent.futures import ProcessPoolExecutor

    # Identificar los archivos a comprimir
    if isinstance(file_list, (list, tuple)):
        filenames = file_list
    elif isinstance(file_list, str):
        filenames = glob.glob(file_list)

    # Extraer en paralelo
    results = []
    with ProcessPoolExecutor() as executor:
        for result in executor.map(extract_zip_to_memory, filenames):
            results.append(result)

    # Escribir en un solo zip de forma secuencial
    with zipfile.ZipFile(output_zip, mode='w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as archive:
        for original_file, extracted_files in results:
            for name, data in extracted_files:
                archive.writestr(name, data)

    # Eliminar los zips originales
    for f in filenames:
        os.remove(f)
