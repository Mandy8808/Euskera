"""DATA SAVING UTILITIES FOR EUSKERA SIMULATIONS"""

import os
import h5py
import glob
import zipfile
import numpy as np

########### Save data
#############################################################################
def data_Objgenerator(data_save, address, format):
    """
    Generates a dictionary of StoreSolution objects for enabled data saving options.

    Parameters:
    - data_save (dict): Dictionary of boolean flags indicating which data to save.
    - address (str): Directory where data will be stored.
    - format (str): File format (e.g., "npz", "hdf5").

    Returns:
    - dict: Dictionary with keys as data types and values as StoreSolution objects.
    """
    data_save_obj = {
                    name: StoreSolution(address=address, filename=name, format=format)
                          for name, save in data_save.items() if save}
    return data_save_obj

def fdata_save(ti, data, data_save_obj, resol, end=False):
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
                obj.save_file(eval(save_map[name]), ti=None if name == "grid" else ti)
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

def JoinFilesInOneZip(file_list, output_zip):
    """
    Uniendo todos los archivos en un .zip

    IN:
    output_zip -> nombre que daremos al .zip
    file_list -> lista o tupla de datos, o dirección donde encontraremos los datos

    Out:
    Crea un archivo: archive_name.zip
    """

    # identificando los archivos a comprimir
    if isinstance(file_list, (list, tuple)):
        filenames = file_list
    elif isinstance(file_list, str):
        filenames = glob.glob(file_list)

    with zipfile.ZipFile(output_zip, mode='w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as archive:
        for filename in filenames:
            with zipfile.ZipFile(filename, mode='r') as f:
                for name in f.namelist():
                    # Save under name without .npy
                    archive.writestr(name[:-4], f.read(name))  # leer y escribir directamente

            # Remove original files
            os.remove(filename)

def JoinFilesInOneHDF5(file_list, output_hdf5):
    """
    Uniendo todos los archivos HDF5 en un unico archivo HDF5.

    Unlike .npz, HDF5 files are not zip archives, so they cannot be merged
    with zipfile. Each dataset is copied into a single consolidated file.

    IN:
    output_hdf5 -> nombre del archivo HDF5 final
    file_list -> lista o tupla de archivos .h5, o patron para localizarlos

    Out:
    Crea un archivo: output_hdf5 con todos los datasets de los archivos originales.
    """

    # identificando los archivos a combinar
    if isinstance(file_list, (list, tuple)):
        filenames = file_list
    elif isinstance(file_list, str):
        filenames = glob.glob(file_list)

    with h5py.File(output_hdf5, mode='w') as archive:
        for filename in filenames:
            with h5py.File(filename, mode='r') as f:
                for name in f.keys():
                    f.copy(name, archive)

            # Remove original files
            os.remove(filename)

class StoreSolution:
    """
    Class to save the results.
    """

    def __init__(self, address, filename, format="npz", info=False):
        """
        Initialize the storage class.

        Parameters:
        - address (str): Directory to save files.
        - filename (str): Base filename for saved data.
        - file_format (str): Format to save files ("npz" or "hdf5").
        - info (bool): Print information about file saving.
        """
        if info:
            print(f"Using address: {address}, Filename: {filename}, Format: {format}")

        # Instance attributes
        self.address = address
        self.name = filename
        self.format = format.lower()  # Normalize format to lowercase
        self.time = []

        # Ensure directory exists
        os.makedirs(self.address, exist_ok=True)

    def save_file(self, data, ti):
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
                fname = os.path.join(self.address, f"{self.name}_u_{ti}.h5")
                with h5py.File(fname, "w") as f:
                    f.create_dataset(f"{self.name}_{ti}", data=data)

            # Store time step
            self.time.append(ti)

    def close_file(self, zip_name):
        """
        Combine all .npz files into a single archive.

        Parameters:
        - zip_name (str): Name of the final zip archive (without extension).
        """
        if self.address is not None:
            # Save all the time points where solutions are saved
            if self.format == "npz" and self.name != "grid":
                np.savez(os.path.join(self.address, f"{self.name}_t.dat"), t=np.array(self.time, dtype=float))
            elif self.format == "hdf5" and self.name != "grid":
                with h5py.File(os.path.join(self.address, f"{self.name}_t.h5"), "w") as f:
                    f.create_dataset("t", data=np.array(self.time, dtype=float))

            # Find all files of the selected format
            file_extension = "*.npz" if self.format == "npz" else "*.h5"
            filenames = glob.glob(os.path.join(self.address, f"{self.name}_*{file_extension}"))

            # Archive files
            archive_name = os.path.join(self.address, f"{zip_name}.{self.format}")
            if self.format == "npz":
                JoinFilesInOneZip(filenames, archive_name)
            else:
                JoinFilesInOneHDF5(filenames, archive_name)


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