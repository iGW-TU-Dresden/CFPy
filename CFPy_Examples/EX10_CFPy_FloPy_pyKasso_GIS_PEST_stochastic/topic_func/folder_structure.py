# folder_structure.py v0.1
# Last updated: 2023-12-21

import gzip
import os
import pickle
from functools import wraps

import fiona
import geopandas as gpd
import numpy as np
import openpyxl
import pandas as pd
import rasterio


class Folders_c:

    """
    This is the class to generate the order/structure of the project folders

    Args:
    main_dir (string): Path of the main_dir (i.e. os.getcwd())
    folders (list): List with the names (string) of the main folders
    subfolders (list): List with the names (string) of the subfolders
    """

    def __init__(self, main_dir, folders, subfolders):
        # ------------------------------------------------------------------------------------------
        self.main_dir = main_dir
        self.folders = folders
        self.subfolders = subfolders
        # ------------------------------------------------------------------------------------------

    def create_folders(self):
        """This creates the defined folders in the defined path

        Use this function to create the folders defined in:
        Folders(main_dir = "...", folders = [...], subfolders = [...])
        This will not overwrite or delete existing folders with identical names

        Returns:
            string: Folder structure were generated

        Examples:
            >>> Folders = Folders_c(main_dir = os.getcwd(),
                                    folders = ['01_Data', '02_Model', '03_Results'],
                                    subfolders = [['Boundaries', 'Climate', 'Geology', 'Hydrogeology', 'Hydrology', 'Land_use', 'Topography', 'Water_demand', 'Water_quality'],
                                                  ['FloPy_model', 'FloPy_model_cali', 'GIS', 'Pest'],
                                                  ['shp', 'img', 'pdf', 'csv', 'rst', 'xlsx', 'LaTeX']])

                Folders.create_folders()

            'Folders were generated'

        """

        # ------------------------------------------------------------------------------------------
        for folder in self.folders:
            if not os.path.exists(os.path.join(self.main_dir, folder)):
                os.mkdir(os.path.join(self.main_dir, folder))
        # ------------------------------------------------------------------------------------------
        for i_folder, folder in enumerate(self.folders):
            for subfolder in self.subfolders[i_folder]:
                if not os.path.exists(os.path.join(self.main_dir, folder, subfolder)):
                    os.mkdir(os.path.join(self.main_dir, folder, subfolder))
        # ------------------------------------------------------------------------------------------

        return print("Folders were generated")

    def path_folders(self):
        """This function creates a dictionary with the paths of the defined folders

        Returns
        -------
        paths : dict
            Dictionary with the paths of the defined folders
        """

        # ------------------------------------------------------------------------------------------
        paths = {}
        for folder in self.folders:
            paths[f"{folder}"] = os.path.join(self.main_dir, folder)
        # ------------------------------------------------------------------------------------------
        for i_folder, folder in enumerate(self.folders):
            for subfolder in self.subfolders[i_folder]:
                paths[f"{folder}_{subfolder}"] = os.path.join(
                    self.main_dir, folder, subfolder
                )
        # ------------------------------------------------------------------------------------------

        return paths


class Files_c:

    """
    This is the class listing and reading all files placed in the defined folders

    Args:
    main_dir (string): Path of the main_dir (i.e. os.getcwd())
    folders (list): List with the names (string) of the main folders
    subfolders (list): List with the names (string) of the subfolders
    """

    def __init__(self, main_dir, folders, subfolders):
        # ------------------------------------------------------------------------------------------
        self.main_dir = main_dir
        self.folders = folders
        self.subfolders = subfolders
        # ------------------------------------------------------------------------------------------

    def list_all_files(self):
        """This function puts all files into a list

        Use this function to create a list with all files in the selected path (folder and subfolders) defined in:
        Files_c(main_dir = "...", folders = [...], subfolders = [...])

        Returns:
            list: Filenames as string

        Examples:
            >>> Files_D = Files_c(main_dir = os.getcwd(),
                                  folders = folders = ['01_Data'],
                                  subfolders = [['Boundaries', 'Climate', 'Geology', 'Hydrogeology', 'Hydrology', 'Land_use', 'Topography', 'Water_demand', 'Water_quality']])

                Files_D.list_all_files()

            ['file_1','file_2'...]

        """

        # ------------------------------------------------------------------------------------------
        files = []

        for i_folder, folder in enumerate(self.folders):
            for subfolder in self.subfolders[i_folder]:
                if os.path.exists(os.path.join(self.main_dir, folder, subfolder)):
                    files.append(
                        os.listdir(os.path.join(self.main_dir, folder, subfolder))
                    )

        files = [sublist for sublist in files if sublist]
        files = [file for sublist in files for file in sublist]
        # ------------------------------------------------------------------------------------------

        return files

    def dict_all_files(self):
        """This function puts all files (key) and its path (value) into a dictionary

        Use this function to create a dictionary with all files and their associated path defined in:
        Files_c(main_dir = "...", folders = [...], subfolders = [...])

        Returns:
            dict: Filenames (key) with paths (value) as string

        Examples:
            >>> Files_D = Files_c(main_dir = os.getcwd(),
                                  folders = folders = ['01_Data'],
                                  subfolders = [['Boundaries', 'Climate', 'Geology', 'Hydrogeology', 'Hydrology', 'Land_use', 'Topography', 'Water_demand', 'Water_quality']])

                Files_D.dict_all_files()

            {key ('file_1') : value ('path_1'), ...}

        """

        # ------------------------------------------------------------------------------------------
        files = []
        paths = []
        # ------------------------------------------------------------------------------------------
        for i_folder, folder in enumerate(self.folders):
            for subfolder in self.subfolders[i_folder]:
                if os.path.exists(os.path.join(self.main_dir, folder, subfolder)):
                    folder_files = os.listdir(
                        os.path.join(self.main_dir, folder, subfolder)
                    )

                    files.append(folder_files)
                    paths.append(
                        [os.path.join(self.main_dir, folder, subfolder)]
                        * len(folder_files)
                    )

        # ------------------------------------------------------------------------------------------
        files = [sublist for sublist in files if sublist]
        files = [file for sublist in files for file in sublist]

        paths = [sublist for sublist in paths if sublist]
        paths = [path for sublist in paths for path in sublist]
        # ------------------------------------------------------------------------------------------
        file_paths = []
        for path, file in zip(paths, files):
            file_paths.append(os.path.join(path, file))
        # ------------------------------------------------------------------------------------------
        files_and_paths = dict(zip(files, file_paths))
        # ------------------------------------------------------------------------------------------

        return files_and_paths


def read_all_files(file_path_dict, **kwargs):
    """This function reads all files and puts the information in a dictionary

    Use this function to read all files and store the data within a dictionary
    Filetypes such as:
        (.xlsx, .csv, .shp, .rst, gpkg)
    Args:
        file_path_dict (dict): Dictionary with the filenames (key)(string) and associating path (value)(string) of the files

    Returns:
        dict: Filenames (key) with paths (value) as string

    Examples:
        >>> read_all_files(file_path_dict = Files_D.dict_all_files())

        {key ('file_1'+'...') : value , ...}

    """

    # ------------------------------------------------------------------------------------------
    # Create a dictionary
    data_dict = {}
    # ------------------------------------------------------------------------------------------
    for file, path in file_path_dict.items():
        # ------------------------------------------------------------------------------------------
        if file.endswith(".xlsx"):  # Excel spreadsheet files / .xlsx
            print(file)
            # Read the Excel file into a DataFrame
            df = pd.read_excel(path)
            # Append DataFrames as values
            data_dict[f"{file.replace('.xlsx','')}"] = df
            # ------------------------------------------------------------------------------------------
        elif file.endswith(".shp"):  # Shapefiles / .shp
            print(file)
            df = gpd.read_file(path)
            # Append DataFrames as values
            data_dict[f"{file.replace('.shp','')}"] = df

            for nr in range(len(df.index)):
                data_dict[f"{file.replace('.shp','')}" + f"-{1+nr:03d}"] = df.loc[
                    nr:nr
                ].reset_index(inplace=False)
            # ------------------------------------------------------------------------------------------
        elif file.endswith(".gpkg"):  # Geopackage / .gpkg
            print(file)
            # List all Shapefile layer names

            try:
                shp_layer_names = fiona.listlayers(path)
                for name in shp_layer_names:
                    df = gpd.read_file(path, layer=name)
                    # Append DataFrames as values
                    data_dict[f"{file.replace(file,'')}{name}"] = df
                    for nr in range(len(df.index)):
                        data_dict[
                            f"{file.replace(file,'')}{name}" + f"-{1+nr:03d}"
                        ] = df.loc[nr:nr].reset_index(inplace=False)
            except:
                print(f"--> No shapefile within {file}")

            rst_layer_names = []
            # List all raster layer names
            try:
                with rasterio.open(path) as gpkg:
                    if gpkg.subdatasets:  # Multiple raster
                        for raster in gpkg.subdatasets:
                            with rasterio.open(raster) as src:
                                rst_layer_names.append(src.tags()["IDENTIFIER"])
                    else:  # Just one raster
                        rst_layer_names.append(gpkg.tags()["IDENTIFIER"])
                for name in rst_layer_names:
                    with rasterio.open(path, table=f"{name}") as src:
                        df_data = src.read()  # raster data
                        df_bounds = src.bounds  # raster bounding box
                        df_res = src.res  # raster resolution
                        df_profile = src.profile  # raster profile/description

                        # Append Data to the dictionary as values
                        data_dict[f"{file.replace(file,'')}_{name}_data"] = df_data
                        data_dict[f"{file.replace(file,'')}_{name}_bounds"] = df_bounds
                        data_dict[f"{file.replace(file,'')}_{name}_res"] = df_res
                        data_dict[
                            f"{file.replace(file,'')}_{name}_profile"
                        ] = df_profile
            except:
                print(f"--> No rasterfile within {file}")
            # ------------------------------------------------------------------------------------------
        elif file.endswith(".csv"):  # Comma separated variable / .csv
            print(file)
            df = pd.read_csv(path, sep=",")
            data_dict[f"{file.replace('.csv','')}"] = df
            # ------------------------------------------------------------------------------------------
        elif file.endswith(".tif"):  # Rasterfile / .tif
            print(file)
            with rasterio.open(path) as src:
                df_data = src.read()  # raster data
                df_bounds = src.bounds  # raster bounding box
                df_res = src.res  # raster resolution
                df_profile = src.profile  # raster profile/description

                # Append Data to the dictionary as values
                data_dict[f"{file.replace('.tif','')}_data"] = df_data
                data_dict[f"{file.replace('.tif','')}_bounds"] = df_bounds
                data_dict[f"{file.replace('.tif','')}_res"] = df_res
                data_dict[f"{file.replace('.tif','')}_profile"] = df_profile
                # ------------------------------------------------------------------------------------------

    return data_dict


def save_compress_dict(data, path, name):
    """This function puts the information of a dictionary in a zipped and compressed file format at a specific path

    Use this function to turn the dictionary into a compressed and zipped file (file: xxx.pkl.gz)

    Args:
        data (dict): Dictionary with the data to save
        path (string): Path to save the compresed and zipped file
        name (string): Name of the compresed and zipped file (file: xxx.pkl.gz)
    Returns:
        Saves a "name.pkl.gz" file in the specified path containing the provided data in the dictionary

    """
    # ------------------------------------------------------------------------------------------
    # Your dictionary
    my_dict = data

    # Specify the path for the compressed file
    compressed_file_path = os.path.join(path, name)

    if os.path.exists(compressed_file_path):
        print(f'Compressed file: "{name}" already exists')

        # Prompt the user for input
        user_input = input(f"Type 'yes' to overwrite '{name}': ")
        # Check the user's input
        if user_input.lower() == "yes":
            with gzip.open(compressed_file_path, "wb") as compressed_file:
                pickle.dump(my_dict, compressed_file)
            print(f"Dictionary saved to: {compressed_file_path}")

    else:
        # Save the dictionary to a compressed file
        with gzip.open(compressed_file_path, "wb") as compressed_file:
            pickle.dump(my_dict, compressed_file)
        print(f"Dictionary saved to: {compressed_file_path}")

    # ------------------------------------------------------------------------------------------
    return


def read_compress_dict(path, name):
    """This function reads the information of a zipped and compressed file format at a specific path

    Use this function to read the data of a compressed and zipped file (file: xxx.pkl.gz)

    Args:
        path (string): Path of the the compresed and zipped file
        name (string): Name of the compresed and zipped file
    Returns:
        Reads a "name.pkl.gz" file in the specified path

    """
    # ------------------------------------------------------------------------------------------
    # Specify the path for the compressed file
    compressed_file_path = os.path.join(path, name)

    # Load the dictionary from the compressed file
    with gzip.open(compressed_file_path, "rb") as compressed_file:
        loaded_dict = pickle.load(compressed_file)
    # ------------------------------------------------------------------------------------------
    return loaded_dict
