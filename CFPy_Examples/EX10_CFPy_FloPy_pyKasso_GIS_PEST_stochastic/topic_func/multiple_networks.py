# This file contains the functions used to create the 2D multiple networks in "02_Stat_model.ipynb"

import glob
import importlib
import os
import shutil

import CFPy as cfpy
import numpy as np
import pykasso as pk


# Check that there are no gaps inside the network
def is_network_connected(matrixx_2, start_points, end_points):
    """Check if all ones in a matrix are linked to each other.

    Returns
    -------
    matrixx_2 :

    start_points : list of tuples
        List of tuples with the start points (row, col)
    end_points : list of tuples
        List of tuples with the end points (row, col)

    """

    rows_2 = len(matrixx_2)
    cols_2 = len(matrixx_2[0])
    visited = [[False] * cols_2 for _ in range(rows_2)]

    def dfs(row2, col2):
        if (
            row2 < 0
            or col2 < 0
            or row2 >= rows_2
            or col2 >= cols_2
            or visited[row2][col2]
            or matrixx_2[row2][col2] == 0
        ):
            return

        visited[row2][col2] = True
        dfs(row2 - 1, col2)  # Up
        dfs(row2 + 1, col2)  # Down
        dfs(row2, col2 - 1)  # Left
        dfs(row2, col2 + 1)  # Right

    for start_point in start_points:
        dfs(start_point[0], start_point[1])

    for end_point in end_points:
        if not visited[end_point[0]][end_point[1]]:
            return False

    return True


def clean_results(path):
    """
    Remove all run-dictionaries and the network plots from the model directory

    Keyword Arguments: -

    Return: -
    """

    # Remove all *jpg files in the model directory
    rm_files = glob.glob(os.path.join(path, "*.jpg"))
    for rmf in rm_files:
        os.remove(rmf)

    rm_dirs = glob.glob(os.path.join(path, "run_*"))
    if os.path.exists("network_plots"):
        rm_dirs.append("network_plots")
    for rmd in rm_dirs:
        shutil.rmtree(rmd)


def create_network(yaml_settings_file, elev_nodes, model):
    """
    Create the network with pyKasso and process/export it with cfpy

    Keyword Arguments: yaml_file -- absolute path to settings file of pyKasso, str

    Return: validator object, valid network array
    """

    # read in settings file
    catchment = pk.SKS(yaml_settings_file=yaml_settings_file)

    catchment.update_all()

    # compute karst networks from the given information
    catchment.compute_karst_network()

    # generate elevation data
    # NOTE: elevation data has to have the same shape as the node network array!
    shp = np.array(catchment.karst_simulations[-1].maps["karst"][-1]).shape

    # Set elevation of all nodes
    elevs = np.ones(shp) * elev_nodes

    # validate the network from pyKasso
    validator = cfpy.preprocessing.pyKassoValidator(
        flopymodel=model, network=catchment, elevations=elevs
    )
    valid_network = validator.validate_network()

    # export the network
    # the exported information can directly be included in the .nbr-file as input for CFPy
    # notes on how to use the generated data with CFPy is given at the end of the notebook
    # per default, the network is exported to "CFPy_exported_network_for_NBR.txt" in the active directory
    # validator.export_network()

    return validator, valid_network


def store_results(number, path, modelname):
    """
    Create a directory and put the current results into the directory.
    Also make a copy of the network plot and put it into the plot directory

    Keyword Arguments:
      number -- identifier of the current iteration, float or int

    Return: -
    """

    # Move results

    # define directory where to store the results
    if number < 10:
        target_dir = os.path.join(path, f"run_00{number}")
    elif number < 100:
        target_dir = os.path.join(path, f"run_0{number}")
    else:
        target_dir = os.path.join(path, f"run_{number}")

    # make sure the directory does not exist
    if os.path.exists(target_dir):
        raise Exception(
            f"Directory {target_dir} exists already! Unable to write results"
        )

    # create directory
    os.makedirs(target_dir)

    # get alle files that will be moved
    files_to_move = (
        glob.glob("NODE*")
        + glob.glob("TUBE*")
        + glob.glob("{}*".format(modelname))
        + glob.glob("*.nbr")
    )

    # + glob.glob("*.coc") + glob.glob("*.cfp") + glob.glob("*.crch") \
    # + glob.glob("*.lpf") + glob.glob("*.bas")

    # move files to directory
    for f in files_to_move:
        source = os.path.join(path, f)
        destination = os.path.join(target_dir, f)
        shutil.move(source, destination)


def check_for_network_duplicates(
    mf, valid_network, networks_collected, count_networks, start_points, end_points
):
    test_network = np.zeros((mf.dis.nrow, mf.dis.ncol)).flatten()

    test_netw = valid_network  # change_greater_than_zero_to_one(valid_network)

    if is_network_connected(test_netw, start_points, end_points):
        print("Network is correctly linked")
    else:
        raise ValueError("Network is not correctly linked")

    network_new = valid_network.flatten()

    if count_networks == 0:
        networks_collected = np.stack((test_network, network_new))

    else:
        networks_collected = np.vstack((networks_collected, network_new))

    for n_network in range(0, len(networks_collected) - 1):
        if np.array_equal(network_new, networks_collected[n_network]):
            raise ValueError("network already exist")

    return networks_collected
