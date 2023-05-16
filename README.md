<img src="https://github.com/iGW-TU-Dresden/CFPy/blob/main/logo.png" width="400">

# CFPy
CFPy - a Python package for the generation of MODFLOW-CFP specific input files (.cfp, .crch, .coc). CFPy is the FloPy-equivalent package for MODFLOW-CFP.

### Contributing:
If you want to contribute to the package (e.g., changing the code or adding examples), please generally use the `dev` branch for commits (especially for changes to the source code!). New examples can also be directly committed to the `main` branch.

### Dependencies:
- `python` >= 3.9, < 3.11 (without `pyKasso`), == 3.9 (with `pyKasso`)
- `numpy` >= 1.18.5, <1.25.0
- `matplotlib` >= 3.3.4, <3.6.0
- `pandas` >= 1.2.1
- `flopy` >= 3.3.3
- `watermark` >= 2.3.0
- `pykasso` == 0.1.0 (OPTIONAL)

See the beginning of the individual example notebooks for more information on `Python` and package versions.

# Installation
The installation is described specifically for using the [Anaconda distribution](https://www.anaconda.com/products/distribution) of Python / using `conda` environments. If you encounter an OpenSSL-related error during one of the following steps, go to `.../anaconda3/lib/bin`, copy `libcrypto-1_1-x64.dll` and `libssl-1_1-x64.dll` and paste them to `.../anaconda3/dlls`, and continue regularly afterwards.
## Installation of `CFPy` only (without `pyKasso`)
#### For experienced users
Install `CFPy` from source in a (new) environment with Python >= 3.9.

#### For inexperienced users
- **Download** the `CFPy` source code [here](https://github.com/iGW-TU-Dresden/CFPy/tree/main)
    + Make sure that you are in the `main` branch (see upper left of the page)
    + **Download** the files by pressing the green "Code" button in the upper right of the page and selecting "Download ZIP"
    + **Unpack** the ZIP (it should be called `CFPy-main.zip`) on your machine
    + **Get the directory path** where the unpacked `CFPy` source files are now stored (e.g., `C:/Users/.../CFPy-main`) by checking the directory of the file `setup.py` in the `CFPy-main` directory
- **Create a new environment** for `CFPy`
    + **Open Anaconda PowerShell Prompt** on your machine
    + To **create a new environment**, type `conda create -n cfpy_env python=3.9` and press `[Enter]` (if user input is asked - typically for confirming the installation, press `[y]` and `[Enter]` to confirm)
    + **Activate the environment** by typing  `conda activate cfpy_env` and pressing `[Enter]` (and **do NOT close** Anaconda PowerShell Prompt)
- **Install** `CFPy`
    + Remain in Anaconda PowerShell Prompt with the `cfpy_env` environment activated
    + To **install** `CFPy`, type `pip install -e <YourPathToCFPy>` and replace `<YourPathToCFPy>` with your path to `CFPy` source files from earlier, e.g., `C:/Users/.../CFPy-main`, and press `[Enter]` (if user input is asked - typically for confirming the installation, press `[y]` and `[Enter]` to confirm)
    + Optionally, you can install `Jupyter Notebooks` / `Jupyter Lab` (all `CFPy` examples are in this format) in the `cfpy_env` environment by typing `conda install -c conda-forge jupyterlab` and pressing `[Enter]` (if user input is asked - typically for confirming the installation, press `[y]` and `[Enter]` to confirm); make sure that the `cfpy_env` is still activated

You are now ready to use `CFPy`! To use `CFPy` in a `Jupyter Notebook` (typical application) at any time:
- Open Anaconda PowerShell prompt
- Activate the environment via the command `conda activate cfpy_env`
- Open `Jupyter Lab` via the command `jupyter lab`

## Installation of `CFPy` with `pyKasso`
#### For experienced users
- create a new / activate an existing environment environment with Python == 3.9.
- Download `pyKasso` source files from the `cfpy` branch [here](https://github.com/randlab/pyKasso/tree/cfpy)
- Create a new environment from the `environment.yml`
- Install [`karstnet`](https://github.com/karstnet/karstnet) in this environment (from source)
- Install `CFPy` in this environment (from source)

#### For inexperienced users
- **Download** the `pyKasso` source code from [here](https://github.com/randlab/pyKasso/tree/cfpy)
    + Make sure that you are in the `main` branch (see upper left of the page)
    + **Download** the files by pressing the green "Code" button in the upper right of the page and selecting "Download ZIP"
    + **Unpack** the ZIP (it should be called `pyKasso-cfpy.zip`) on your machine
    + **Get the file path** of the `ènvironment.yml` from the unpacked `pyKasso` source directory (e.g., `C:/Users/.../pyKasso-cfpy/environment.yml`
    + **Get the directory path** where the unpacked `pyKasso` source files are now stored (e.g., `C:/Users/.../pyKasso-cfpy`) by checking the directory of the file `setup.py` in the `pyKasso-cfpy` directory
- **Create a new environment** for `pyKasso`
    + **Open Anaconda PowerShell Prompt** on your machine
    + To **create a new environment** from the `environment.yml` file, type `conda create -f <YourPathTo/environment.yml>` (where you replace `<YourPathTo/environment.yml>` with the path you got before) and press `[Enter]` (if user input is asked - typically for confirming the installation, press `[y]` and `[Enter]` to confirm)
    + **Activate the environment** by typing  `conda activate pykasso2D` and pressing `[Enter]` (and **do NOT close** Anaconda PowerShell Prompt)
- **Download and install** `karstnet`
    + **Download** the `karstnet` yource files from [here](https://github.com/karstnet/karstnet) by pressing the green "Code" button in the upper right of the page and selecting "Download ZIP"
    + **Unpack** the ZIP (it should be called `karstnet-master.zip`) on your machine
    + **Get the directory path** where the unpacked `karstnet` source files are now stored (e.g., `C:/Users/.../karstnet-master`) by checking the directory of the file `setup.py` in the `karstnet-master` directory
    + Remain in Anaconda PowerShell Prompt with the `pykasso2D` environment activated
    + To **install** `karstnet`, type `pip install -e <YourPathTokarstnet>` and replace `<YourPathTokarstnet>` with your path to `karstnet` source files from earlier, e.g., `C:/Users/.../karstnet-master`, and press `[Enter]` (if user input is asked - typically for confirming the installation, press `[y]` and `[Enter]` to confirm) and **do NOT close** Anaconda PowerShell Prompt
- **Install** `pyKasso`
    + Remain in Anaconda PowerShell Prompt with the `pykasso2D` environment activated
    + To **install** `pyKasso`, type `pip install -e <YourPathTopyKasso>` and replace `<YourPathTopyKasso>` with your path to `pyKasso` source files from earlier, e.g., `C:/Users/.../pyKasso-cfpy` and press `[Enter]` (if user input is asked - typically for confirming the installation, press `[y]` and `[Enter]` to confirm) and **do NOT close** Anaconda PowerShell Prompt
- **Download and install** `CFPy`
    + **Download** the `CFPy` yource files from [here](https://github.com/iGW-TU-Dresden/CFPy/tree/main) by pressing the green "Code" button in the upper right of the page and selecting "Download ZIP"
    + **Unpack** the ZIP (it should be called `CFPy-main.zip`) on your machine
    + **Get the directory path** where the unpacked `CFPy` source files are now stored (e.g., `C:/Users/.../CFPy-main`) by checking the directory of the file `setup.py` in the `CFPy-main` directory
    + Remain in Anaconda PowerShell Prompt with the `pykasso2D` environment activated
    + To **install** `CFPy`, type `pip install -e <YourPathToCFPy>` and replace `<YourPathToCFPy>` with your path to `CFPy` source files from earlier, e.g., `C:/Users/.../CFPy-main`, and press `[Enter]` (if user input is asked - typically for confirming the installation, press `[y]` and `[Enter]` to confirm)
    + Optionally, you can install `Jupyter Notebooks` / `Jupyter Lab` (all `CFPy` examples are in this format) in the `pykasso2D` environment by typing `conda install -c conda-forge jupyterlab` and pressing `[Enter]` (if user input is asked - typically for confirming the installation, press `[y]` and `[Enter]` to confirm); make sure that the `pykasso2D` is still activated

You are now ready to use `CFPy` together with `pyKasso`! To use `CFPy` with `pyKasso` in a `Jupyter Notebook` (typical application) at any time:
- Open Anaconda PowerShell prompt
- Activate the environment via the command `conda activate pykasso2D`
- Open `Jupyter Lab` via the command `jupyter lab`

# CFPy_Examples:
Note: make sure you have pyKasso (and its dependencies) installed if you want to use the full functionality of CFPy!

## Example 1
- Notebook: EX01_CFPy_FloPy_I.ipynb
- Files: CFPy_EX01_RUN
- Description:
    + simple example for CFP mode 1 from the MODFLOW-CFP documentation, coupling `CFPy` (MODFLOW-CFP input-file generation) and `FloPy` (MODFLOW input-file generation)
    + generating CFP input files from a user-defined node network structure
    + model computation with `cfpv2`

## Example 2
- Notebook: EX02_stochastic_pyKasso_FloPy_I.ipynb
- Files: CFPy_EX02_RUN
- Description:
    + complex example coupling `pyKasso` (network generation), `CFPy` (MODFLOW-CFP input-file generation), and `FloPy` (MODFLOW input-file generation)
    + **multiple** node networks are generated with pyKasso, given a single set of inlet and outlet locations
    + generating CFP input files from the validated networks generated by `pyKasso`
    + model computation with `cfpv2`
    + evaluation / initial assessment of simulation uncertainty, taking multiple network realizations into account

## Example 3
- Notebook: EX03_pyKasso_CFPy_coupling_I.ipynb
- MODFLOW files: none
- Description:
    + demonstration of the coupling process between `pyKasso` and `CFPy`
    + a network is generated with `pyKasso` and validated with the `CFPy.preprocessing` module
    + the structure of the `CFPy` input file (the .nbr-file) is described

## Example 4
- Notebook: EX04_CFPy_FloPy_PreExisting_FloPyModel.ipynb
- Files: CFPy_EX04_RUN
- Description:
    + EX04 is very similar to EX01
    + it is shown how one can use a pre-existing `FloPy` model in `CFPy`
    + first, a simple MODFLOW model is created with `FloPy`
    + passing the `FloPy` model instance to the `CFPy` preprocessor automatically extracts all the usable information from the `FloPy` model (such as layer elevations, discretization etc.)
    + this greatly improves the usability of CFPy for pre-existing MODFLOW models (which do not yet include `CFP`)

## Example 5
- Notebook: EX05_CFPy_FloPy_ConduitWell.ipynb
- Files: CFPy_EX05_RUN
- Description:
    + like example EX01 but with a conduit network flow boundary condition (FBC)

## Example 6
- Notebook: EX06_CFPy_FloPy_II.ipynb
- Files: CFPy_EX06_RUN
- Description:
    + like example EX01 but with an already existing .nbr file (two node planes)

## Example 7
- Notebook: EX07_pyKasso_CFPy_coupling_II.ipynb
- Files: none
- Description:
    + like example EX03 but with some additions (i.e., making the pyKasso-generated network quasi-3D and validating a network with two node planes)

## Example 8
- Notebook: EX08_stochastic_pyKasso_FloPy_II.ipynb
- Files: CFPy_EX08_RUN
- Description:
    + like example EX02 but more streamlined with less additional functions (more straightforward application of CFPy)

## NOTE
To run the examples, make sure that the `CFPy` package directory (the `CFPy` folder) is inside the example directory or that `CFPy` is installed as a package in the active environment.

# Workflow Description

### Workflow Description for Coupling Stochastic Network Generation (e.g., `pyKasso`), `CFPy` and `FloPy` to Create a Script-Based Framework for Spatially Distributed Karst Simulation with MODFLOW-CFP

1. Definition of Model Domain
    - define a rectangular (in the 2D case) or cuboid (in the 3D case) model domain together with a spatial discretization (i.e., number and width of rows, number and width of columns, number of layers together with layer elevation information)
    - if a non-rectangular model domain should be used, still start with a rectangular domain; later, after all input data are generated, the model domain can be altered in FloPy by setting cells inactive

2. `pyKasso` Network Generation
    - include the domain definition from step (1) in the settings.yaml-file for `pyKasso`
    - define other stochastic network generation options in the `settings.yaml`
    - generate one (or multiple) karst networks with `pyKasso`

    - **NOTE**: Right now, `pyKasso` can only generate 2D karst networks. However, multiple (vertically connected) node planes can be given to `CFPy`. To give the 2D network some 3D structure nonetheless, node elevations need to be given to `CFPy`. Those elevations can be uniform or non- uniform. Also, multiple 2D node planes (with each different uniform or non-uniform node elevations) can be connected vertically. See the `pyKasso_CFPy_coupling.ipynb` notebook for more information.

3. Coupling `pyKasso` and `CFPy`
    - **NOTE**: `CFPy` generally relies on one single input-file, where all neccessary information is summarized (domain discretization, number of node planes, MODFLOW layer elevations, node network and elevations). This file can either be generated automatically with the `generate_nbr` method of the `CFPy.preprocessing` module (preferred option) or can also be generated manually
    - use the `pyKassoValidator` object in the `CFPy.preprocessing` module to validate the node network generated with `pyKasso` in step (2); provide node elevation information in this step
    - option 1 (preferred option): generate the `CFPy` input information (the `.nbr`-file) automatically via the `generate_nbr` method in the `CFPy.preprocessing` module (providing the `pyKasso.SKS` catchment as well as domain discretization and layer elevation information)
    - option 2: or generate the `.nbr`-file manually (also see the `pyKasso_CFPy_coupling.ipynb` notebook for more information) by exporting the validated network and copying to the `.nbr`-file

4. Set Up the MODFLOW-CFP Model with `CFPy` and `FloPy`
    - create the model with `FloPy` (using a MODFLOW-CFP distribution, preferrably the `cfpv2.exe` version distributed by TU Dresden)
    - create all input-files for CFP with `CFPy` based on the `.nbr`-file generated in step (3)
    - create the remaining input-files for MODFLOW with `FloPy`
    - simulate the model with `FloPy` (`model.run_model` method)

5. Post-Processing of Results
    - use `FloPy` post-processing methods to process classical MODFLOW results such as matrix head information etc.
    - use the methods in the `CFPy.postprocessing` module to obtain and process CFP-specific results (node- and tube- related data)

**NOTE**: all steps are shown in the example notebooks!
