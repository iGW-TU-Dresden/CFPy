![CFPy - a Python package for the generation of MODFLOW-CFP specific input files (.cfp, .crch, .coc). CFPy is the FloPy-equivalent package for MODFLOW-CFP.]

# CFPy

### Contributing:
If you want to contribute to the package (e.g., changing the code or adding examples), please generally use the `dev` branch for commits (especially for changes to the source code!). New examples can also be directly committed to the `main` branch.

### Dependencies:
- `python` >= 3.9, < 3.11
- `numpy`>= 1.18.5, <1.25.0
- `matplotlib`>= 3.3.4, <3.6.0
- `pykasso`= 0.1.0
- `pandas`>= 1.2.1
- `flopy`>= 3.3.3

See the beginning of the individual example notebooks for more information on `Python` and package versions.

# Installation

You can use / set up a new virtual environment with the environment file (`cfpy_env.yml`) prior to installing `CFPy`, which contains all required dependencies. Note, however, that **`pyKasso` and `karstnet` are not available on `PyPI` and have to be installed from source**! If you have an environment with all the dependencies available, activate it (`conda activate myEnvironment`, replace `myEnvironment` with your environment name). Afterwards you can install `CFPy`.
Download the source code as `.zip` and unpack it on your machine. In the command line (or in the `Anaconda Prompt` / `Anaconda PowerShell Prompt`) navigate to the unpacked folder (e.g., `C:\User\...\CFPy-main`). Install the package via `pip install .` (don't forget the "." at the end).

# CFPy_Examples:
Note: make sure you have pyKasso (and its dependencies) installed if you want to use the full functionality of CFPy!

## Example 1
- Notebook: EX01_CFPy_FloPy.ipynb (updated version) or EX01_CFPy_FloPy_OLD.ipynb (old but functioning version)
- Files: CFPy_EX01_RUN
- Description:
    + simple example for CFP mode 1 from the MODFLOW-CFP documentation, coupling `CFPy` (MODFLOW-CFP input-file generation) and `FloPy` (MODFLOW input-file generation)
    + generating CFP input files from a user-defined node network structure
    + model computation with `cfpv2`

## Example 2
- Notebook: EX02_CFPy_idealized_EASY.ipynb
- Files: CFPy_EX02_RUN
- Description:
    + complex example coupling `pyKasso` (network generation), `CFPy` (MODFLOW-CFP input-file generation), and `FloPy` (MODFLOW input-file generation)
    + **multiple** node networks are generated with pyKasso, given a single set of inlet and outlet locations
    + generating CFP input files from the validated networks generated by `pyKasso`
    + model computation with `cfpv2`
    + evaluation / initial assessment of simulation uncertainty, taking multiple network realizations into account

## Example 3
- Notebook: EX03_pyKasso_CFPy_coupling.ipynb
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
- Notebook: EX05_CFPy_FloPy_ConduitWell
- Files: CFPy_EX05_RUN
- Description:
    + like example EX01 but with a conduit network flow boundary condition (FBC)

## NOTE
To run the examples, make sure that the `CFPy` package directory (the `CFPy` folder) is inside the example directory

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
