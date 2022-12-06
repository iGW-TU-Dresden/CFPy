# CFPy_TUD
CFPy - a Python package for the generation of MODFLOW-CFP specific input files (.cfp, .crch, .coc). CFPy is the FloPy-equivalent package for MODFLOW-CFP.

### Contributing:
If you want to contribute to the package (e.g., changing the code or adding examples), please generally use the `dev` branch for commits (especially for changes to the source code!). New examples can also be directly committed to the `main` branch.

### Dependencies:
- `numpy`>= 1.18.5
- `matplotlib`>= 3.3.4
- `pykasso`>= 0.1.0
- `re`>= 2.2.1
- `pandas`>= 1.2.1
- `flopy`>= 3.3.3
- `sys`>= 3.8.5

See the beginning of the individual example notebooks for more information on `Python` and package versions.

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
