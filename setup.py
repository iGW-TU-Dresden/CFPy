from setuptools import find_packages, setup

setup(
    name = "CFPy",
    author = "M.G. Rudolph, T. Noffz, L. Grabow, T. Reimann",
    version = "0.1",
    description = "a Python package for the generation of MODFLOW-CFP "
    "specific input files (.cfp, .crch, .coc). CFPy is the "
    "FloPy-equivalent package for MODFLOW-CFP.",
    url = "https://github.com/iGW-TU-Dresden/CFPy",
    author_email = "max_gustav.rudolph@tu-dresden.de, thomas.reimann@tu-dresden.de",
    license = "CC0 1.0 Universal",
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Other Audience",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Hydrology"
        ],
    platforms = "Windows, Linux",
    install_requires = [
        "numpy>=1.18.5, <1.25.0",
        "matplotlib>=3.3.4",
        "pandas>=1.2.1",
        "flopy>=3.3.3"
        ],
    packages = find_packages(exclude=[])
    )
