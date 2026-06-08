import os
import shutil
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pykasso as pk
import CFPy as cfpy
import flopy
import flopy.utils.binaryfile as bf
from scipy.stats import qmc

from .recharge import FlexModelFrac


# ── Default fracture families ─────────────────────────────────────────────────
# Three fracture sets used by pyKasso to generate the karst conduit network.
# Each family has a density (fractures per m²), an orientation range (degrees
# from north), a dip range (degrees), a length range (m), and an alpha exponent
# that controls the shape of the length distribution (power law).

_DEFAULT_FRACTURES = {
    "family_01": {
        "density": 0.000007,
        "orientation": [135, 160],
        "dip": [90, 160],
        "length": [300, 400],
        "alpha": 1.5,
    },
    "family_02": {
        "density": 0.000007,
        "orientation": [45, 80],
        "dip": [90, 120],
        "length": [500, 700],
        "alpha": 1.5,
    },
    "family_03": {
        "density": 0.000004,
        "orientation": [0, 60],
        "dip": [90, 120],
        "length": [500, 700],
        "alpha": 1.5,
    },
}

# Absolute path to the vlab/ directory — used to locate default CSV data files
_VLAB_DIR = Path(__file__).parent


# ── Data ──────────────────────────────────────────────────────────────────────

class Data:
    """
    Container for precipitation and evapotranspiration time series.

    Both series must be pandas Series objects with a DatetimeIndex. The
    expected time step and units depend on the temporal_resolution setting
    of the System that owns this Data object:

    - 'daily'  (default): daily frequency, values in mm/d
    - 'hourly': hourly frequency, values in mm/h

    Attributes
    ----------
    prec : pd.Series
        Precipitation time series, DatetimeIndex.
        Units: mm/d (daily) or mm/h (hourly).
    evap : pd.Series
        Potential evapotranspiration time series, DatetimeIndex.
        Units: mm/d (daily) or mm/h (hourly).
    """

    def __init__(self, prec: pd.Series, evap: pd.Series):
        """
        Parameters
        ----------
        prec : pd.Series
            Precipitation time series (mm/d for daily, mm/h for hourly).
        evap : pd.Series
            Potential evapotranspiration time series (mm/d for daily, mm/h for hourly).
        """
        self.prec = prec
        self.evap = evap

    @classmethod
    def from_csv(cls, prec_path, evap_path, temporal_resolution="daily"):
        """
        Load precipitation and evapotranspiration time series from CSV files.

        Both files must be semicolon-separated with a date column as the first
        column (day-first format, e.g. 01.01.2000). For precipitation, rain
        values are read from the third column (index 2). Precipitation is
        first resampled to daily totals; if temporal_resolution='hourly', the
        daily values are then uniformly disaggregated to hourly by distributing
        each day's total evenly across 24 hours (each hour receives 1/24 of
        the daily amount).

        Parameters
        ----------
        prec_path : str or Path
            Path to the precipitation CSV file.
        evap_path : str or Path
            Path to the evapotranspiration CSV file.
        temporal_resolution : str
            'daily' (default) or 'hourly'. Controls the time step of the
            returned Series and the units (mm/d or mm/h respectively).

        Returns
        -------
        Data
            A new Data instance with prec and evap Series at the requested
            temporal resolution.
        """
        # Read precipitation — only the date column (0) and rain amount column (2)
        prec_raw = pd.read_csv(
            prec_path,
            parse_dates=True,
            index_col=[0],
            sep=";",
            usecols=[0, 2],
            dayfirst=True,
        )

        # Read evapotranspiration — full file, date as index
        evap_raw = pd.read_csv(
            evap_path,
            parse_dates=True,
            index_col=[0],
            sep=";",
            dayfirst=True,
        )

        # Resample to the target time step by summing within each period.
        # For daily mode this aggregates any sub-daily input to daily totals.
        # For hourly mode this aggregates any sub-hourly input to hourly totals,
        # and leaves truly hourly input unchanged.
        if temporal_resolution == "hourly":
            resample_rule = "h"
        else:
            resample_rule = "D"

        prec = prec_raw.resample(resample_rule).sum().squeeze()
        evap = evap_raw.resample(resample_rule).sum().squeeze()

        return cls(prec, evap)

    @classmethod
    def default(cls, temporal_resolution="daily"):
        """
        Load the default climate data bundled with the vlab package.

        Reads prec.csv and evap.csv from the vlab/ directory. The bundled
        files contain daily data. For temporal_resolution='hourly', each
        daily value is disaggregated uniformly to 24 hourly values so that
        summing any 24 consecutive hours recovers the original daily total.

        Parameters
        ----------
        temporal_resolution : str
            'daily' (default) or 'hourly'.

        Returns
        -------
        Data
        """
        # Always load as daily first — that is the native resolution of the
        # bundled CSVs.
        daily = cls.from_csv(
            _VLAB_DIR / "prec.csv",
            _VLAB_DIR / "evap.csv",
            temporal_resolution="daily",
        )

        if temporal_resolution == "hourly":
            # Disaggregate: forward-fill each day's value to all 24 hours,
            # then divide by 24 so each hour carries 1/24 of the daily total.
            prec = daily.prec.resample("h").ffill() / 24
            evap = daily.evap.resample("h").ffill() / 24
            return cls(prec, evap)

        return daily


# ── Epikarst ──────────────────────────────────────────────────────────────────

class Epikarst:
    """
    Soil and epikarst recharge model.

    Wraps FlexModelFrac to simulate the partitioning of precipitation into
    soil storage, slow drainage (matrix recharge), and quick bypass flow
    (conduit recharge) through the epikarst zone.

    To remove the influence of the unknown initial soil moisture state, the
    model is run with a warmup period: the full climate time series is tiled
    n_tiles times and the model is run over the whole tiled sequence. Only
    the last tile — the actual study period — is returned.

    Attributes
    ----------
    n_tiles : int
        Number of times the input time series is repeated. With n_tiles=5,
        the model sees four warmup repetitions before the study period.
    """

    def __init__(self, n_tiles: int = 5):
        """
        Parameters
        ----------
        n_tiles : int
            How many times to tile the climate input. Higher values improve
            the warmup but increase computation time linearly.
        """
        self.n_tiles = n_tiles
        self._model = FlexModelFrac()

    def simulate(self, data: Data, rch_params,
                 temporal_resolution="daily") -> np.ndarray:
        """
        Simulate total recharge (slow + quick) for the study period.

        The FlexModelFrac model operates internally with all fluxes expressed
        as mm/d rates and a time step dt expressed as a fraction of a day.
        This keeps the model parameters (ks, kv, etc.) in mm/d regardless of
        the chosen temporal resolution, and the output is always in mm/d.

        For hourly input (mm/h), values are first multiplied by 24 to convert
        them to mm/d rates, and the model is integrated with dt = 1/24 so that
        each hourly sub-step advances the soil storage by the correct amount.
        The output recharge values are also in mm/d (instantaneous daily-rate
        equivalents at each hourly time step).

        Parameters
        ----------
        data : Data
            Climate input. For temporal_resolution='daily', prec and evap are
            expected in mm/d. For 'hourly', they are expected in mm/h.
        rch_params : array-like, shape (7,)
            Model parameters in order:
            [srmax, lp, ks, gamma, simax, kv, omega].
            All rate parameters (ks, kv) are in mm/d regardless of resolution.
            See RechargeParams for units and meaning.
        temporal_resolution : str
            'daily' (default) or 'hourly'. Controls the time step passed to
            FlexModelFrac and the unit conversion applied to the climate input.

        Returns
        -------
        rch_ts : np.ndarray, shape (len(data.prec),)
            Total recharge (slow + quick) for the study period.
            Always in mm/d (instantaneous daily-rate equivalent per time step).
        """
        # Length of the actual study period (before tiling)
        n = len(data.prec)

        # Tile the climate input to create the warmup repetitions
        prec_tiled = np.tile(data.prec.values.ravel(), self.n_tiles)
        evap_tiled = np.tile(data.evap.values.ravel(), self.n_tiles)

        if temporal_resolution == "hourly":
            # Convert mm/h rates to mm/d so all FlexModelFrac parameters
            # (ks, kv in mm/d) remain in consistent units.
            # dt=1/24 ensures each hourly step integrates 1/24 of a day's
            # worth of flux into the soil storage.
            prec_tiled = prec_tiled * 24
            evap_tiled = evap_tiled * 24
            dt = 1.0 / 24.0
        else:
            dt = 1.0

        # Run the epikarst model over the full tiled input
        rch_slow, rch_quick = self._model.simulate(
            prec=prec_tiled,
            evap=evap_tiled,
            p=np.array(rch_params),
            dt=dt,
        )

        # Add slow and quick components to get total recharge
        rch_total = rch_slow + rch_quick

        # Discard the warmup tiles; return only the last n time steps
        return rch_total[n * (self.n_tiles - 1):]


# ── Parameter groups ──────────────────────────────────────────────────────────

class RechargeParams:
    """
    Parameter group for the epikarst recharge model (FlexModelFrac).

    Holds default values and sampling bounds for the seven recharge model
    parameters. When stochastic=True, parameter values are drawn via Latin
    Hypercube Sampling (LHS) during System.generate_samples(). When
    stochastic=False, the fixed values in self.values are used for every run.

    Class attributes
    ----------------
    NAMES : list of str
        Parameter names in the order expected by FlexModelFrac.simulate().
    _DEFAULTS : dict
        Default fixed values used when stochastic=False.
    _BOUNDS : dict
        Default (low, high) sampling bounds used when stochastic=True.

    Parameters
    ----------
    stochastic : bool
        If True, parameters are sampled via LHS. If False, fixed values are used.
    values : dict, optional
        Override one or more default fixed values. Only used when stochastic=False.
        Example: values={"srmax": 300.} sets soil reservoir capacity to 300 mm
        while keeping all other defaults unchanged.
    bounds : dict, optional
        Override one or more default sampling bounds. Each entry must be a
        (low, high) tuple. Only used when stochastic=True.
    """

    NAMES = ["srmax", "lp", "ks", "gamma", "simax", "kv", "omega"]

    # Default fixed parameter values (used when stochastic=False)
    _DEFAULTS = {
        "srmax":  250.,   # maximum soil reservoir capacity (mm)
        "lp":     0.25,   # fraction of srmax below which ET is limited (-)
        "ks":     1000.,  # slow drainage coefficient from the soil store (mm/d)
        "gamma":  2.,     # non-linearity exponent for slow drainage (-)
        "simax":  1.,     # maximum interception storage capacity (mm)
        "kv":     5.,     # fast bypass drainage coefficient through the epikarst (mm/d)
        "omega":  0.2,    # fraction of soil drainage routed as quick flow to conduits (-)
    }

    # Default LHS sampling bounds (used when stochastic=True)
    _BOUNDS = {
        "srmax":  (100.,  500.),
        "lp":     (0.1,   0.8),
        "ks":     (100.,  5000.),
        "gamma":  (1.,    5.),
        "simax":  (0.5,   10.),
        "kv":     (0.1,   10.),
        "omega":  (0.1,   0.9),
    }

    def __init__(self, stochastic=True, values=None, bounds=None):
        self.stochastic = stochastic

        # Start from the class-level defaults, then overwrite with any user-provided values
        self.values = self._DEFAULTS.copy()
        if values is not None:
            self.values.update(values)

        # Start from the class-level bounds, then overwrite with any user-provided bounds
        self.bounds = self._BOUNDS.copy()
        if bounds is not None:
            self.bounds.update(bounds)

        # Will hold the LHS sample array once System.generate_samples() has been called
        self._samples = None


class NetworkStructure:
    """
    Parameter group for the structural configuration of the karst network.

    Only one parameter is stochastic here: sks_seed, the random seed for the
    karst skeleton generation algorithm. The other properties — number of
    inlets, placement seeds, fracture families — define the physical setup of
    the aquifer and are fixed for a given System instance.

    Class attributes
    ----------------
    NAMES : list of str
        Stochastic parameter names (only sks_seed).
    _DEFAULTS : dict
        Default fixed value for sks_seed.
    _BOUNDS : dict
        Default sampling range for sks_seed when stochastic=True.

    Parameters
    ----------
    stochastic : bool
        If True, sks_seed is drawn from LHS. If False, the default seed is used.
    values : dict, optional
        Override the default sks_seed value, e.g. values={"sks_seed": 999}.
    bounds : dict, optional
        Override the sampling range, e.g. bounds={"sks_seed": (0, 50000)}.
    n_inlets : int
        Number of sinkhole inlet points generated by pyKasso.
    inlets_seed : int
        Random seed controlling the spatial placement of the sinkhole inlets.
    outlets_seed : int
        Random seed controlling the placement of the spring outlet.
    fracture_families : dict, optional
        pyKasso fracture family definitions. If None, the built-in default
        fracture families (_DEFAULT_FRACTURES) are used.
    """

    NAMES = ["sks_seed"]

    # Default fixed value for the karst skeleton random seed
    _DEFAULTS = {"sks_seed": 12345}

    # Default LHS sampling range for sks_seed
    _BOUNDS = {"sks_seed": (0, 99999)}

    def __init__(self, stochastic=True, values=None, bounds=None,
                 n_inlets=5, inlets_seed=456, outlets_seed=7,
                 fracture_families=None):
        self.stochastic = stochastic

        # Start from the class-level defaults, then apply any user overrides
        self.values = self._DEFAULTS.copy()
        if values is not None:
            self.values.update(values)

        self.bounds = self._BOUNDS.copy()
        if bounds is not None:
            self.bounds.update(bounds)

        # Will hold the LHS sample array once System.generate_samples() has been called
        self._samples = None

        # Fixed structural properties — these define the aquifer geometry and
        # are not varied stochastically
        self.n_inlets    = n_inlets
        self.inlets_seed = inlets_seed
        self.outlets_seed = outlets_seed

        # Use the provided fracture families, or fall back to the module defaults
        if fracture_families is not None:
            self.fracture_families = fracture_families
        else:
            self.fracture_families = _DEFAULT_FRACTURES


class NetworkProperties:
    """
    Parameter group for hydraulic properties of the conduit network.

    These parameters control how water flows through the conduit pipes
    (diameter, wall roughness) and how the conduit exchanges water with the
    surrounding porous matrix (k_exchange, cad). The crch_fraction parameter
    controls what fraction of the recharge at conduit inlet cells is routed
    into the conduit rather than the matrix.

    Class attributes
    ----------------
    NAMES : list of str
        Parameter names.
    _DEFAULTS : dict
        Default fixed values.
    _BOUNDS : dict
        Default LHS sampling bounds.

    Parameters
    ----------
    stochastic : bool
        If True, parameters are sampled via LHS. If False, fixed values are used.
    values : dict, optional
        Override one or more default fixed values.
        Keys: k_exchange, diameter, rheight, cad, crch_fraction.
    bounds : dict, optional
        Override one or more default sampling bounds.

    Notes
    -----
    When stochastic=True (the default), crch_fraction is included in LHS
    sampling with bounds (0.0, 1.0). To fix it at a specific value while
    keeping other parameters stochastic, pass either:
        bounds={"crch_fraction": (1.0, 1.0)}   # fix at 1.0
    or construct with stochastic=False and set values={"crch_fraction": 1.0}.
    """

    NAMES = ["k_exchange", "diameter", "rheight", "cad", "crch_fraction"]

    # Default fixed parameter values
    _DEFAULTS = {
        "k_exchange":    1e-4,  # conduit–matrix exchange coefficient (m/s)
        "diameter":      0.5,   # conduit pipe diameter, applied to all pipes (m)
        "rheight":       0.02,  # conduit wall roughness height (m)
        "cad":           0.01,  # conduit storage coefficient (-)
        "crch_fraction": 1.0,   # fraction of inlet-cell recharge entering the conduit (-)
    }

    # Default LHS sampling bounds
    _BOUNDS = {
        "k_exchange":    (1e-5, 1e-4),
        "diameter":      (0.1,  1.0),
        "rheight":       (0.001, 0.1),
        "cad":           (0.01, 0.5),
        "crch_fraction": (0.0, 1.0),
    }

    def __init__(self, stochastic=True, values=None, bounds=None):
        self.stochastic = stochastic

        self.values = self._DEFAULTS.copy()
        if values is not None:
            self.values.update(values)

        self.bounds = self._BOUNDS.copy()
        if bounds is not None:
            self.bounds.update(bounds)

        self._samples = None


class AquiferParams:
    """
    Parameter group for matrix aquifer hydraulic properties.

    These parameters describe the porous rock matrix surrounding the conduit
    network and are used in the MODFLOW LPF (Layer-Property Flow) package.

    Class attributes
    ----------------
    NAMES : list of str
        Parameter names.
    _DEFAULTS : dict
        Default fixed values.
    _BOUNDS : dict
        Default LHS sampling bounds.

    Parameters
    ----------
    stochastic : bool
        If True, parameters are sampled via LHS. If False, fixed values are used.
    values : dict, optional
        Override one or more default fixed values. Keys: hk, ss, sy.
    bounds : dict, optional
        Override one or more default sampling bounds.
    """

    NAMES = ["hk", "ss", "sy"]

    # Default fixed parameter values
    _DEFAULTS = {
        "hk": 1e-4,  # horizontal hydraulic conductivity of the rock matrix (m/s)
        "ss": 1e-4,  # specific storage (1/m) — elastic storage under confined conditions
        "sy": 0.01,  # specific yield (-) — drainable porosity under unconfined conditions
    }

    # Default LHS sampling bounds
    _BOUNDS = {
        "hk": (1e-5, 1e-4),
        "ss": (1e-5, 1e-3),
        "sy": (0.01, 0.3),
    }

    def __init__(self, stochastic=True, values=None, bounds=None):
        self.stochastic = stochastic

        self.values = self._DEFAULTS.copy()
        if values is not None:
            self.values.update(values)

        self.bounds = self._BOUNDS.copy()
        if bounds is not None:
            self.bounds.update(bounds)

        self._samples = None


# ── Network ───────────────────────────────────────────────────────────────────

class Network:
    """
    Container for the outputs of a completed pyKasso network generation run.

    Passing a Network object to Groundwater.run_modflow() (or
    System.simulate_from_network()) allows the MODFLOW + CFP simulation to be
    run multiple times on the same network geometry, without repeating the
    computationally expensive network generation step.

    Attributes
    ----------
    validator : cfpy.preprocessing.GeneralValidator
        Validated network object from CFPy. Contains the 2-D conduit network
        map and elevation array, and exposes methods to write them to the
        MODFLOW-readable formats expected by CFP.
    idxs_spring : list of int
        1-based MODFLOW [col, row, lay] index of the spring outlet node.
    idxs_inlet : list of list of int
        List of 1-based MODFLOW [col, row, lay] indices, one per sinkhole inlet.
    """

    def __init__(self, validator, idxs_spring, idxs_inlet):
        """
        Parameters
        ----------
        validator : GeneralValidator
            CFPy validator instance with validated network and elevation data.
        idxs_spring : list of int
            1-based [col, row, lay] MODFLOW index of the spring node.
        idxs_inlet : list of list of int
            1-based [col, row, lay] MODFLOW indices for each inlet node.
        """
        self.validator    = validator
        self.idxs_spring  = idxs_spring
        self.idxs_inlet   = idxs_inlet


# ── Groundwater ───────────────────────────────────────────────────────────────

class Groundwater:
    """
    Karst network generation and MODFLOW + CFP groundwater simulation.

    This class holds all fixed solver settings and structural parameters for
    the pyKasso network generation and the MODFLOW-CFP simulation. It does
    not store any of the stochastic parameters — those are passed at runtime
    through a params dict from System._get_params().

    The main public methods are:
        simulate()         — run the full pipeline (network + MODFLOW + CFP)
        generate_network() — run only the pyKasso network generation
        run_modflow()      — run only the MODFLOW + CFP step on an existing network

    Parameters
    ----------
    exe_name : str
        Full path to the CFPv2.exe executable.
    nz : int
        Number of vertical layers in the pyKasso 3-D grid. This is the
        internal pyKasso resolution, not the MODFLOW layer count.
    vka_ratio : float
        Ratio of vertical to horizontal hydraulic conductivity (Kv / Kh).
        A value of 1.0 means isotropic conditions.
    H_init : float
        Initial hydraulic head (m) assigned to all MODFLOW cells at the
        start of the simulation.
    nstp : int
        Number of time steps in each stress period. By default, each stress
        period has a length of one day.
    mftol : float
        Head convergence criterion for the MODFLOW PCG solver (m). The
        solver stops iterating when the maximum head change between
        iterations falls below this value.
    mfrelax : float
        Relaxation factor for the MODFLOW PCG solver (dimensionless, 0–1).
        Values close to 1 converge faster but may oscillate for difficult
        problems.
    mxiter : int
        Maximum number of outer iterations for the MODFLOW PCG solver.
    tortuosity : float
        Conduit tortuosity factor (-). Ratio of the actual conduit path
        length to the straight-line distance between nodes. A value of 1.0
        means a perfectly straight conduit.
    lcritrey : float
        Lower critical Reynolds number for the transition from laminar to
        turbulent conduit flow.
    hcritrey : float
        Upper critical Reynolds number for fully turbulent conduit flow.
    cfptol : float
        Convergence tolerance for the CFP conduit flow solver.
    cfprelax : float
        Relaxation factor for the CFP solver (0–1).
    perlen : float
        Length of each stress period in model time units. Default is 86400
        seconds (= 1 day), matching the daily recharge time series.
    time_unit : int
        MODFLOW time unit flag. 1 = seconds.
    """

    def __init__(
        self,
        exe_name="C:/WRDAPP/CFPv2.exe",
        nz=10,
        vka_ratio=1.0,
        H_init=75.0,
        nstp=1,
        mftol=1e-8,
        mfrelax=0.98,
        mxiter=5000,
        tortuosity=1.0,
        lcritrey=500,
        hcritrey=5000,
        cfptol=1e-9,
        cfprelax=0.98,
        perlen=86400,
        time_unit=1,
    ):
        self.exe_name    = exe_name
        self.nz          = nz
        self.vka_ratio   = vka_ratio
        self.H_init      = H_init
        self.nstp        = nstp
        self.mftol       = mftol
        self.mfrelax     = mfrelax
        self.mxiter      = mxiter
        self.tortuosity  = tortuosity
        self.lcritrey    = lcritrey
        self.hcritrey    = hcritrey
        self.cfptol      = cfptol
        self.cfprelax    = cfprelax
        self.perlen      = perlen
        self.time_unit   = time_unit

    def plot_network(self, network, inlet_indices, spring_index, n_rows, n_cols):
        """
        Plot the 2-D conduit network map with inlet and spring locations.

        Conduit node elevations are shown as a colour map. Background cells
        (no conduit present) are shown as white. Sinkhole inlets are marked
        with red downward triangles; the spring outlet is marked with a cyan
        right-pointing triangle.

        Parameters
        ----------
        network : np.ndarray, shape (n_rows, n_cols)
            2-D conduit network elevation map. Cells with no conduit have
            a value of 0.
        inlet_indices : list of list of int
            1-based [col, row, lay] indices for each sinkhole inlet.
        spring_index : list of int
            1-based [col, row, lay] index of the spring outlet node.
        n_rows : int
            Number of grid rows.
        n_cols : int
            Number of grid columns.
        """
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 5))

        # Replace zeros (no conduit) with NaN so they appear white on the plot
        network_plot = network.copy()
        network_plot[network_plot == 0.] = np.nan

        # Plot conduit node elevations as a colour map
        im = ax.matshow(network_plot, cmap="cividis", vmin=2., vmax=11.)
        plt.colorbar(im, label="Elevation [m]", shrink=0.7)

        # Plot inlet locations — convert 1-based col/row to 0-based for matshow
        inlets_array = np.array(inlet_indices)
        ax.scatter(
            inlets_array[:, 0] - 1, inlets_array[:, 1] - 1,
            c="r", marker="v", edgecolor="k", s=70,
        )

        # Plot spring location (col and row are already in matshow coordinates)
        ax.scatter(
            [spring_index[0]], [spring_index[1]],
            c="cyan", marker=">", edgecolor="k", s=70,
            clip_on=False, zorder=1000,
        )

        ax.set_xlabel("Col")
        ax.set_ylabel("Row")
        ax.set_xlim(-0.5, n_cols - 0.5)
        ax.set_ylim(-0.5, n_rows - 0.5)
        plt.show()

    def simulate(self, system, recharge_ts, ss_rch, params: dict, proj_name, working_dir,
                 n_stress_periods=None, plot_network=False, routing_mode="both"):
        """
        Run the full pipeline: pyKasso network generation + MODFLOW + CFP.

        This is a convenience method that calls _generate_network() followed
        by _run_modflow(). Use generate_network() and run_modflow() separately
        if you want to reuse the same network geometry for multiple MODFLOW runs.

        Parameters
        ----------
        system : System
            The parent System object (provides grid dimensions, elevations, etc.).
        recharge_ts : np.ndarray
            Recharge time series (mm/d). One value per transient stress period.
        ss_rch : float
            The steady-state recharge flux (mm/d). If None, the mean value
            of the recharge time series is used instead.
        params : dict
            Full parameter dict from System._get_params().
        proj_name : str
            Base name for output directories and files.
        working_dir : str
            Path to the directory where model files will be written.
        n_stress_periods : int or None
            If given, truncate recharge_ts to this many time steps.
        plot_network : bool
            If True, show a plot of the generated conduit network.

        Returns
        -------
        success : bool
            True if MODFLOW converged successfully.
        Q_out : pd.Series
            Spring discharge time series (m³/s).
        heads : np.ndarray, shape (n_pers, n_lays, n_rows, n_cols)
            Simulated hydraulic head for each stress period and layer.
        """
        # Truncate the recharge time series if requested
        if n_stress_periods is not None:
            recharge_ts = recharge_ts[:n_stress_periods]

        # Step 1: generate the karst conduit network with pyKasso
        network = self._generate_network(
            proj_name=proj_name,
            sks_seed=int(params["sks_seed"]),
            n_inlets=int(params["n_inlets"]),
            inlets_seed=int(params["inlets_seed"]),
            outlets_seed=int(params["outlets_seed"]),
            fracture_families=params["fracture_families"],
            n_rows=system.n_rows,
            n_cols=system.n_cols,
            delr=system.delr,
            delc=system.delc,
            lay_elevs=system.lay_elevs,
            chb_spring=system.chb_spring,
            working_dir=str(working_dir),
            plot_network=plot_network,
        )

        # Step 2: run MODFLOW + CFP on the generated network
        return self._run_modflow(
            network=network,
            proj_name=proj_name,
            recharge_ts=recharge_ts,
            ss_rch=ss_rch,
            hk=float(params["hk"]),
            ss=float(params["ss"]),
            sy=float(params["sy"]),
            k_exchange=float(params["k_exchange"]),
            diameter=float(params["diameter"]),
            rheight=float(params["rheight"]),
            cad=float(params["cad"]),
            n_rows=system.n_rows,
            n_cols=system.n_cols,
            n_lays=system.n_lays,
            delr=system.delr,
            delc=system.delc,
            lay_elevs=system.lay_elevs,
            chb_spring=system.chb_spring,
            working_dir=str(working_dir),
            routing_mode=routing_mode,
            crch_fraction=float(params.get("crch_fraction", 1.0)),
        )

    def generate_network(self, system, params: dict, proj_name, working_dir,
                         plot_network=False) -> "Network":
        """
        Run only the pyKasso network generation phase.

        Use this when you want to generate a network once and then run
        multiple MODFLOW simulations on it with different hydraulic parameters.
        Pass the returned Network object to run_modflow().

        Parameters
        ----------
        system : System
            The parent System object (provides grid dimensions, etc.).
        params : dict
            Full parameter dict from System._get_params().
        proj_name : str
            Base name for the pyKasso output directory.
        working_dir : str
            Path to the directory where the pyKasso project will be created.
        plot_network : bool
            If True, show a plot of the generated network.

        Returns
        -------
        Network
            Validated network object. Pass to run_modflow() or
            System.simulate_from_network().
        """
        return self._generate_network(
            proj_name=proj_name,
            sks_seed=int(params["sks_seed"]),
            n_inlets=int(params["n_inlets"]),
            inlets_seed=int(params["inlets_seed"]),
            outlets_seed=int(params["outlets_seed"]),
            fracture_families=params["fracture_families"],
            n_rows=system.n_rows,
            n_cols=system.n_cols,
            delr=system.delr,
            delc=system.delc,
            lay_elevs=system.lay_elevs,
            chb_spring=system.chb_spring,
            working_dir=str(working_dir),
            plot_network=plot_network,
        )

    def run_modflow(self, system, network: "Network", recharge_ts, ss_rch, params: dict,
                    proj_name, working_dir, n_stress_periods=None,
                    routing_mode="both"):
        """
        Run MODFLOW + CFP given a pre-generated Network object.

        Parameters
        ----------
        system : System
            The parent System object (provides grid dimensions, etc.).
        network : Network
            Output of generate_network() or System.generate_network().
        recharge_ts : np.ndarray
            Recharge time series (mm/d). One value per transient stress period.
        params : dict
            Full parameter dict from System._get_params(). Must contain:
            hk, ss, sy, k_exchange, diameter, rheight, cad, crch_fraction.
        proj_name : str
            Base name for the MODFLOW output directory.
        working_dir : str
            Path to the directory where MODFLOW files will be written.
        n_stress_periods : int or None
            If given, truncate recharge_ts to this many time steps.
        routing_mode : str
            Controls where recharge is applied. See System docstring for details.
            'both' (default), 'matrix_only', or 'conduit_only'.

        Returns
        -------
        success : bool
            True if MODFLOW converged.
        Q_out : pd.Series
            Spring discharge time series (m³/s).
        heads : np.ndarray, shape (n_pers, n_lays, n_rows, n_cols)
            Simulated hydraulic head for each stress period and layer.
        """
        # Truncate the recharge time series if requested
        if n_stress_periods is not None:
            recharge_ts = recharge_ts[:n_stress_periods]

        return self._run_modflow(
            network=network,
            proj_name=proj_name,
            ss_rch=ss_rch,
            recharge_ts=recharge_ts,
            hk=float(params["hk"]),
            ss=float(params["ss"]),
            sy=float(params["sy"]),
            k_exchange=float(params["k_exchange"]),
            diameter=float(params["diameter"]),
            rheight=float(params["rheight"]),
            cad=float(params["cad"]),
            n_rows=system.n_rows,
            n_cols=system.n_cols,
            n_lays=system.n_lays,
            delr=system.delr,
            delc=system.delc,
            lay_elevs=system.lay_elevs,
            chb_spring=system.chb_spring,
            working_dir=str(working_dir),
            routing_mode=routing_mode,
            crch_fraction=float(params.get("crch_fraction", 1.0)),
        )

    def _generate_network(self, proj_name, sks_seed, n_inlets, inlets_seed,
                          outlets_seed, fracture_families, n_rows, n_cols,
                          delr, delc, lay_elevs, chb_spring, working_dir,
                          plot_network=False) -> "Network":
        """
        Internal method: run pyKasso to generate and validate a karst network.

        Generates a stochastic 3-D karst conduit network using the Riemann
        skeletonisation algorithm. The 3-D network is then collapsed to 2-D
        by taking the maximum elevation in each vertical column. The 2-D
        result is validated and converted to MODFLOW cell indices.

        pyKasso requires the working directory to be set before creating a
        project. The caller (System.generate_network) is responsible for
        restoring the original working directory afterwards.

        Parameters
        ----------
        proj_name : str
            Base name for the pyKasso project directory.
        sks_seed : int
            Random seed for the karst skeleton (sks) generation algorithm.
            Changing this produces a different network geometry.
        n_inlets : int
            Number of sinkhole inlet points on the land surface.
        inlets_seed : int
            Random seed controlling the spatial placement of the inlets.
        outlets_seed : int
            Random seed controlling the placement of the spring outlet.
        fracture_families : dict
            pyKasso fracture family definitions (density, orientation, etc.).
        n_rows : int
            Number of grid rows (north–south, MODFLOW convention).
        n_cols : int
            Number of grid columns (east–west).
        delr : float
            Cell size in the column (east–west) direction (m).
        delc : float
            Cell size in the row (north–south) direction (m).
        lay_elevs : list of float
            [top_elevation, bottom_elevation] of the aquifer layer (m).
        chb_spring : float
            Fixed hydraulic head at the spring outlet (m). Used as the
            water table elevation in the pyKasso domain definition.
        working_dir : str
            Directory where the pyKasso project folder will be created.
        plot_network : bool
            If True, show a map of the generated conduit network.

        Returns
        -------
        Network
            Validated network with spring and inlet MODFLOW cell indices.
        """
        # pyKasso requires the working directory to be set to the project directory
        os.chdir(working_dir)

        # --- Set up the pyKasso 3-D grid ---
        # The vertical extent of the aquifer is divided equally across nz layers
        dz = (lay_elevs[0] - lay_elevs[-1]) / self.nz
        grid_parameters = {
            "x0": 0, "y0": 0, "z0": 0,
            "nx": n_cols, "ny": n_rows, "nz": self.nz,
            "dx": delr,   "dy": delc,   "dz": dz,
        }

        # --- Initialise pyKasso and create a new project ---
        app = pk.pykasso()
        app.new_project(name=proj_name + "_pyKasso", grid_parameters=grid_parameters)
        # Disable interactive visualiser output (not needed in scripted/batch runs)
        app.visualizer.notebook = True

        # --- Define all model generation parameters ---
        model_parameters = {
            "sks": {
                "seed": sks_seed,
                "algorithm": "Riemann3",  # 3-D Riemann skeletonisation
                "ratio": 0.1,             # fraction of domain nodes used for the skeleton
            },
            "domain": {
                # Bedrock surface: flat at z=0 (bottom of the domain)
                "bedrock": np.zeros((n_cols, n_rows)),
                # Water table: flat at the spring head elevation
                "water_table": np.full((n_cols, n_rows), chb_spring),
            },
            "outlets": {
                "seed": outlets_seed,
                "number": 1,
                "subdomain": "domain_borders_bottom",  # spring is at the aquifer boundary
            },
            "inlets": {
                "number": n_inlets,
                "seed": inlets_seed,
                "subdomain": "domain_surface",  # sinkholes are at the land surface
            },
            "fractures": {
                "generate": fracture_families,
            },
        }

        # Fix the numpy random seed for full reproducibility of the pyKasso run
        np.random.seed(sks_seed)
        app.model.generate(model_parameters=model_parameters)

        # --- Extract 3-D arrays from pyKasso and reorder axes ---
        # pyKasso stores arrays in (nx, ny, nz) = (n_cols, n_rows, nz) order.
        # np.moveaxis(..., 1, 0) swaps the first two axes to give (n_rows, n_cols, nz),
        # which matches the MODFLOW row-major convention.
        network_3d    = np.moveaxis(np.array(app.model.maps["karst"][0]), 1, 0)
        elevations_3d = np.moveaxis(app.model.node_elev_arr, 1, 0)
        outlets_3d    = np.moveaxis(app.model.maps["outlets"], 1, 0)

        # --- Collapse 3-D arrays to 2-D ---
        # Take the maximum across the vertical (nz) axis to project onto a 2-D plan view.
        # The [::-1, :] flip reverses the row axis to reconcile the north-south orientation
        # between pyKasso (row 0 = south) and MODFLOW (row 0 = north).
        network_2d    = np.max(network_3d,    axis=-1)[::-1, :]
        elevations_2d = np.max(elevations_3d, axis=-1)[::-1, :]
        outlets_2d    = np.nanmax(outlets_3d, axis=-1)[::-1, :]

        # --- Validate the network using CFPy ---
        # GeneralValidator checks that conduit node elevations are consistent
        # and prepares the network for export to MODFLOW-readable files
        validator = cfpy.preprocessing.GeneralValidator(
            network=network_2d, elevations=elevations_2d
        )
        validator.validate_network_elevations()

        # --- Find the spring outlet MODFLOW cell index ---
        # argwhere returns positions [row, col] of non-NaN cells in the outlet map
        spring_position = np.argwhere(~np.isnan(outlets_2d))[0]
        # Convert 0-based array indices to 1-based MODFLOW [col, row, lay] indices
        idxs_spring = [
            int(spring_position[1] + 1),  # col: 0-based → 1-based
            int(spring_position[0] + 1),  # row: 0-based → 1-based
            1,                             # lay: single-layer model
        ]

        # --- Find all inlet (sinkhole) MODFLOW cell indices ---
        # Extract the node coordinate table from pyKasso as a DataFrame
        node_df = pd.DataFrame.from_dict(
            app.model.vectors["nodes"],
            orient="index",
            columns=["x", "y", "z", "type"],
        )

        # Keep only the rows that pyKasso classified as inlet (sinkhole) nodes
        inlet_df = node_df[node_df["type"] == "inlet"]

        idxs_inlet = []
        for _, row in inlet_df.iterrows():
            # Convert spatial coordinates (m) to grid indices using pyKasso's grid
            col_idx, row_idx = app.model.grid.get_indices((row["x"], row["y"]))
            # Convert to 1-based MODFLOW [col, row, lay] indices.
            # Note: pyKasso rows increase southward, so row must be flipped.
            modflow_col = int(col_idx[0]) + 1
            modflow_row = n_rows - int(row_idx[0])
            idxs_inlet.append([modflow_col, modflow_row, 1])

        # --- Optionally show a map of the network ---
        if plot_network:
            self.plot_network(
                validator.network, idxs_inlet, idxs_spring, n_rows, n_cols
            )

        return Network(validator, idxs_spring, idxs_inlet)

    def _run_modflow(self, network: "Network", proj_name, ss_rch, recharge_ts, hk, ss, sy,
                     k_exchange, diameter, rheight, cad, n_rows, n_cols, n_lays,
                     delr, delc, lay_elevs, chb_spring, working_dir,
                     routing_mode="both", crch_fraction=1.0):
        """
        Internal method: set up and run a MODFLOW-CFP simulation.

        Builds all required MODFLOW packages (DIS, BAS, LPF, OC, PCG, RCH),
        writes the CFP input files (cfp, coc, crch), runs the coupled model
        via CFPv2.exe, and reads back hydraulic heads and spring discharge.

        The simulation always starts with one steady-state stress period (using
        mean recharge) to spin up the head field. All subsequent stress periods
        are transient, one per day of the recharge time series.

        Parameters
        ----------
        network : Network
            Validated conduit network from _generate_network().
        proj_name : str
            Base name for the MODFLOW output directory.
        recharge_ts : np.ndarray
            Recharge time series (mm/d). Length sets the number of transient
            stress periods.
        ss_rch : float
            The steady-state recharge flux (mm/d). If None, the mean value
            of the recharge time series is used instead.
        hk : float
            Horizontal hydraulic conductivity of the rock matrix (m/s).
        ss : float
            Specific storage of the matrix (1/m).
        sy : float
            Specific yield of the matrix (-).
        k_exchange : float
            Conduit–matrix exchange coefficient, applied to all nodes (m/s).
        diameter : float
            Conduit pipe diameter, applied to all pipes (m).
        rheight : float
            Conduit wall roughness height (m).
        cad : float
            Conduit storage coefficient (-).
        n_rows : int
            Number of grid rows.
        n_cols : int
            Number of grid columns.
        n_lays : int
            Number of model layers.
        delr : float
            Cell size in the column direction (m).
        delc : float
            Cell size in the row direction (m).
        lay_elevs : list of float
            [top_elevation, bottom_elevation] of the aquifer layer (m).
        chb_spring : float
            Fixed hydraulic head assigned to the spring outlet node (m).
        working_dir : str
            Directory where MODFLOW files will be written.

        Returns
        -------
        success : bool
            True if MODFLOW converged.
        Q_out : pd.Series
            Spring discharge time series (m³/s).
        heads : np.ndarray, shape (n_pers, n_lays, n_rows, n_cols)
            Simulated hydraulic head for each stress period and layer.
        """
        # --- Create the MODFLOW output directory and change into it ---
        modflow_dir = os.path.join(working_dir, proj_name + "_modflow")
        os.makedirs(modflow_dir, exist_ok=True)
        os.chdir(modflow_dir)

        modelname = "pyKasso_example"

        # --- Build layer elevation arrays (FloPy needs 2-D arrays) ---
        # FloPy's DIS package and GeneralValidator.generate_nbr() expect top and
        # bottom as (n_rows × n_cols) arrays, not scalars
        top_array = np.ones((n_rows, n_cols)) * lay_elevs[0]
        bot_array = np.ones((n_rows, n_cols)) * lay_elevs[1]
        lay_elevs_array = [top_array, bot_array]

        # --- Initialise the MODFLOW model object ---
        mf = flopy.modflow.Modflow(modelname, exe_name=self.exe_name)

        # Total stress periods = 1 steady-state warm-up + 1 per recharge time step
        n_pers = len(recharge_ts) + 1

        # --- DIS package: spatial and temporal discretisation ---
        # nstp: one internal time step per stress period
        # steady: period 0 is steady-state (1), all others are transient (0)
        nstp_list   = [1] + [self.nstp] * (n_pers - 1)
        steady_list = [1] + [0] * (n_pers - 1)

        flopy.modflow.ModflowDis(
            mf, n_lays, n_rows, n_cols, n_pers,
            delr, delc,
            top=lay_elevs[0],
            botm=lay_elevs[1],
            perlen=self.perlen,
            nstp=nstp_list,
            steady=steady_list,
            itmuni=self.time_unit,
            lenuni=2,  # length units: metres
        )

        # --- BAS package: boundary conditions and initial head ---
        # ibound = 1 everywhere means all cells are active (no fixed-head or
        # no-flow boundaries in the matrix; the spring CHB is handled by CFP)
        ibound_array = np.ones((n_lays, n_rows, n_cols), dtype=np.int32)
        flopy.modflow.ModflowBas(mf, ibound=ibound_array, strt=self.H_init)

        # --- LPF package: layer-property flow (hydraulic conductivity and storage) ---
        # laytyp=1: unconfined (convertible) layer
        # layvka=0: vka is interpreted as the ratio Kv/Kh (not absolute Kv)
        flopy.modflow.ModflowLpf(
            mf,
            laytyp=1, chani=1, layvka=0, laywet=0,
            ipakcb=50, hdry=999.,
            hk=hk, hani=1,
            vka=hk * self.vka_ratio,
            wetdry=0, ss=ss, sy=sy,
        )

        # --- OC package: output control ---
        # Save heads and print the water budget for every stress period
        oc_stress_period_data = {}
        for i in range(n_pers):
            for j in range(nstp_list[i]):
                oc_stress_period_data[(i, j)] = ["save head", "print budget"]

        flopy.modflow.ModflowOc(mf, stress_period_data=oc_stress_period_data)

        # --- PCG package: preconditioned conjugate-gradient solver ---
        flopy.modflow.ModflowPcg(
            mf,
            mxiter=self.mxiter, iter1=self.mxiter, npcond=1,
            hclose=self.mftol, rclose=self.mftol, relax=self.mfrelax,
            nbpol=2, iprpcg=5, mutpcg=0, damp=0.99, ihcofadd=9999,
        )

        # --- RCH package: distributed recharge ---
        # How recharge is spatially distributed depends on routing_mode:
        #
        #   'both' and 'matrix_only':
        #       A single scalar rate is applied uniformly to the top active cell
        #       of every column (nrchop=1). This is the standard MODFLOW approach.
        #       In 'matrix_only' mode the CRCH fractions will be set to 0, so
        #       no recharge reaches the conduit despite being in the RCH array.
        #
        #   'conduit_only':
        #       A 2-D array (n_rows × n_cols) is used (nrchop=2), with non-zero
        #       values only at the grid cells that contain a conduit inlet node.
        #       All other cells receive zero recharge. The CRCH fraction at inlet
        #       nodes is then set to 1.0 so the full cell recharge enters the conduit.

        # Steady-state uses the given intensity.
        if ss_rch is not None:
            rch_mean_for_ss = ss_rch if ss_rch > 0 else 0.0
            rch_mean_m_per_s = rch_mean_for_ss / 86400 / 1000
        else:
            rch_mean_for_ss = np.mean(recharge_ts)
            rch_mean_for_ss = rch_mean_for_ss if rch_mean_for_ss > 0. else 0.0
            rch_mean_m_per_s = rch_mean_for_ss / 86400 / 1000

        if routing_mode in ("both", "matrix_only"):
            # Uniform scalar rate applied to the whole grid
            rech_dict = {}
            rech_dict[0] = rch_mean_m_per_s  # stress period 0: steady-state warm-up

            for i, r in enumerate(recharge_ts):
                rch_m_per_s = r / 86400 / 1000 # convert from mm/d to m/s
                rech_dict[i + 1] = rch_m_per_s  # stress periods 1…n: daily transient

            # nrchop=1: one uniform rate per stress period over the whole top layer
            flopy.modflow.mfrch.ModflowRch(mf, nrchop=1, ipakcb=50, rech=rech_dict)

        else:
            # routing_mode == 'conduit_only':
            # Build a 2-D recharge array (n_rows × n_cols) where all cells are zero
            # except the cells that contain a conduit inlet node.
            # Each element of network.idxs_inlet is [col_1based, row_1based, lay_1based].
            # Convert to 0-based numpy indices: row_0 = inlet[1]-1, col_0 = inlet[0]-1.

            def _inlet_rch_array(rate_m_per_s):
                arr = np.zeros((n_rows, n_cols))
                for inlet in network.idxs_inlet:
                    arr[inlet[1] - 1, inlet[0] - 1] = rate_m_per_s
                return arr

            rech_dict = {}
            rech_dict[0] = _inlet_rch_array(rch_mean_m_per_s)  # steady-state

            for i, r in enumerate(recharge_ts):
                rech_dict[i + 1] = _inlet_rch_array(r / 86400 / 1000)  # transient

            # nrchop=2: MODFLOW reads one recharge value per top-layer cell
            flopy.modflow.mfrch.ModflowRch(mf, nrchop=2, ipakcb=50, rech=rech_dict)

        # --- Write all MODFLOW input files to disk ---
        mf.write_input()

        # --- Export the conduit network and generate the NBR file ---
        # export_network() writes the conduit map in the format CFPy requires
        network.validator.export_network()
        # generate_nbr() creates the node-branch-reach (NBR) file, which
        # describes the topology of the conduit network for MODFLOW-CFP
        network.validator.generate_nbr(
            path=modflow_dir,
            nrows=n_rows, ncols=n_cols, nlays=n_lays, nplanes=1,
            layer_elevations=lay_elevs_array,
        )

        # --- Read the NBR file to get node and pipe connectivity data ---
        nbr = cfpy.nbr()
        bot_elev, cond_elev = nbr.nbr_read()
        nbr_data = nbr.nbr(bot_elev, cond_elev)
        # nbr_data contents (by index):
        #   [0] node numbers (list of ints)
        #   [1] MODFLOW cell flat indices
        #   [2] MODFLOW [col, row, lay] indices for each node
        #   [3] node bottom elevations
        #   [4] ...
        #   [5] pipe (tube) numbers

        n_nodes = len(nbr_data[0])
        n_pipes = len(nbr_data[5])

        # --- Build CFP pipe data ---
        # All pipes are assigned the same uniform diameter, tortuosity,
        # roughness height, and Reynolds number thresholds
        pipe_data = [
            nbr_data[5],                                    # pipe numbers
            (np.ones(n_pipes) * diameter).tolist(),         # diameter (m)
            (np.ones(n_pipes) * self.tortuosity).tolist(),  # tortuosity (-)
            (np.ones(n_pipes) * rheight).tolist(),          # roughness height (m)
            (np.ones(n_pipes) * self.lcritrey).tolist(),    # lower critical Reynolds number
            (np.ones(n_pipes) * self.hcritrey).tolist(),    # upper critical Reynolds number
        ]

        # --- Assign nodal head boundary conditions ---
        # Default is -1 (free head); only the spring node gets a fixed head
        n_head_list = (np.ones(n_nodes) * -1).tolist()

        # Locate the spring in the node list and assign the constant head
        if network.idxs_spring not in nbr_data[2]:
            raise ValueError(f"Spring node {network.idxs_spring} not found in network.")
        spring_node_position = nbr_data[2].index(network.idxs_spring)
        n_head_list[spring_node_position] = chb_spring

        # --- Build conduit exchange and storage coefficient arrays ---
        # k_exchange: uniform exchange coefficient between conduit and matrix
        kex_data  = [nbr_data[0], (np.ones(n_nodes) * k_exchange).tolist()]
        # cads: uniform conduit storage coefficient
        cads_data = (np.ones(n_nodes) * cad).tolist()

        # --- Write the CFP input file ---
        cfp_str = cfpy.cfp(
            mode=1, nnodes=n_nodes, npipes=n_pipes, nlay=n_lays,
            nbr_data=nbr_data, geoheight=cond_elev, sa_exchange=0,
            epsilon=self.cfptol, niter=2000, relax=self.cfprelax,
            p_nr=0, cond_data=pipe_data,
            n_head=[nbr_data[0], n_head_list],
            k_exchange=kex_data, ncl=0, cl=0, ltemp=10,
            condl_data=0, cads=cads_data,
        ).cfp()

        # --- Write the COC input file (conduit output control) ---
        coc_str = cfpy.coc(
            nnodes=n_nodes, node_numbers=nbr_data[0], n_nts=1,
            npipes=n_pipes, pipe_numbers=nbr_data[5], t_nts=1,
        ).coc()

        # --- Build the CRCH recharge fraction array ---
        # p_crch[i] is the fraction of the MODFLOW RCH value in the cell that hosts
        # conduit node i which is routed directly into that conduit node.
        #
        # Routing mode logic:
        #   'both':         inlet nodes get crch_fraction; all other nodes get 0.
        #                   crch_fraction=1.0 (the default) reproduces the original
        #                   hard-coded behaviour.
        #   'matrix_only':  all nodes get 0. No recharge reaches the conduit.
        #   'conduit_only': inlet nodes get 1.0; all other nodes get 0.
        #                   crch_fraction is ignored here because the RCH array
        #                   already concentrates all recharge at the inlet cells,
        #                   so p_crch=1.0 routes 100 % of that cell recharge into
        #                   the conduit.

        p_crch = np.zeros(n_nodes).tolist()  # start with all zeros

        if routing_mode == "matrix_only":
            # All p_crch remain 0: the conduit receives no recharge at all.
            pass

        elif routing_mode == "both":
            # Each inlet node gets crch_fraction of its cell's RCH value.
            for idx in network.idxs_inlet:
                if idx not in nbr_data[2]:
                    raise ValueError(
                        f"Inlet node {idx} not found in the conduit network. "
                        "This inlet was generated by pyKasso but has no matching "
                        "conduit node in the NBR file."
                    )
                inlet_node_position = nbr_data[2].index(idx)
                p_crch[inlet_node_position] = float(crch_fraction)

        else:
            # routing_mode == 'conduit_only':
            # p_crch=1.0 at inlet nodes routes all of the spatially concentrated
            # RCH into the conduit.
            for idx in network.idxs_inlet:
                if idx not in nbr_data[2]:
                    raise ValueError(
                        f"Inlet node {idx} not found in the conduit network."
                    )
                inlet_node_position = nbr_data[2].index(idx)
                p_crch[inlet_node_position] = 1.0

        # --- Write the CRCH input file ---
        crch_str = cfpy.crch(
            iflag_crch=1, nper=n_pers,
            node_numbers=nbr_data[0], p_crch=p_crch,
        ).crch()

        # --- Write all CFP-related input files to disk ---
        cfpy.write_input(
            modelname=modelname,
            data_strings=[coc_str, crch_str, cfp_str],
            file_extensions=["coc", "crch", "cfp"],
        ).write_input()

        # --- Patch the MODFLOW .nam file to include the CFP package unit numbers ---
        cfpy.update_nam(
            modelname=modelname, mode=1,
            cfp_unit_num=52, crch_unit_num=53, coc_unit_num=54,
        ).update_nam()

        # --- Run the coupled MODFLOW-CFP model ---
        success, _ = mf.run_model(silent=True)

        # --- Read hydraulic head output ---
        hds   = bf.HeadFile(os.path.join(modflow_dir, f"{modelname}.hds"))
        heads = hds.get_alldata()  # shape: (n_pers, n_lays, n_rows, n_cols)

        # --- Read spring discharge from the CFP node output file ---
        # Find which 1-based node number corresponds to the spring cell index
        node_arr = np.array(nbr_data[2])
        spring_match       = np.where((node_arr == network.idxs_spring).all(axis=1))[0]
        spring_node_number = spring_match[0] + 1  # convert 0-based position to 1-based node number

        # CFP writes one output file per node: NODE########.OUT (8-digit zero-padded)
        node_filename = f"NODE{str(spring_node_number).zfill(8)}.OUT"
        flows = pd.read_fwf(os.path.join(modflow_dir, node_filename))

        # Spring discharge is the sum of all flow components at the spring node:
        #   QMAT:  exchange flux from the surrounding rock matrix
        #   QCADS: release from conduit storage
        #   QTUB1–QTUB6: inflow from up to 6 connected conduit pipe segments
        Q_out = (
            flows["QMAT"]
            + flows["QCADS"]
            + flows["QTUB1"]
            + flows["QTUB2"]
            + flows["QTUB3"]
            + flows["QTUB4"]
            + flows["QTUB5"]
            + flows["QTUB6"]
        )

        return success, Q_out, heads


# ── System ────────────────────────────────────────────────────────────────────

class System:
    """
    Central object for a virtual karst laboratory simulation.

    System manages all components of the simulation pipeline:
      - Climate data (Data)
      - Epikarst recharge model (Epikarst)
      - Groundwater flow model (Groundwater)
      - Four parameter groups (RechargeParams, NetworkStructure,
        NetworkProperties, AquiferParams)

    Each parameter group can be set to stochastic=True (values drawn from a
    Latin Hypercube Sample table) or stochastic=False (fixed default values).
    The LHS table is generated once and stored in self._param_samples. An
    integer random_state selects one row from that table; random_state=None
    uses the default (fixed) parameter values for all groups.

    Example usage
    -------------
    # Run with all default parameters (no randomness):
    s = vlab.System()
    rch = s.simulate_recharge(random_state=None)
    success, Q_out, heads = s.simulate_full(random_state=None)

    # Run the 42nd LHS realisation:
    success, Q_out, heads = s.simulate_full(random_state=42)

    # Mix stochastic and fixed groups:
    s = vlab.System(
        recharge_params=vlab.RechargeParams(stochastic=True),
        aquifer_params=vlab.AquiferParams(stochastic=False, values={"hk": 5e-5}),
    )

    Parameters
    ----------
    name : str
        Label for this simulation (used in output directory names).
    n_rows : int
        Number of model grid rows (north–south direction).
    n_cols : int
        Number of model grid columns (east–west direction).
    n_lays : int
        Number of model layers (vertical).
    delr : float
        Cell size in the column (east–west) direction (m).
    delc : float
        Cell size in the row (north–south) direction (m).
    lay_elevs : list of float, optional
        [top_elevation, bottom_elevation] of the model layer (m).
        Defaults to [50.0, 0.0].
    chb_spring : float
        Fixed hydraulic head at the spring outlet node (m).
    data : Data, optional
        Climate input. Defaults to Data.default() (bundled CSVs).
    epikarst : Epikarst, optional
        Epikarst recharge model. Defaults to Epikarst() with n_tiles=5.
    groundwater : Groundwater, optional
        Groundwater flow model. Defaults to Groundwater() with standard settings.
    n_stress_periods : int or None
        Number of transient stress periods to simulate. None uses the full
        length of the recharge time series. Must be set (not None) when
        recharge_source='block'.
    recharge_params : RechargeParams, optional
        Defaults to RechargeParams() with stochastic=True.
    network_structure : NetworkStructure, optional
        Defaults to NetworkStructure() with stochastic=True.
    network_properties : NetworkProperties, optional
        Defaults to NetworkProperties() with stochastic=True.
    aquifer_params : AquiferParams, optional
        Defaults to AquiferParams() with stochastic=True.
    temporal_resolution : str
        'daily' (default) or 'hourly'. Sets the time step for the full
        simulation pipeline. With 'daily', each MODFLOW stress period is
        86400 s (one day) and climate data are expected in mm/d. With
        'hourly', each stress period is 3600 s (one hour) and climate data
        are expected in mm/h. The MODFLOW recharge package always receives
        rates in m/s; the conversion factor (/ 86400 / 1000) is the same
        for both resolutions because recharge is internally always kept as
        an instantaneous mm/d rate.
    recharge_source : str
        'epikarst' (default) runs the FlexModelFrac epikarst model.
        'block' uses a simple rectangular pulse of fixed intensity followed
        by zeros, with total length equal to n_stress_periods.
    ss_rch : float
        The steady-state recharge intensity (mm/d). If it is None, the mean
        of the recharge time series (or block) is used instead.
    block_intensity : float
        Recharge rate during the block pulse (mm/d). Only used when
        recharge_source='block'.
    block_length : int
        Number of time steps the block pulse lasts. Must be <= n_stress_periods.
        Only used when recharge_source='block'.
    routing_mode : str
        Controls where recharge is applied in the model:
        'both' (default): uniform RCH to all matrix cells, plus a fraction
            (crch_fraction from NetworkProperties) of each inlet cell's
            recharge is routed directly into the conduit.
        'matrix_only': uniform RCH to all matrix cells, no conduit recharge
            (CRCH fractions all set to 0).
        'conduit_only': RCH is zero everywhere except the inlet (sinkhole)
            grid cells; all of that recharge goes to the conduit (p_crch=1).
    """

    def __init__(
        self,
        name="vlab",
        n_rows=50,
        n_cols=100,
        n_lays=1,
        delr=25.0,
        delc=25.0,
        lay_elevs=None,
        chb_spring=10.0,
        data=None,
        epikarst=None,
        groundwater=None,
        n_stress_periods=None,
        recharge_params=None,
        ss_rch=None,
        network_structure=None,
        network_properties=None,
        aquifer_params=None,
        temporal_resolution="daily",
        recharge_source="epikarst",
        block_intensity=10.0,
        block_length=30,
        routing_mode="both",
    ):
        # --- Grid geometry ---
        self.name             = name
        self.n_rows           = n_rows
        self.n_cols           = n_cols
        self.n_lays           = n_lays
        self.delr             = delr
        self.delc             = delc
        self.chb_spring       = chb_spring
        self.n_stress_periods = n_stress_periods

        # Use provided layer elevations, or fall back to the default (50 m thick aquifer)
        if lay_elevs is not None:
            self.lay_elevs = lay_elevs
        else:
            self.lay_elevs = [50.0, 0.0]

        # --- Temporal resolution ---
        _valid_resolutions = {"daily", "hourly"}
        if temporal_resolution not in _valid_resolutions:
            raise ValueError(
                f"temporal_resolution must be one of {_valid_resolutions}, "
                f"got '{temporal_resolution}'."
            )
        self.temporal_resolution = temporal_resolution

        # --- Sub-model components ---
        if data is not None:
            self.data = data
        else:
            # Load the default climate data at the correct temporal resolution
            self.data = Data.default(temporal_resolution=self.temporal_resolution)

        if epikarst is not None:
            self.epikarst = epikarst
        else:
            self.epikarst = Epikarst()

        if groundwater is not None:
            self.groundwater = groundwater
        else:
            self.groundwater = Groundwater()

        # Hourly mode requires 3600 s stress periods instead of the default 86400 s.
        # Override perlen here so the user does not have to set it manually.
        if self.temporal_resolution == "hourly":
            self.groundwater.perlen = 3600

        # --- Parameter groups ---
        if recharge_params is not None:
            self.recharge_params = recharge_params
        else:
            self.recharge_params = RechargeParams()

        if network_structure is not None:
            self.network_structure = network_structure
        else:
            self.network_structure = NetworkStructure()

        if network_properties is not None:
            self.network_properties = network_properties
        else:
            self.network_properties = NetworkProperties()

        if aquifer_params is not None:
            self.aquifer_params = aquifer_params
        else:
            self.aquifer_params = AquiferParams()

        # --- Recharge source ---
        _valid_sources = {"epikarst", "block"}
        if recharge_source not in _valid_sources:
            raise ValueError(
                f"recharge_source must be one of {_valid_sources}, "
                f"got '{recharge_source}'."
            )
        self.recharge_source = recharge_source
        self.block_intensity = float(block_intensity)

        if ss_rch is not None:
            self.ss_rch = float(ss_rch)
        else:
            self.ss_rch = None

        self.block_length    = int(block_length)

        # --- Recharge routing mode ---
        _valid_modes = {"both", "matrix_only", "conduit_only"}
        if routing_mode not in _valid_modes:
            raise ValueError(
                f"routing_mode must be one of {_valid_modes}, "
                f"got '{routing_mode}'."
            )
        self.routing_mode = routing_mode

        # LHS parameter table — generated lazily the first time param_samples is accessed
        self._param_samples: pd.DataFrame | None = None

    def _groups(self):
        """
        Return all four parameter groups as an ordered list.

        Used internally to iterate over groups during LHS sampling and parameter
        assembly. The order here determines the column order in param_samples.
        """
        return [
            self.recharge_params,
            self.network_structure,
            self.network_properties,
            self.aquifer_params,
        ]

    # ── Parameter sampling ────────────────────────────────────────────────────

    def generate_samples(self, n_samples=500, seed=1) -> pd.DataFrame:
        """
        Generate a joint Latin Hypercube Sampling (LHS) parameter table.

        Draws n_samples parameter sets across all stochastic dimensions at once
        (joint sampling), so the space-filling property of LHS holds across all
        parameters together. Fixed groups (stochastic=False) are not included
        in the table; their values are handled separately in _get_params().

        After sampling, the full table is split back into per-group arrays
        and stored in each group's ._samples attribute so that _get_params()
        can quickly retrieve the value for a given random_state.

        Parameters
        ----------
        n_samples : int
            Number of parameter sets (rows) to generate.
        seed : int
            Random seed for the LHS sampler (ensures reproducibility).

        Returns
        -------
        param_samples : pd.DataFrame, shape (n_samples, n_stochastic_dims)
            One column per stochastic parameter, one row per realisation.
        """
        # Collect only the groups that are marked as stochastic
        stochastic_groups = []
        for g in self._groups():
            if g.stochastic:
                stochastic_groups.append(g)

        # Build ordered lists of parameter names and their lower/upper bounds
        all_names = []
        all_lows  = []
        all_highs = []
        for g in stochastic_groups:
            for name in g.NAMES:
                all_names.append(name)
                low, high = g.bounds[name]
                all_lows.append(low)
                all_highs.append(high)

        n_dims = len(all_names)

        # If every group is fixed, return an empty DataFrame and exit early
        if n_dims == 0:
            self._param_samples = pd.DataFrame(index=range(n_samples))
            for g in self._groups():
                g._samples = np.empty((n_samples, 0))
            return self._param_samples

        # --- Draw joint LHS samples in the unit hypercube [0, 1]^n_dims ---
        lhs_sampler = qmc.LatinHypercube(d=n_dims, rng=seed)
        raw_samples = lhs_sampler.random(n=n_samples)

        # Scale each dimension from [0, 1] to its physical [low, high] range
        scaled_samples = qmc.scale(raw_samples, all_lows, all_highs)

        # --- Distribute the scaled columns back to each stochastic group ---
        # Each group receives a slice of columns equal to its number of parameters
        col_start = 0
        for g in stochastic_groups:
            n_params  = len(g.NAMES)
            g._samples = scaled_samples[:, col_start:col_start + n_params].copy()
            col_start += n_params

        # sks_seed must be a whole number (it is used as a random seed in pyKasso)
        if self.network_structure.stochastic:
            seed_col = self.network_structure.NAMES.index("sks_seed")
            self.network_structure._samples[:, seed_col] = np.round(
                self.network_structure._samples[:, seed_col]
            ).astype(int)

        # --- Assemble and store the combined DataFrame ---
        self._param_samples = pd.DataFrame(scaled_samples, columns=all_names)
        if "sks_seed" in self._param_samples.columns:
            self._param_samples["sks_seed"] = (
                self._param_samples["sks_seed"].round().astype(int)
            )

        return self._param_samples

    @property
    def param_samples(self) -> pd.DataFrame:
        """
        The LHS parameter table.

        Generated automatically with default settings (n_samples=500, seed=1)
        on first access. Call generate_samples() explicitly beforehand to
        control the number of samples or the random seed.
        """
        if self._param_samples is None:
            self.generate_samples()
        return self._param_samples

    def _get_params(self, random_state) -> dict:
        """
        Assemble a full parameter dictionary for a single simulation run.

        For stochastic groups: uses the LHS sample at row random_state if an
        integer is given, or the group's default values if random_state=None.
        For fixed groups (stochastic=False): always uses the default values.

        Fixed structural parameters from NetworkStructure (n_inlets, seeds,
        fracture families) are always included, regardless of stochastic setting.

        Parameters
        ----------
        random_state : int or None
            Row index into param_samples. None uses default values for all groups.

        Returns
        -------
        params : dict
            Full parameter dictionary with one key per parameter name.
        """
        # Trigger LHS table generation if needed and if we have a stochastic run
        if random_state is not None:
            _ = self.param_samples

        params = {}
        for g in self._groups():
            if g.stochastic and random_state is not None:
                # Use the sampled value from the LHS table for this realisation
                for i, name in enumerate(g.NAMES):
                    val = g._samples[random_state, i]
                    # sks_seed must be an integer; all other parameters are floats
                    if name == "sks_seed":
                        params[name] = int(round(val))
                    else:
                        params[name] = float(val)
            else:
                # Use the fixed default values stored in the group
                for name in g.NAMES:
                    params[name] = g.values[name]

        # Always append the fixed structural parameters from NetworkStructure.
        # These define the aquifer geometry and are never part of the LHS table.
        params["n_inlets"]          = self.network_structure.n_inlets
        params["inlets_seed"]       = self.network_structure.inlets_seed
        params["outlets_seed"]      = self.network_structure.outlets_seed
        params["fracture_families"] = self.network_structure.fracture_families

        return params

    # ── Simulate recharge ─────────────────────────────────────────────────────

    def simulate_recharge(self, random_state=None) -> np.ndarray:
        """
        Simulate the recharge time series for one parameter realisation.

        When recharge_source='epikarst', runs the FlexModelFrac epikarst model
        on the loaded climate data and returns total (slow + quick) recharge.

        When recharge_source='block', returns a rectangular pulse: block_intensity
        (mm/d) for block_length steps, then zeros for the remaining
        n_stress_periods - block_length steps. n_stress_periods must be set and
        must be >= block_length. The random_state parameter is ignored in block
        mode because there are no stochastic recharge parameters.

        Parameters
        ----------
        random_state : int or None
            Row index into param_samples. Used only when recharge_source='epikarst'.

        Returns
        -------
        rch_ts : np.ndarray, shape (n_steps,)
            Total recharge (mm/d). Length equals n_stress_periods in block mode
            and len(data.prec) in epikarst mode (before any truncation by
            n_stress_periods in the groundwater model).
        """
        if self.recharge_source == "block":
            # Validate that the full simulation length has been specified
            if self.n_stress_periods is None:
                raise ValueError(
                    "n_stress_periods must be set when recharge_source='block'. "
                    "Set it in System(..., n_stress_periods=N)."
                )
            if self.block_length > self.n_stress_periods:
                raise ValueError(
                    f"block_length ({self.block_length}) exceeds "
                    f"n_stress_periods ({self.n_stress_periods}). "
                    "The pulse cannot be longer than the simulation period."
                )
            # Rectangular pulse followed by zeros to fill the full simulation period
            block = np.full(self.block_length, self.block_intensity)
            tail  = np.zeros(self.n_stress_periods - self.block_length)
            return np.concatenate([block, tail])

        else:
            # recharge_source == "epikarst" (the default)
            # Get the full parameter dict for this realisation
            params = self._get_params(random_state)

            # Extract only the recharge model parameters in the order FlexModelFrac expects
            rch_params = [params[name] for name in RechargeParams.NAMES]

            return self.epikarst.simulate(
                self.data,
                rch_params,
                temporal_resolution=self.temporal_resolution,
            )

    # ── Simulate network only ─────────────────────────────────────────────────

    def generate_network(self, random_state=None, working_dir=None, cleanup=True,
                         plot_network=False) -> "Network":
        """
        Run only the pyKasso network generation phase.

        Use this to generate a conduit network separately from the MODFLOW run.
        The returned Network object can then be passed to simulate_from_network()
        to run multiple MODFLOW simulations on the same network geometry.

        Parameters
        ----------
        random_state : int or None
            Row index into param_samples. None uses the default network parameters.
        working_dir : str or Path, optional
            Base directory for model files. Defaults to the current working directory.
        cleanup : bool
            If True, delete the pyKasso run directory after generation.
        plot_network : bool
            If True, show a map of the generated conduit network.

        Returns
        -------
        network : Network
            Validated network with spring and inlet MODFLOW cell indices.
        """
        params = self._get_params(random_state)

        # Build a unique run tag from the simulation name and random_state
        tag       = str(random_state) if random_state is not None else "default"
        proj_name = f"{self.name}_{tag}"

        # Use the provided working directory, or fall back to the current directory
        if working_dir is not None:
            base_dir = str(working_dir)
        else:
            base_dir = str(Path.cwd())

        # Store the current directory so we can restore it after the run
        orig_dir = os.getcwd()

        try:
            os.chdir(base_dir)
            _flush_pykasso_loggers()
            network = self.groundwater.generate_network(
                system=self,
                params=params,
                proj_name=proj_name,
                working_dir=base_dir,
                plot_network=plot_network,
            )
        finally:
            # Always restore the working directory and release logger file locks,
            # even if an exception was raised inside the try block
            _flush_pykasso_loggers()
            os.chdir(orig_dir)
            if cleanup:
                pykasso_dir = os.path.join(base_dir, proj_name + "_pyKasso")
                if os.path.exists(pykasso_dir):
                    shutil.rmtree(pykasso_dir)

        return network

    # ── Simulate from pre-generated network ───────────────────────────────────

    def simulate_from_network(self, network: "Network", random_state=None,
                              working_dir=None, cleanup=True):
        """
        Run MODFLOW + CFP on a pre-generated conduit network.

        Recharge is simulated internally using the same random_state that was
        used during generate_network(), keeping the recharge and aquifer
        parameters consistent with the network that was generated.

        Parameters
        ----------
        network : Network
            Output of System.generate_network().
        random_state : int or None
            Row index into param_samples. None uses the default parameter values.
        working_dir : str or Path, optional
            Base directory for model files. Defaults to the current working directory.
        cleanup : bool
            If True, delete the MODFLOW run directory after the simulation.

        Returns
        -------
        success : bool
            True if MODFLOW converged.
        Q_out : pd.Series
            Spring discharge time series (m³/s).
        heads : np.ndarray, shape (n_pers, n_lays, n_rows, n_cols)
            Simulated hydraulic head for each stress period and layer.
        """
        params = self._get_params(random_state)

        # Simulate the recharge time series for this realisation.
        # Using simulate_recharge() here (rather than calling epikarst directly)
        # ensures that block mode is respected when simulate_from_network() is
        # called on its own, without going through simulate_full().
        rch_ts = self.simulate_recharge(random_state)

        # Build run tag and directory paths
        tag       = str(random_state) if random_state is not None else "default"
        proj_name = f"{self.name}_{tag}"

        if working_dir is not None:
            base_dir = str(working_dir)
        else:
            base_dir = str(Path.cwd())

        # Store the current directory so we can restore it after the run
        orig_dir = os.getcwd()

        try:
            os.chdir(base_dir)
            _flush_pykasso_loggers()
            success, Q_out, heads = self.groundwater.run_modflow(
                system=self,
                network=network,
                recharge_ts=rch_ts,
                ss_rch=self.ss_rch,
                params=params,
                proj_name=proj_name,
                working_dir=base_dir,
                n_stress_periods=self.n_stress_periods,
                routing_mode=self.routing_mode,
            )
        finally:
            # Always restore the working directory and release logger file locks
            _flush_pykasso_loggers()
            os.chdir(orig_dir)
            if cleanup:
                modflow_dir = os.path.join(base_dir, proj_name + "_modflow")
                if os.path.exists(modflow_dir):
                    shutil.rmtree(modflow_dir)

        return success, Q_out, heads

    # ── Simulate full ─────────────────────────────────────────────────────────

    def simulate_full(self, random_state=None, working_dir=None, cleanup=True,
                      plot_network=False):
        """
        Run the complete simulation pipeline for one parameter realisation.

        Steps performed:
          1. Simulate epikarst recharge (FlexModelFrac)
          2. Generate a stochastic karst conduit network (pyKasso)
          3. Run the coupled groundwater flow model (MODFLOW + CFP)

        Parameters
        ----------
        random_state : int or None
            Row index into param_samples. None uses the default parameter values.
        working_dir : str or Path, optional
            Base directory for model files. Defaults to the current working directory.
        cleanup : bool
            If True, remove pyKasso and MODFLOW run directories after the simulation.
        plot_network : bool
            If True, show a map of the generated conduit network.

        Returns
        -------
        success : bool
            True if MODFLOW converged.
        Q_out : pd.Series
            Spring discharge time series (m³/s).
        heads : np.ndarray, shape (n_pers, n_lays, n_rows, n_cols)
            Simulated hydraulic head for each stress period and layer.
        """
        # Steps 1 + 2: generate the conduit network (recharge is simulated inside)
        network = self.generate_network(
            random_state=random_state,
            working_dir=working_dir,
            cleanup=cleanup,
            plot_network=plot_network,
        )

        # Step 3: run MODFLOW + CFP on the generated network
        return self.simulate_from_network(
            network=network,
            random_state=random_state,
            working_dir=working_dir,
            cleanup=cleanup,
        )


# ── Helpers ───────────────────────────────────────────────────────────────────

def _flush_pykasso_loggers():
    """
    Close and remove all handlers from the root logger.

    pyKasso registers file handlers when it creates a project. On Windows,
    these handlers hold open file locks that prevent the temporary pyKasso
    directory from being deleted (e.g. during cleanup with shutil.rmtree()).
    Removing the handlers here releases those locks before the directory
    removal is attempted.
    """
    root = logging.getLogger()
    for h in root.handlers[:]:
        h.close()
        root.removeHandler(h)
