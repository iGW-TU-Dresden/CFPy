import os
import shutil
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pykasso as pk
import CFPy as cfpy
import flopy
import flopy.utils.binaryfile as bf
from scipy.stats import qmc

from recharge import FlexModelFrac

# ── Module-level constants ────────────────────────────────────────────────────

PARAM_NAMES = [
    "srmax", "lp", "ks", "gamma", "simax", "kv", "omega",
    "hk", "k_exchange", "cad", "diameter",
]
PARAM_LOW  = [100., 0.1, 100., 1., 0.5, 0.1, 0.1, 1e-5, 1e-5, 0.01, 0.1]
PARAM_HIGH = [500., 0.8, 1000., 5., 10., 10., 0.9, 1e-4, 1e-4, 0.5, 1.0]
VALID_SEEDS = [12345, 1234, 234, 345]

# Column order: [sks_seed, srmax, lp, ks, gamma, simax, kv, omega,
#                hk, k_exchange, cad, diameter]
DEFAULT_PARAMS = np.array(
    [12345, 250., 0.25, 1000., 2., 1., 5., 0.2, 1e-4, 1e-4, 0.01, 0.5]
)

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

_VLAB_DIR = Path(__file__).parent


# ── Data ──────────────────────────────────────────────────────────────────────

class Data:
    """Daily precipitation and evaporation time series."""

    def __init__(self, prec: pd.Series, evap: pd.Series):
        self.prec = prec
        self.evap = evap

    @classmethod
    def from_csv(cls, prec_path, evap_path):
        prec = pd.read_csv(
            prec_path, parse_dates=True, index_col=[0],
            sep=";", usecols=[0, 2], dayfirst=True,
        ).resample("D").sum().squeeze()
        evap = pd.read_csv(
            evap_path, parse_dates=True, index_col=[0],
            sep=";", dayfirst=True,
        ).squeeze()
        return cls(prec, evap)

    @classmethod
    def default(cls):
        return cls.from_csv(_VLAB_DIR / "prec.csv", _VLAB_DIR / "evap.csv")


# ── Epikarst ──────────────────────────────────────────────────────────────────

class Epikarst:
    """Soil/epikarst recharge model (FlexModelFrac by default)."""

    def __init__(self, n_tiles: int = 5):
        self.n_tiles = n_tiles
        self._model = FlexModelFrac()

    def simulate(self, data: Data, rch_params) -> np.ndarray:
        """
        Simulate total recharge (slow + quick) for the study period.

        Parameters
        ----------
        data : Data
        rch_params : array-like, shape (7,)
            [srmax, lp, ks, gamma, simax, kv, omega]

        Returns
        -------
        rch_ts : np.ndarray, shape (len(data.prec),), mm/d
        """
        n = len(data.prec)
        prec_ml = np.tile(data.prec.values.ravel(), self.n_tiles)
        evap_ml = np.tile(data.evap.values.ravel(), self.n_tiles)

        rch_slow, rch_quick = self._model.simulate(
            prec=prec_ml, evap=evap_ml, p=np.array(rch_params)
        )
        # discard warmup repetitions
        return (rch_slow + rch_quick)[n * (self.n_tiles - 1):]


# ── Groundwater ───────────────────────────────────────────────────────────────

class Groundwater:
    """pyKasso network generation + MODFLOW + CFP simulation."""

    def __init__(
        self,
        exe_name="C:/WRDAPP/CFPv2.exe",
        nz=10,
        valid_seeds=None,
        fracture_families=None,
        ss=1e-4,
        sy=1e-2,
        vka_ratio=1.0,
        H_init=75.0,
        mftol=1e-4,
        mfrelax=0.98,
        mxiter=5000,
        tortuosity=1.0,
        lcritrey=500,
        hcritrey=5000,
        cfptol=1e-5,
        cfprelax=0.98,
        rheight=0.02,
        perlen=86400,
        time_unit=1,
    ):
        self.exe_name = exe_name
        self.nz = nz
        self.valid_seeds = valid_seeds if valid_seeds is not None else VALID_SEEDS
        self.fracture_families = (
            fracture_families if fracture_families is not None else _DEFAULT_FRACTURES
        )
        self.ss = ss
        self.sy = sy
        self.vka_ratio = vka_ratio
        self.H_init = H_init
        self.mftol = mftol
        self.mfrelax = mfrelax
        self.mxiter = mxiter
        self.tortuosity = tortuosity
        self.lcritrey = lcritrey
        self.hcritrey = hcritrey
        self.cfptol = cfptol
        self.cfprelax = cfprelax
        self.rheight = rheight
        self.perlen = perlen
        self.time_unit = time_unit

    def simulate(self, system, recharge_ts, params, proj_name, working_dir,
                 n_stress_periods=None):
        """
        Run pyKasso + MODFLOW + CFP.

        Parameters
        ----------
        system : System
        recharge_ts : np.ndarray  full recharge time series (mm/d)
        params : array-like, shape (5,)  [hk, k_exchange, cad, diameter, sks_seed]
        proj_name : str
        working_dir : str or Path
        n_stress_periods : int or None  truncate recharge to this length

        Returns
        -------
        success : bool
        Q_out   : pd.Series  spring discharge time series
        heads   : np.ndarray  shape (n_pers, n_lays, n_rows, n_cols)
        """
        hk         = float(params[0])
        k_exchange = float(params[1])
        cad        = float(params[2])
        diameter   = float(params[3])
        sks_seed   = int(params[4])

        if n_stress_periods is not None:
            recharge_ts = recharge_ts[:n_stress_periods]

        return self._run(
            proj_name=proj_name,
            recharge_ts=recharge_ts,
            sks_seed=sks_seed,
            hk=hk,
            k_exchange=k_exchange,
            cad=cad,
            diameter=diameter,
            n_rows=system.n_rows,
            n_cols=system.n_cols,
            n_lays=system.n_lays,
            delr=system.delr,
            delc=system.delc,
            lay_elevs=system.lay_elevs,
            chb_spring=system.chb_spring,
            working_dir=str(working_dir),
        )

    def _run(self, proj_name, recharge_ts, sks_seed, hk, k_exchange, cad,
             diameter, n_rows, n_cols, n_lays, delr, delc, lay_elevs,
             chb_spring, working_dir):

        # ── pyKasso ───────────────────────────────────────────────────────────
        os.chdir(working_dir)

        dz = (lay_elevs[0] - lay_elevs[-1]) / self.nz
        grid_parameters = {
            "x0": 0, "y0": 0, "z0": 0,
            "nx": n_cols, "ny": n_rows, "nz": self.nz,
            "dx": delr, "dy": delc, "dz": dz,
        }

        app = pk.pykasso()
        app.new_project(name=proj_name + "_pyKasso", grid_parameters=grid_parameters)
        app.visualizer.notebook = True

        model_parameters = {
            "sks": {"seed": sks_seed, "algorithm": "Riemann3", "ratio": 0.1},
            "domain": {
                "bedrock": np.zeros((n_cols, n_rows)),
                "water_table": np.full((n_cols, n_rows), chb_spring),
            },
            "outlets": {"seed": 7, "number": 1, "subdomain": "domain_borders_bottom"},
            "inlets":  {"number": 5, "seed": 456, "subdomain": "domain_surface"},
            "fractures": {"generate": self.fracture_families},
        }
        app.model.generate(model_parameters=model_parameters)

        # ── Post-process pyKasso → CFPy inputs ────────────────────────────────
        network    = np.moveaxis(np.array(app.model.maps["karst"][0]), 1, 0)
        elevations = np.moveaxis(app.model.node_elev_arr, 1, 0)
        outlets    = np.moveaxis(app.model.maps["outlets"], 1, 0)

        network    = np.max(network, axis=-1)[::-1, :]
        elevations = np.max(elevations, axis=-1)[::-1, :]
        outlets    = np.nanmax(outlets, axis=-1)[::-1, :]

        validator = cfpy.preprocessing.GeneralValidator(
            network=network, elevations=elevations
        )
        validator.validate_network_elevations()

        # spring location (1-based MODFLOW indices: [col, row, lay])
        sp = np.argwhere(~np.isnan(outlets))[0]
        idxs_spring = [int(sp[1] + 1), int(sp[0] + 1), 1]

        # inlet locations
        df = pd.DataFrame.from_dict(
            app.model.vectors["nodes"], orient="index",
            columns=["x", "y", "z", "type"],
        )
        idxs_inlet = []
        for _, row in df[df["type"] == "inlet"].iterrows():
            col, r = app.model.grid.get_indices((row["x"], row["y"]))
            idxs_inlet.append([int(col[0]) + 1, n_rows - int(r[0]), 1])

        lay_elevs_array = [
            np.ones((n_rows, n_cols)) * lay_elevs[0],
            np.ones((n_rows, n_cols)) * lay_elevs[1],
        ]

        # ── MODFLOW setup ─────────────────────────────────────────────────────
        modflow_dir = os.path.join(working_dir, proj_name + "_modflow")
        os.makedirs(modflow_dir, exist_ok=True)
        os.chdir(modflow_dir)

        modelname = "pyKasso_example"
        mf = flopy.modflow.Modflow(modelname, exe_name=self.exe_name)

        n_pers = len(recharge_ts) + 1  # 1 steady + n transient
        flopy.modflow.ModflowDis(
            mf, n_lays, n_rows, n_cols, n_pers,
            delr, delc, top=lay_elevs[0], botm=lay_elevs[1],
            perlen=self.perlen,
            nstp=[1] * n_pers,
            steady=[1] + [0] * (n_pers - 1),
            itmuni=self.time_unit, lenuni=2,
        )

        flopy.modflow.ModflowBas(
            mf, ibound=np.ones((n_lays, n_rows, n_cols), dtype=np.int32),
            strt=self.H_init,
        )

        flopy.modflow.ModflowLpf(
            mf, laytyp=1, chani=1, layvka=0, laywet=0,
            ipakcb=50, hdry=999., hk=hk, hani=1,
            vka=hk * self.vka_ratio, wetdry=0, ss=self.ss, sy=self.sy,
        )

        flopy.modflow.ModflowOc(
            mf,
            stress_period_data={
                (i, 0): ["save head", "print budget"] for i in range(n_pers)
            },
        )

        flopy.modflow.ModflowPcg(
            mf, mxiter=self.mxiter, iter1=self.mxiter, npcond=1,
            hclose=self.mftol, rclose=self.mftol, relax=self.mfrelax,
            nbpol=2, iprpcg=5, mutpcg=0, damp=0.99, ihcofadd=9999,
        )

        rch_mean = np.mean(recharge_ts) / 86400 / 1000  # mm/d → m/s
        rech = {0: rch_mean}
        for i, r in enumerate(recharge_ts):
            rech[i + 1] = r / 86400 / 1000
        flopy.modflow.mfrch.ModflowRch(mf, nrchop=1, ipakcb=50, rech=rech)

        # ── CFP setup ─────────────────────────────────────────────────────────
        mf.write_input()
        validator.export_network()
        validator.generate_nbr(
            path=modflow_dir, nrows=n_rows, ncols=n_cols,
            nlays=n_lays, nplanes=1, layer_elevations=lay_elevs_array,
        )

        nbr = cfpy.nbr()
        bot_elev, cond_elev = nbr.nbr_read()
        nbr_data = nbr.nbr(bot_elev, cond_elev)

        n_nodes = len(nbr_data[0])
        n_pipes = len(nbr_data[5])

        pipe_data = [
            nbr_data[5],
            (np.ones(n_pipes) * diameter).tolist(),
            (np.ones(n_pipes) * self.tortuosity).tolist(),
            (np.ones(n_pipes) * self.rheight).tolist(),
            (np.ones(n_pipes) * self.lcritrey).tolist(),
            (np.ones(n_pipes) * self.hcritrey).tolist(),
        ]

        n_head = (np.ones(n_nodes) * -1).tolist()
        if idxs_spring not in nbr_data[2]:
            raise ValueError(f"Spring node {idxs_spring} not found in network.")
        n_head[nbr_data[2].index(idxs_spring)] = chb_spring

        kex_data  = [nbr_data[0], (np.ones(n_nodes) * k_exchange).tolist()]
        cads_data = (np.ones(n_nodes) * cad).tolist()

        cfp_str = cfpy.cfp(
            mode=1, nnodes=n_nodes, npipes=n_pipes, nlay=n_lays,
            nbr_data=nbr_data, geoheight=cond_elev, sa_exchange=0,
            epsilon=self.cfptol, niter=2000, relax=self.cfprelax,
            p_nr=0, cond_data=pipe_data,
            n_head=[nbr_data[0], n_head],
            k_exchange=kex_data, ncl=0, cl=0, ltemp=10,
            condl_data=0, cads=cads_data,
        ).cfp()

        coc_str = cfpy.coc(
            nnodes=n_nodes, node_numbers=nbr_data[0], n_nts=1,
            npipes=n_pipes, pipe_numbers=nbr_data[5], t_nts=1,
        ).coc()

        p_crch = np.zeros(n_nodes).tolist()
        for idx in idxs_inlet:
            if idx not in nbr_data[2]:
                raise ValueError(f"Inlet node {idx} not found in network.")
            p_crch[nbr_data[2].index(idx)] = 1.0

        crch_str = cfpy.crch(
            iflag_crch=1, nper=n_pers,
            node_numbers=nbr_data[0], p_crch=p_crch,
        ).crch()

        cfpy.write_input(
            modelname=modelname,
            data_strings=[coc_str, crch_str, cfp_str],
            file_extensions=["coc", "crch", "cfp"],
        ).write_input()

        cfpy.update_nam(
            modelname=modelname, mode=1,
            cfp_unit_num=52, crch_unit_num=53, coc_unit_num=54,
        ).update_nam()

        # ── Run ───────────────────────────────────────────────────────────────
        success, _ = mf.run_model(silent=True)

        # ── Read heads ────────────────────────────────────────────────────────
        hds = bf.HeadFile(os.path.join(modflow_dir, f"{modelname}.hds"))
        heads = hds.get_alldata()

        # ── Read spring discharge ─────────────────────────────────────────────
        node_arr = np.array(nbr_data[2])
        idx = np.where((node_arr == idxs_spring).all(axis=1))[0] + 1
        flows = pd.read_fwf(
            os.path.join(modflow_dir, f"NODE{str(idx[0]).zfill(8)}.OUT")
        )
        Q_out = (
            flows["QMAT"] + flows["QCADS"]
            + flows["QTUB1"] + flows["QTUB2"] + flows["QTUB3"]
            + flows["QTUB4"] + flows["QTUB5"] + flows["QTUB6"]
        )

        return success, Q_out, heads


# ── System ────────────────────────────────────────────────────────────────────

class System:
    """
    Central object for a virtual karst laboratory.

    Parameters
    ----------
    name : str
        Project name prefix used for run directories.
    n_rows, n_cols, n_lays : int
        Model grid dimensions.
    delr, delc : float
        Cell width along rows / columns (m).
    lay_elevs : list[float]
        [top, bottom] elevation of the single model layer (m).
    chb_spring : float
        Fixed head boundary at the spring node (m).
    data : Data or None
        Precipitation/evaporation. None loads the default CSVs.
    epikarst : Epikarst or None
        Recharge model. None uses FlexModelFrac with n_tiles=5.
    groundwater : Groundwater or None
        Network and flow model config. None uses all defaults.
    n_stress_periods : int or None
        Number of transient stress periods fed to MODFLOW. None uses
        the full recharge time series length.

    Examples
    --------
    >>> s = vlab.System()
    >>> rch = s.simulate_recharge(random_state=0)
    >>> success, Q, heads = s.simulate_full(random_state=0)
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
        n_stress_periods=5,
    ):
        self.name = name
        self.n_rows = n_rows
        self.n_cols = n_cols
        self.n_lays = n_lays
        self.delr = delr
        self.delc = delc
        self.lay_elevs = lay_elevs if lay_elevs is not None else [50.0, 0.0]
        self.chb_spring = chb_spring
        self.n_stress_periods = n_stress_periods

        self.data = data if data is not None else Data.default()
        self.epikarst = epikarst if epikarst is not None else Epikarst()
        self.groundwater = groundwater if groundwater is not None else Groundwater()

        self._param_samples = None

    # ── Parameter sampling ────────────────────────────────────────────────────

    def generate_samples(self, n_samples=500, seed=1) -> np.ndarray:
        """
        Generate LHS parameter samples and store in self._param_samples.

        Returns
        -------
        param_samples : np.ndarray, shape (n_samples, 12)
            Column order: [sks_seed, srmax, lp, ks, gamma, simax, kv, omega,
                           hk, k_exchange, cad, diameter]
        """
        raw = qmc.LatinHypercube(d=11, rng=seed).random(n=n_samples)
        scaled = qmc.scale(raw, PARAM_LOW, PARAM_HIGH)
        seeds = np.random.default_rng(seed).choice(VALID_SEEDS, n_samples)
        self._param_samples = np.column_stack([seeds, scaled])
        return self._param_samples

    @property
    def param_samples(self) -> np.ndarray:
        if self._param_samples is None:
            self.generate_samples()
        return self._param_samples

    def _get_params(self, random_state) -> np.ndarray:
        if random_state is None:
            return DEFAULT_PARAMS.copy()
        return self.param_samples[random_state]

    # ── Simulate recharge ─────────────────────────────────────────────────────

    def simulate_recharge(self, random_state=None) -> np.ndarray:
        """
        Simulate total recharge for one parameter realization.

        Parameters
        ----------
        random_state : int or None
            Row index into param_samples. None uses default parameters.

        Returns
        -------
        rch_ts : np.ndarray, shape (len(data.prec),), mm/d
        """
        params = self._get_params(random_state)
        return self.epikarst.simulate(self.data, params[1:8])

    # ── Simulate full ─────────────────────────────────────────────────────────

    def simulate_full(self, random_state=None, working_dir=None, cleanup=True):
        """
        Run the complete pipeline: recharge → karst network → MODFLOW + CFP.

        Parameters
        ----------
        random_state : int or None
            Row index into param_samples. None uses default parameters.
        working_dir : str, Path, or None
            Base directory for model files. Defaults to the vlab/ directory.
        cleanup : bool
            Remove pyKasso and modflow run directories after the run.

        Returns
        -------
        success : bool
        Q_out   : pd.Series  spring discharge time series
        heads   : np.ndarray shape (n_pers, n_lays, n_rows, n_cols)
        """
        params = self._get_params(random_state)
        # params: [sks_seed(0), rch_params(1-7), hk(8), k_ex(9), cad(10), diam(11)]
        rch_ts = self.epikarst.simulate(self.data, params[1:8])
        gw_params = [params[8], params[9], params[10], params[11], params[0]]

        tag = str(random_state) if random_state is not None else "default"
        proj_name = f"{self.name}_{tag}"

        base_dir = str(working_dir) if working_dir is not None else str(_VLAB_DIR)
        orig_dir = os.getcwd()

        try:
            os.chdir(base_dir)
            _flush_pykasso_loggers()

            success, Q_out, heads = self.groundwater.simulate(
                system=self,
                recharge_ts=rch_ts,
                params=gw_params,
                proj_name=proj_name,
                working_dir=base_dir,
                n_stress_periods=self.n_stress_periods,
            )
        finally:
            _flush_pykasso_loggers()
            os.chdir(orig_dir)
            if cleanup:
                for suffix in ("_pyKasso", "_modflow"):
                    d = os.path.join(base_dir, proj_name + suffix)
                    if os.path.exists(d):
                        shutil.rmtree(d)

        return success, Q_out, heads


# ── Helpers ───────────────────────────────────────────────────────────────────

def _flush_pykasso_loggers():
    """Close and remove all root logger handlers (Windows file-lock workaround)."""
    root = logging.getLogger()
    for h in root.handlers[:]:
        h.close()
        root.removeHandler(h)
