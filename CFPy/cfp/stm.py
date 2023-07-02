"""
Class to create data list for the Solute Transport Simulation Package (STM)

Documentation: Unpublished, please get in touch with one of the authors
Thomas Reimann (Thomas.Reimann@tu-dresden.de)

"""

class stm():
    def __init__(
        self,
        nsections,
        nnodes,
        dr,
        conc_rmw,
        dens_rock,
        por_rock,
        kd,
        dc,
        sherwood_lam,
        sherwood_turb,
        disp,
        conc_rch={},
        conc_well={},
        conc_cauchy={},
        conc_chd={},
        conc_fs={},
        cconc={},
        cds_mfr={},
        cds_rch={}
        ):

        self.nsections = str(nsections)
        self.nnodes = str(nnodes)
        self.dr = str(dr)
        self.conc_rmw = str(conc_rmw)
        self.dens_rock = str(dens_rock)
        self.por_rock = str(por_rock)
        self.kd = str(kd)
        self.dc = str(dc)
        self.sherwood_lam = str(sherwood_lam)
        self.sherwood_turb = str(sherwood_turb)
        self.disp = str(disp)
        self.conc_rch = conc_rch
        self.conc_well = conc_well
        self.conc_cauchy = conc_cauchy
        self.conc_chd = conc_chd
        self.conc_fs = conc_fs
        self.cconc = cconc
        self.cds_mfr = cds_mfr
        self.cds_rch = cds_rch

        return

    def stm(self):
        self.stm = []
        self.stm.append("#NUMBER OF CELLS PER TUBE")
        self.stm.append(self.nsections)
        self.stm.append("#NUMBER OF NODES IN THE ROCK MATRIX")
        self.stm.append(self.nnodes)
        self.stm.append("#CELL INCREMENT OF ROCK MATRIX IN M")
        self.stm.append(self.dr)
        self.stm.append("#INITIAL CONCENTRATION IN ROCK MATRIX")
        self.stm.append(self.conc_rmw)
        self.stm.append("#DENSITY OF ROCK IN KG/L")
        self.stm.append(self.dens_rock)
        self.stm.append("#EFFECTIVE POROSITY")
        self.stm.append(self.por_rock)
        self.stm.append("#KD IN L/KG")
        self.stm.append(self.kd)
        self.stm.append("#DIFFUSION COEFFICIENT (IN WATER) IN M**2/S")
        self.stm.append(self.dc)
        self.stm.append("#SHERWOOD NUMBER IN LAMINAR FLOW (FLAG)")
        self.stm.append(self.sherwood_lam)
        self.stm.append("#SHERWOOD NUMBER IN TURBULENT FLOW (FLAG)")
        self.stm.append(self.sherwood_turb)
        self.stm.append("#DISPERSION (FLAG): DISPERSION COEFFICIENT IS CALCULATED")
        self.stm.append("# <0 => BASED ON FLOW DATA ACCORDING TO TAYLOR (1953, 1954)")
        self.stm.append("# >=0 => D_dis= VALUE * VELOCITY, I.E., VALUE=DISPERSIVITY")
        self.stm.append(self.disp)
        self.stm.append("#\n#----------------STRESS PERIOD DATA----------------")
        itms = [("#CONCENTRATION OF DIRECT RECHARGE", self.conc_rch),
                ("#CONCENTRATION OF WELL RECHARGE", self.conc_well),
                ("#CONCENTRATION OF CAUCHY BC", self.conc_cauchy),
                ("#CONCENTRATION OF FIXED HEAD / FHLQ", self.conc_chd),
                ("#CONCENTRATION OF INFLOW FROM FISSURED SYSTEM", self.conc_fs),
                ("#FIXED CONCENTRATION NODES?", self.cconc),
                ("#MASS FLOW RATE CADS", self.cds_mfr),
                ("#CADS RECHARGE CONCENTRATION", self.cds_rch)]
        # Restructure stress period data
        self.spd = {}
        for itm in itms:
            for key, val in itm[1].items():
                add = [itm[0]]
                if type(val) is list:
                    for node, node_val in zip(val[0], val[1]):
                        add.append(f"{node} {node_val}")
                else:
                    add.append(str(val))
                if key in self.spd.keys():
                    self.spd[key] += add
                else:
                    self.spd[key] = add

        for stressperiod in sorted(self.spd.keys()):
            self.stm += self.spd[stressperiod]
            self.stm.append("#STRESS PERIOD DATA VALID UNTIL STRESS PERIOD")
            self.stm.append(str(stressperiod+1))# !TR: 2023 06 29 THIS TO ACCOUNT PYTHON START COUNTING ZERO-BASED

        return self.stm
