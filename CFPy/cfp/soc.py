"""
Class to create data list for the Solute Transport Output Control (SOC)

Documentation: Unpublished, please get in touch with one of the authors
Thomas Reimann (Thomas.Reimann@tu-dresden.de)

"""

class soc():
    def __init__(
        self,
        unitnumber_trans,
        t_nts_trans,
        unitnumber_conc,
        node_numbers,
        t_nts_conc,
        unitnumber_bgt=None,
        t_nts_bgt=None,
        unitnumber_std=None,
        t_nts_std=None,
        unitnumber_smb=None,
        t_nt_smb=None,
        unitnumber_cads=None,
        t_nt_cads=None
        ):
        self.unitnumber_trans = unitnumber_trans
        self.t_nts_trans = t_nts_trans
        self.unitnumber_conc = unitnumber_conc
        self.node_numbers = node_numbers
        self.t_nts_conc = t_nts_conc
        self.unitnumber_bgt = unitnumber_bgt
        self.t_nts_bgt = t_nts_bgt
        self.unitnumber_std = unitnumber_std
        self.t_nts_std = t_nts_std
        self.unitnumber_smb = unitnumber_smb
        self.t_nt_smb = t_nt_smb
        self.unitnumber_cads = unitnumber_cads
        self.t_nt_cads = t_nt_cads

        return

    def soc(self):
        self.soc = []
        self.soc.append("#unit number for solute transport output after flow time steps")
        self.soc.append(str(self.unitnumber_trans))
        self.soc.append("#output each ... flow time step")
        self.soc.append(str(self.t_nts_trans))
        self.soc.append("#UNIT NUMBER FOR OUTPUT OF NODE CONCENTRATION AFTER TRANSPORT TIME STEPS")
        self.soc.append(str(self.unitnumber_conc))
        self.soc.append("#NUMBER OF NODES FOR OUTPUT (-1 = ALL NODES) / SUBSEQUENT LINES NODE NUMBER")
        self.soc.append(str(len(self.node_numbers)))
        for node in self.node_numbers:
            self.soc.append(str(node))
        self.soc.append("#OUTPUT EACH ... TRANSPORT TIME STEP")
        self.soc.append(str(self.t_nts_conc))
        itms = [("#UNIT NUMBER FOR OUTPUT OF STM BUDGET AFTER TRANSPORT TIME STEPS", self.unitnumber_bgt),
                ("#OUTPUT EACH ... TRANSPORT TIME STEP",self.t_nts_bgt),
                ("#UNIT NUMBER FOR SOLUTE TRANSPORT OUTPUT AFTER FLOW TIME STEPS",self.unitnumber_std),
                ("#OUTPUT EACH ... FLOW TIME STEP",self.t_nts_std),
                ("#UNIT NUMBER FOR STM NODE MASS BUDGET OUTPUT AFTER FLOW TIME STEPS",self.unitnumber_smb),
                ("#OUTPUT EACH ... FLOW TIME STEP",self.t_nt_smb),
                ("#UNIT NUMBER FOR OUTPUT OF CADS CONCENTRATION AFTER FLOW TIME STEPS",self.unitnumber_cads),
                ("#OUTPUT EACH ... FLOW TIME STEP",self.t_nt_cads)]
        for itm in itms:
            if itm[1]:
                self.soc.append(itm[0])
                self.soc.append(str(itm[1]))

        return self.soc
