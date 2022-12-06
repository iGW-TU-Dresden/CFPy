"""
Class to create data list for the Conduit Recharge Package (CRCH) (only used with mode 1 or 3)

Documentation: https://pubs.usgs.gov/tm/tm6a24/    
"""

class crch():
    """
    Dependencies: None

    Input Variables of the MODFLOW CFP CRCH Module:
        
        0: Comment lines.
        1: iflag_crch - is an integer value that activates or deactivates the reading of CRCH data.
            If IFLAG_CRCH is not equal to -1, NODE_NUMBERS and P_CRCH values are read for the total number of nodes (NNODES) in the simulation. 
            Each node must be listed with NODE_NUMBERS and P_CRCH values.
            If IFLAG_CRCH equals -1, NODE_NUMBERS and P_CRCH from the last stress period are used for the current stress period.
        2: nper - number of stress periods
        3: node_numbers - is a list of integer values indicating the node numbers. 
        4: p_crch - is a list of real numbers equal to a fraction of diffuse areal recharge (entered in the MODFLOW-2005 RCH Package) partitioned directly into the conduit node NODE_NUMBERS. 
            If the user, for example, wants the direct conduit recharge to equal the diffuse recharge rate assigned for the model cell in which the pipe is located, the user would enter a value of 1.0 for P_CRCH. 
            In this case, the diffuse areal recharge for the model cell would equal zero in MODFLOW water-budget calculations.

        # INPUT VARIABLES of the CFPy.cfp.crch module
        # iflag_crch : is an integer value that activates or deactivates the reading of CRCH data.
            If IFLAG_CRCH is not equal to -1, NODE_NUMBERS and P_CRCH values are read for the total number of nodes (NNODES) in the simulation. 
            Each node must be listed with NODE_NUMBERS and P_CRCH values.
            If IFLAG_CRCH equals -1, NODE_NUMBERS and P_CRCH from the last stress period are used for the current stress period.
        # nper : number of stress periods, int
        # node_numbers : is a list of integer values indicating the node numbers
        # p_crch : is a list of floats equal to a fraction of diffuse areal recharge (entered in the MODFLOW-2005 RCH Package) partitioned directly into the conduit node NODE_NUMBERS. 
            If the user, for example, wants the direct conduit recharge to equal the diffuse recharge rate assigned for the model cell in which the pipe is located, the user would enter a value of 1.0 for P_CRCH. 
            In this case, the diffuse areal recharge for the model cell would equal zero in MODFLOW water-budget calculations
        # fbc_well : is a list with tuples containing node number and specified pumping
            rate during all stress periods [(node_num, rate), (...), ...]; list of tuples
            example: [(1, -0.15), (5, -0.1)]

    Input Files:
    
    Output: List of strings.
    
    """
    def __init__(self, iflag_crch, nper, node_numbers, p_crch=0, fbc_well=None):
        
        self.iflag_crch = iflag_crch
        self.nper = nper
        self.node_numbers = node_numbers
        self.p_crch = p_crch
        self.fbc_well = fbc_well
        
        return
    
    def crch(self):
        
        in0 = '# CRCH file; Stress Period '
        
        self.frac = []
        
        if self.fbc_well is None:
            # if there is no fbc_well data, only append node numbers and
            #   recharge fractions
            for i in range(len(self.node_numbers)):
                self.frac.append(str(self.node_numbers[i]) + ' ' + str(self.p_crch[i]))
        elif self.fbc_well is not None:
            # if there is fbc_well data, handle it
            node_iter = 0
            for i in range(len(self.node_numbers)):
                # check if there is fbc_well data corresponding
                #   to the current node that is checked
                if i+1 in [list(a)[0] for a in self.fbc_well]:
                    # if there is such fbc_data for the current node,
                    #   append it together with node number and recharge
                    #   fraction
                    self.frac.append(str(self.node_numbers[i]) + 
                        " " + str(self.p_crch[i]) + " " +
                        str(self.fbc_well[node_iter][1]))
                    node_iter += 1
                else:
                    # if there is no fbc_well data for the current node,
                    #   only append node number and recharge fraction
                    self.frac.append(str(self.node_numbers[i]) + ' ' + str(self.p_crch[i]))
        
        self.frac = '\n'.join(self.frac)
        
        self.crch = []
        
        for i in range(self.nper):
            if self.iflag_crch != -1:
                self.crch.append(in0 + str(i+1) + '\n' + str(self.iflag_crch) + '\n' + self.frac)
            else:
                self.crch.append(in0 + str(i+1) + '\n' + str(self.iflag_crch))
        
        return self.crch
