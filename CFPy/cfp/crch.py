"""
Class to create data list for the Conduit Recharge Package (CRCH) (only used
    with mode 1 or 3)

Documentation: https://pubs.usgs.gov/tm/tm6a24/    
"""

class crch():
    """
    This class handles the writing of the .crch input file for MODFLOW CFP.

    Dependencies: None

    Parameters
    ----------
    iflag_crch : is an integer value that activates or deactivates the reading
        of CRCH data. If IFLAG_CRCH is not equal to -1, NODE_NUMBERS and P_CRCH
        values are read for the total number of nodes (NNODES) in the
        simulation. Each node must be listed with NODE_NUMBERS and P_CRCH
        values. If IFLAG_CRCH equals -1, NODE_NUMBERS and P_CRCH from the last
        stress period are used for the current stress period.
    nper : number of stress periods; int
    node_numbers : is a list-like of integer values indicating the node numbers;
        list-like of ints
    p_crch : is a list-like of floats of length npipes representing the
        fraction of diffuse areal recharge (entered in the MODFLOW-2005 RCH
        Package) partitioned directly into the conduit node NODE_NUMBERS. If the
        user, for example, wants the direct conduit recharge to equal the
        diffuse recharge rate assigned for the model cell in which the pipe is
        located, the user would enter a value of 1.0 for P_CRCH. In this case,
        the diffuse areal recharge for the model cell would equal zero in
        MODFLOW water-budget calculations.
    fbc_well : is a list with tuples containing node number and specified
        pumping rate during all / different stress periods
        [(node_num, rate), (...), ...], if the pumping rate is uniform over all
        stress periods, use a list with a single float value as secont tuple
        item and if the pumping rate is non-uniform, use a list of length nper;
        list of tuples; example: [(1, [-0.15]), (5, [-0.1])] or
        [(1, [-0.1, -0.3])]

    Input Variables / Lines of the MODFLOW CFP Module (.crch file):
        0 : Comment lines.
        1 : iflag_crch - is an integer value that activates or deactivates the
            reading of CRCH data. If IFLAG_CRCH is not equal to -1, NODE_NUMBERS
            and P_CRCH values are read for the total number of nodes (NNODES) in
            the simulation. Each node must be listed with NODE_NUMBERS and
            P_CRCH values. If IFLAG_CRCH equals -1, NODE_NUMBERS and P_CRCH from
            the last stress period are used for the current stress period.
        2 : nper - number of stress periods
        3 : node_numbers - is a list of integer values indicating the node
            numbers. 
        4 : p_crch - is a list of real numbers equal to a fraction of diffuse
            areal recharge (entered in the MODFLOW-2005 RCH Package) partitioned
            directly into the conduit node NODE_NUMBERS. If the user, for
            example, wants the direct conduit recharge to equal the diffuse
            recharge rate assigned for the model cell in which the pipe is
            located, the user would enter a value of 1.0 for P_CRCH. In this
            case, the diffuse areal recharge for the model cell would equal zero
            in MODFLOW water-budget calculations.    
    """
    def __init__(
        self,
        iflag_crch,
        nper,
        node_numbers,
        p_crch=0,
        fbc_well=None
        ):
        
        self.iflag_crch = iflag_crch
        self.nper = nper
        self.node_numbers = node_numbers
        self.p_crch = p_crch
        self.fbc_well = fbc_well

        return
    
    def crch(self):
        """
        Write the .crch input file for MODFLOW-CFP

        Parameters
        ----------
        None

        Returns
        -------
        crch : the list of strings representing the contents of the .crch file
        """
        
        # initialize the crch string list
        self.crch = []
        # create line 0 string
        in0 = '# CRCH file; Stress Period '
        
        # loop over the number of stress periods
        for p in range(self.nper):
            # initialize the list holding the frationation information
            self.frac = []
            
            # check if an fbc well boundary is present
            # if not, only process fractionation information
            if self.fbc_well is None:
                # if there is no fbc_well data, only append node numbers and
                #   recharge fractions
                for i in range(len(self.node_numbers)):
                    # create the string and append it to the fractionation
                    #   information
                    self.frac.append(str(self.node_numbers[i]) + ' ' +
                                         str(self.p_crch[i]))
            # if fbc information is present, handle it
            elif self.fbc_well is not None:
                # initialize the counter
                node_iter = 0
                # loop over the node numbers
                for i in range(len(self.node_numbers)):
                    # check if there is fbc_well data corresponding
                    #   to the current node that is checked
                    if i+1 in [list(a)[0] for a in self.fbc_well]:
                        # if there is such fbc_data for the current node,
                        #   append it together with node number and recharge
                        #   fraction (case where there is a single float value
                        #   for the pumping)
                        if len(self.fbc_well[node_iter][1]) == 1:
                            # append the fractionation and well information
                            self.frac.append(str(self.node_numbers[i]) + 
                                " " + str(self.p_crch[i]) + " " +
                                str(self.fbc_well[node_iter][1]))
                        # if multiple extraction rates are given (for each
                        #   stress period), use only the value for the current
                        #   stress period
                        else:
                            # append the fractionation and well information
                            self.frac.append(str(self.node_numbers[i]) + 
                                " " + str(self.p_crch[i]) + " " +
                                str(self.fbc_well[node_iter][1][p]))
                        # increment the counter
                        node_iter += 1
                    else:
                        # if there is no fbc_well data for the current node,
                        #   only append node number and recharge fraction
                        self.frac.append(str(self.node_numbers[i]) + ' ' +
                                         str(self.p_crch[i]))
            
            # make the string such that each block of information is on a new
            #   line
            self.frac = '\n'.join(self.frac)

            # check if the information from the last stress period should be
            #   used
            if self.iflag_crch != -1:
                self.crch.append((in0 + str(p+1) + '\n' +
                                  str(self.iflag_crch) + '\n' + self.frac))
            else:
                self.crch.append(in0 + str(p+1) + '\n' + str(self.iflag_crch))
        
        return self.crch
