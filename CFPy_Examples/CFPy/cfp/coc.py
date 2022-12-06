"""
Class to create data list for the Conduit Output Control Option (COC) (only used with mdoe 1 or 3)

Documentation: https://pubs.usgs.gov/tm/tm6a24/    
"""

class coc():
    """
    Dependencies: None

    Input Variables:
        
        0, 1, 3, 5, 7, 9, 11: Comment lines.
        2: nnodes - is an integer number equal to the number of nodes for which flow and head values are desired in separate ouput files.
        4: node_numbers - are integer values of the node numbers for which flow and head values are desired in separate output files. List one node number per line.
        6: n_nts - is an integer value equal to the time step interval for output of node head and flow values. As an example, N_NTS equal to 2 activates node output at each NODE_NUMBERS every 2 time steps.
        8: npipes - is an integer value equal to the number of pipes for which flow rates and Reynolds numbers are desired in separate output files.
        10: pipe_numbers - are integer values of the pipe numbers for which flow and head values are desired in separate output files. List one pipe number per line.
        12: t_nts - is an integer value equal to the time step interval for output of pipe flow rates and Reynolds num- bers. As an example, T_NTS equal to 2 activates node output at each PIPE_NUMBERS every 2 time steps.
    
    Input Files:
    
    Output: List of strings.
    
    """
    def __init__(self, nnodes, node_numbers, n_nts, npipes, pipe_numbers, t_nts):
        
        self.nnodes = str(nnodes)
        self.node_numbers = '\n'.join([str(nn) for nn in node_numbers])
        self.n_nts = str(n_nts)
        self.npipes = str(npipes)
        self.pipe_numbers = '\n'.join([str(pn) for pn in pipe_numbers])
        self.t_nts = str(t_nts)
        
        return
    
    def coc(self):
        
        in0 = '# COC file'
        in1 = '# Number of nodes for output (nnodes)'
        in3 = '# Node numbers, one per line (node_numbers)'
        in5 = '# Node output each n time steps (n_nts)'
        in7 = '# Number of conduits for output (npipes)'
        in9 = '# Conduit numbers, one per line (pipe_numbers)'
        in11 = '# Conduit output each n time steps (t_nts)'
        
        self.coc = [in0, in1, self.nnodes, in3, self.node_numbers, in5, self.n_nts, in7, self.npipes, in9, self.pipe_numbers, in11, self.t_nts]
        
        return self.coc
