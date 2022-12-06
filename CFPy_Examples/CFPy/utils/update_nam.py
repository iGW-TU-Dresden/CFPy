"""
Class to update exsiting .nam file
"""

class update_nam():
    """
    Dependencies: None
        
    Input Variables:
        
        modelname: Name of the model (str)
        mode: CFP mode (int)
        unit_num: Fortran unit number (int)
        
    Input Files:
    
    Output: Updated nam file
    
    Note: Firstly, if present, any CFP related entries in the nam file will be removed. 
        Then, the nam file is updated according to the users initialization of the update_nam class. 
        
    """
    
    def __init__(self, modelname, mode, coc_unit_num=23, cfp_unit_num=16, crch_unit_num=14):
        
        self.modelname = modelname
        self.mode = mode
        self.coc_unit_num = coc_unit_num
        self.cfp_unit_num = cfp_unit_num
        self.crch_unit_num = crch_unit_num
        self.update = ['COC' + '%17s'%self.coc_unit_num + '  ' + self.modelname + '.coc' + '\n',
                       'CFP' + '%17s'%self.cfp_unit_num + '  ' + self.modelname + '.cfp' + '\n',
                       'CRCH' + '%16s'%self.crch_unit_num + '  ' + self.modelname + '.crch']
        
        if self.mode == 2:
            del(self.update[-1])

    def update_nam(self):

        ftypes = ["COC", "CFP", "CRCH"]

        # open the existing .nam-file and read the lines
        with open(self.modelname + ".nam", "r") as f:
            lines = f.readlines()

        # create new set of lines and append only those lines
        # from the original .nbr-file, which are not CFP-specific
        # (i.e., COC, CRCH, CFP)
        lines_ = []
        for line in lines:
            for ftype in ftypes:
                if ftype not in line:
                    if line in lines_:
                        continue
                    else:
                        lines_.append(line)

        # append new lines for CFP (i.e., COC, CRCH, CFP)
        for ftype in self.update:
            lines_.append(ftype)

        # open the .nam-file and write the new lines
        with open(self.modelname + ".nam", "w") as f:
            for line in lines_:
                f.write(line)
