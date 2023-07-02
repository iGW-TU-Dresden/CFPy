"""
Class to create COC, CRCH and CFP input files from data strings
"""
import os

class write_input():
    """
    Write the MODFLOW-CFP specific input files

    Dependencies: None
        
    Parameters
    ----------        
    modelname : name of the model; string
    data_strings : list of strings obtained from CFP classes
    file_extensions : list of strings representing file extensions used by CFP    
    """
    
    def __init__(
        self,
        modelname,
        data_strings,
        file_extensions=['coc', 'crch', 'cfp']
        ):
        
        if len(data_strings) != len(file_extensions):
            raise Exception("The length of given file extensions do not match\
                            with the length of given data strings!\n Did you\
                            forget modules?")

        self.modelname = modelname
        self.data_strings = data_strings
        self.file_extensions = file_extensions
        
    def write_input(self):
        """
        Write the input files

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        # iterate over all files
        for i in range(len(self.file_extensions)):
            # try to remove a previously existing file if it exists
            try:
                os.remove(self.modelname + '.' + self.file_extensions[i])
            # catch the possible error and do nothing in that case
            except FileNotFoundError:
                pass
        
        # iterate over all files
        for i in range(len(self.file_extensions)):
            # open the file in writing mode
            with open (self.modelname + '.' + self.file_extensions[i], 'w') as file:
                # iterate over the data strings
                for string in self.data_strings[i]:
                    # write the input to the file line by line
                    self.input = file.write(string + '\n')