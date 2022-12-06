"""
Class to create COC, CRCH and CFP input files from data strings
"""
import os

class write_input():
    """
    Dependencies: None
        
    Input Variables:
        
        modelname: Name of the model (str)
        data_strings: List of strings obtained from CFP classes
        file_extensions: file extensions used by CFP (str)
        
    Input Files:
    
    Output: CFP input files (.coc, .crch, .cfp)
    
    """
    
    def __init__(self, modelname, data_strings, file_extensions=['coc', 'crch', 'cfp']):
        
        self.modelname = modelname
        self.data_strings = data_strings
        self.file_extensions = file_extensions
        
    def write_input(self):

        for i in range(len(self.file_extensions)):
            try:
                os.remove(self.modelname + '.' + self.file_extensions[i])
            except FileNotFoundError:
                pass
        
        for i in range(len(self.file_extensions)):
            with open (self.modelname + '.' + self.file_extensions[i], 'w') as file:
                for string in self.data_strings[i]:
                    self.input = file.write(string + '\n')