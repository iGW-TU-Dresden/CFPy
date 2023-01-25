"""
Class to create nbr data for the Conduit Flow Package (CFP)

Documentation: https://pubs.usgs.gov/tm/tm6a24/    
"""
import os
import sys
import time

class nbr():
    """
    Dependencies: None

    Input Variables:
        
        0: .nbr file containing conduit elevations as positive real numbers; 
            diffuse-recharge-only cells have a negative value (real)

    Input Files: .nbr file
    
    Output: List of strings.
    
    """
    def __init__(self):
        
        self.bot_elev = []
        self.bot_elev2 = []
        self.cond_elev = []
        self.lines = []
        
    def nbr_read(self):
        
        tic1 = time.time()
        
        # check if line list is empty - if not, make a new empty list
        if len(self.lines) > 0:
            self.lines = []

        # Search and list .nbr files found in the working directory.
        self.nbr_file = [file for file in os.listdir('.') if file.endswith('.nbr')]
        
        read_count = 0
        
        # If a single .nbr file is found, read the file and append each line that does not begin with '#' to a list.
        if len(self.nbr_file) == 1:
            with open(self.nbr_file[0], 'r') as file:
                for line in file:
                    if line[0] == '#':
                        # count comment lines
                        read_count += 1
                    else:
                        self.lines.append(line.split())
        
        # Exit, if no .nbr file is found.
        elif len(self.nbr_file) == 0:
                sys.exit('No .nbr file found.')
        
        # Exit, if more than one .nbr file is found.
        elif len(self.nbr_file) > 1:
                sys.exit('More than one .nbr file found.')
        
        # If a line contains a single value, add it to 'bot_elev'.
        # NOTE: for this to work, one has to specify BOTH values in line 2
        #     of the .nbr-file (number of layers, NUMBER OF NODE PLANES)
        #     otherwise, line 2 has only one character and gets appended,
        #     which is not what we want
        for line in self.lines:
            # we only have a single value in a line if it is a spatially
            #   constant elevation
            if len(line) == 1:
                self.bot_elev.append(line[0])
                read_count -= 1
        # the second line (index 1 ) always contains the number of layers
        #   and the number of node planes
        nlays = int(self.lines[1][0])
        nplanes = int(self.lines[1][1])

        # Remove top/bottom elevations that were added to 'bot_elev' from 'lines'.
        # remove already used lines from list of lines
        # PROBLEM: if model dimensions and number of layers is given, this leads
        #     to problems. In this case, len(self.bot_elev) + 2 lines would have
        #     to be removed
        # --> assume that the user specifies the first two lines in .nbr-file
        #     (equal to CONGEN input)
        if len(self.bot_elev) > 0:
            # + 2 because first two lines contain information as well
            #     but are not used
            del self.lines[0:len(self.bot_elev) + 2]

        # NEW
        # if bot_elev does NOT contain values (i.e. no single values as layer elevations are used)
        # delete the first two lines (nrows / ncols & nlays / nplanes)
        else:
            del self.lines[0:2]
        
        # Split 'lines' into lists of equal sizes; total number of lists is defined by 'read_count'.        
        self.lines = [self.lines[i:i+int(len(self.lines)/read_count)] for i in range(0, len(self.lines),
            int(len(self.lines)/read_count))]

        # If list 'bot_elev' has already been defined and removed from 'lines', 'cond_elev' = 'lines'.
        if len(self.bot_elev) > 0:
            self.cond_elev = self.lines
        
        # If list 'bot_elev' has not yet been defined (i.e., top/bottom elevations are not given as single values),
        # append each list element in 'lines' that contains a negative value (i.e., inactive conduit cell) 
        # or produces a ValueError upon check (e.g., for values proceeded by a character) to list 'cond_elev'.
        # Break from iteration, if a condition becomes True.
        else:
            # for lay in range(len(lines)):
            #     for row in range(len(lines[lay])):
            #         for col in range(len(lines[lay][row])):
            #             try:
            #                 if int(lines[lay][row][col]) == -999: # change from < 0
            #                     self.cond_elev.append(lines[lay])
            #                     break
            #             except ValueError:
            #                 self.cond_elev.append(lines[lay])
            #                 break
            #             break
            #         break

            # use number of layers and number of planes information to fill cond_elev and bot_elev
            self.cond_elev = self.lines[-nplanes:]

        # Define bot_elev as all elements of list 'lines' that have not previously been added to 'cond_elev'.
        if len(self.bot_elev) == 0:
            self.bot_elev = self.lines[0:len(self.lines)-len(self.cond_elev)]
            
        # If single values are given as bottom elevation, extent the list element to the size of 'cond_elev'
        if isinstance(self.bot_elev[0], list) == False:
            for lay in range(len(self.bot_elev)):
                for row in range(len(self.cond_elev[0])):
                    self.bot_elev2.append([self.bot_elev[lay]] * len(self.cond_elev[0][0]))
        
        if len(self.bot_elev2) > 0:
            self.bot_elev2 = [self.bot_elev2[i:i+int(len(self.bot_elev2)/len(self.bot_elev))] for i in range(0,len(self.bot_elev2), int(len(self.bot_elev2)/len(self.bot_elev)))]
            
        toc1 = time.time()
        print('Elapsed time (.nbr file read): ' + str(round(toc1-tic1,2)) + ' s')
        
        if len(self.bot_elev2) > 0:
            return self.bot_elev2, self.cond_elev
        else:
            return self.bot_elev, self.cond_elev
        
    def nbr(self, bot_elev, cond_elev):
        
        tic2 = time.time()
        
        self.bot_elev = bot_elev
        self.cond_elev = cond_elev
        self.node_numbers = []
        self.node_loc = []
        self.cond_loc = []
        self.node_nbr = []
        self.plane_number = []
        self.tube_numbers = []
        self.tube_pair = []
        self.tube_nbr = []
            
        node = 0
        tube = 0
        
        #Search for conduit elevations and create (i) a list of node_numbers and (ii) a list of their respective location, node_loc, in the model grid
        # iterate over layers, rows, columns and look for value in cond_elev:
        #     if no vertical connection (no leading "c"), find layer in which the node lies in;
        #     if vertical connection, disregard "c" and find layer;
        #     append locations to node_loc, cond_loc, append plane number to plane_number
        for lay in range(len(self.cond_elev)):
            for row in range(len(self.cond_elev[lay])):
                for col in range(len(self.cond_elev[lay][row])):
                    try:
                        if round(float(self.cond_elev[lay][row][col]),1) >= 0:
                            node += 1
                            self.node_numbers.append(node)

                            for layer in range(len(self.bot_elev)):
                                if round(float(self.cond_elev[lay][row][col]),1) <= round(float(self.bot_elev[layer][row][col]),1):
                                    z = layer
                            
                            # append to cond_loc similarly as for the "except"-case below
                            self.cond_loc.append([col+1, row+1, z+1])    
                            self.node_loc.append([col+1, row+1, z+1])
                            self.plane_number.append([lay+1])
                    
                    except ValueError:
                        if self.cond_elev[lay][row][col][0] == 'c':
                            node += 1
                            self.node_numbers.append(node)
                            
                            for layer in range(len(self.bot_elev)):
                                if round(float(self.cond_elev[lay][row][col][1:]),1) <= round(float(self.bot_elev[layer][row][col]),1):
                                    z = layer
                            
                            self.cond_loc.append([col+1, row+1, z+1])
                            self.node_loc.append([col+1, row+1, z+1])
                            self.plane_number.append([lay+1])
                                
        for i in range(len(self.node_numbers)):
            self.node_nbr.append([])
        
        #Check for adjacent nodes; clockwise in the plane for all locations given by node_loc and plane_number
        #Append neighbour nodes to node_nbr, tube connections to tube_pair and the respective tube to tube_numbers
        # iterate over nodes, nodes again (checking if neighbor);

        # iterate over nodes
        for i in range(len(self.node_loc)):
            # iterate over nodes
            for j in range(len(self.node_loc)):
                # compare (column, row - 1, plane (= location above)) information of original node i with
                #     (column, row, plane) information of node j;
                #     if equal (node j is node above node i): append j to node_nbr of i and
                #     check if connection of nodes (i, j) or (j, i) are already in tube_pair;
                #     if not: append connection and tube number to tube_pair and tube_numbers
                if [self.node_loc[i][0],self.node_loc[i][1]-1,self.plane_number[i]] == [self.node_loc[j][0],self.node_loc[j][1],self.plane_number[j]]:
                    self.node_nbr[i].append(self.node_numbers[j])
                    if [self.node_numbers[i],self.node_numbers[j]] not in self.tube_pair:
                        if [self.node_numbers[j],self.node_numbers[i]] not in self.tube_pair:
                            tube+=1
                            self.tube_pair.append([self.node_numbers[i],self.node_numbers[j]])
                            self.tube_numbers.append(tube)
            
            # same procedure, check for node to the right
            for j in range(len(self.node_loc)):
                if [self.node_loc[i][0]+1,self.node_loc[i][1],self.plane_number[i]] == [self.node_loc[j][0],self.node_loc[j][1],self.plane_number[j]]:
                    self.node_nbr[i].append(self.node_numbers[j])
                    if [self.node_numbers[i],self.node_numbers[j]] not in self.tube_pair:
                        if [self.node_numbers[j],self.node_numbers[i]] not in self.tube_pair:
                            tube+=1
                            self.tube_pair.append([self.node_numbers[i],self.node_numbers[j]])
                            self.tube_numbers.append(tube)

            # same procedure, check for node below
            for j in range(len(self.node_loc)):
                if [self.node_loc[i][0],self.node_loc[i][1]+1,self.plane_number[i]] == [self.node_loc[j][0],self.node_loc[j][1],self.plane_number[j]]:
                    self.node_nbr[i].append(self.node_numbers[j])
                    if [self.node_numbers[i],self.node_numbers[j]] not in self.tube_pair:
                        if [self.node_numbers[j],self.node_numbers[i]] not in self.tube_pair:
                            tube+=1
                            self.tube_pair.append([self.node_numbers[i],self.node_numbers[j]])
                            self.tube_numbers.append(tube)

            # same procedure, check for node to the left
            for j in range(len(self.node_loc)):
                if [self.node_loc[i][0]-1,self.node_loc[i][1],self.plane_number[i]] == [self.node_loc[j][0],self.node_loc[j][1],self.plane_number[j]]:
                    self.node_nbr[i].append(self.node_numbers[j])
                    if [self.node_numbers[i],self.node_numbers[j]] not in self.tube_pair:
                        if [self.node_numbers[j],self.node_numbers[i]] not in self.tube_pair:
                            tube+=1
                            self.tube_pair.append([self.node_numbers[i],self.node_numbers[j]])
                            self.tube_numbers.append(tube)

        #Check for adjacent nodes; up and down from locations given by cond_loc
        #Append neighbour nodes to node_nbr, tube connections to tube_pair and the respective tube to tube_numbers

        ##############################################################################
        # In the following two loops there is the difference to the "old" version:
        # contrary to the old version, here it now gets correctly checked if two nodes
        #   are vertically connected (case happens if one has two node planes which
        #   are vertically connected at one or multiple locations by a leading "c")
        ##############################################################################

        # similar as before but compare ((col, row), plane)
        # check node above (up)
        for i in range(len(self.cond_loc)):
            for j in range(len(self.cond_loc)):
                # original check for plane number:
                # self.plane_number[self.node_loc.index(self.cond_loc[i])][0]

                if [self.cond_loc[i][0:2], self.plane_number[i][0]-1] == [self.cond_loc[j][0:2], self.plane_number[j][0]]:

                    # check if nodes i and j both have a leading "c" in
                    #     their corresponding cond_elev entry;
                    #     if not, don't create connection
                    if self.cond_elev[self.plane_number[i][0] - 1][self.cond_loc[i][1] - 1][self.cond_loc[i][0] - 1][0] == "c" and self.cond_elev[self.plane_number[j][0] - 1][self.cond_loc[j][1] - 1][self.cond_loc[j][0] - 1][0] == "c":
                        # self.node_nbr[self.node_loc.index(self.cond_loc[i])].append(self.node_numbers[self.node_loc.index(self.cond_loc[j])])
                        self.node_nbr[i].append(self.node_numbers[j])
                    
                        if [self.node_numbers[i],self.node_numbers[j]] not in self.tube_pair:
                            if [self.node_numbers[j],self.node_numbers[i]] not in self.tube_pair:
                                tube+=1
                                self.tube_pair.append([self.node_numbers[i],self.node_numbers[j]])
                                self.tube_numbers.append(tube)

            # check node below (down)
            for j in range(len(self.cond_loc)): 
                if [self.cond_loc[i][0:2], self.plane_number[i][0]+1] == [self.cond_loc[j][0:2], self.plane_number[j][0]]:
                    # check if nodes i and j both have a leading "c" in
                    #     their corresponding cond_elev entry;
                    #     if not, don't create connection
                    if self.cond_elev[self.plane_number[i][0] - 1][self.cond_loc[i][1] - 1][self.cond_loc[i][0] - 1][0] == "c" and self.cond_elev[self.plane_number[j][0] - 1][self.cond_loc[j][1] - 1][self.cond_loc[j][0] - 1][0] == "c":
                        # self.node_nbr[self.node_loc.index(self.cond_loc[i])].append(self.node_numbers[self.node_loc.index(self.cond_loc[j])])
                        self.node_nbr[i].append(self.node_numbers[j])
                    
                        if [self.node_numbers[i],self.node_numbers[j]] not in self.tube_pair:
                            if [self.node_numbers[j],self.node_numbers[i]] not in self.tube_pair:
                                tube+=1
                                self.tube_pair.append([self.node_numbers[i],self.node_numbers[j]])
                                self.tube_numbers.append(tube)

        
        #Add 0 to each lsit in node_nbr until length is 6
        for nbr in self.node_nbr:
            while len(nbr) < 6:
                nbr.append(0)
        
        for i in range(len(self.node_numbers)):
            self.tube_nbr.append([])
        
        #Add tube_numbers to tube_nbr
        for i in range(len(self.node_numbers)):
            for j in range(len(self.node_nbr[i])):
                if [self.node_numbers[i], self.node_nbr[i][j]] in self.tube_pair:
                    self.tube_nbr[i].append(self.tube_numbers[self.tube_pair.index([self.node_numbers[i], self.node_nbr[i][j]])])
                elif [self.node_nbr[i][j], self.node_numbers[i]] in self.tube_pair:
                    self.tube_nbr[i].append(self.tube_numbers[self.tube_pair.index([self.node_nbr[i][j], self.node_numbers[i]])])
        
        #Add 0 to each lsit in node_nbr until length is 6
        for nbr in self.tube_nbr:
            while len(nbr) < 6:
                nbr.append(0)
                
        toc2 = time.time()
        print('Elapsed time (write nbr data): ' + str(round(toc2-tic2,2)) + ' s')

        return self.node_numbers, self.plane_number, self.node_loc, self.cond_loc, self.node_nbr, self.tube_numbers, self.tube_pair, self.tube_nbr
        
