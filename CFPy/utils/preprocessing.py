import numpy as np
import os

"""
TODO: at the moment, only 2D information can be processed. Additionally,
if multiple 2D node networks are generated and vertically stacked in the
.nbr-file, the vertical connections still need to be defined manually.
Develop a utility to automatically connect the nodes vertically (e.g.,
if they are at the same ROW, COL and in adjacent layers, add a vertical
connection).
"""

class Preprocessor():

	"""
	A class for pre-processing node network input information.
	It is meant to process, e.g., existing network data from
	pyKasso or similar packages.

	NOTE: For now, only 2D arrays are supported as input data
	types as pyKasso can only create 2D structures.
	If the GeneralValidator is used, use

	network = np.array(catchment.karst_simulations[-1].maps['karst'][-1])
	network = np.flip(network, 0)

	to obtain a data structure that can be processed with this
	class.
	If the pyKasso-specific validator (pyKassoValidator) is used,
	the network argument should be a pyKasso.SKS()-catchment instance. From
	that, all data is extracted automatically.
	"""

	def __init__(self, network=None, elevations=None, flopymodel=None):
		"""
		Arguments
		---------
		network : if the GeneralValidator will be used: a 2D numpy array of shape
			(n_rows, n_cols) representing the node network, numpy ndarray;
			if the pyKassoValidator will be used: a pyKasso.SKS()-catchment
			instance
		elevations : a 2D numpy array of shape (n_rows, n_cols) representing
			the node elevations; note that it may happen that nodes are added
			to the network during pre-processing - give node elevations at
			each cell, if possible, numpy ndarray
		flopymodel : specifying whether CFPy will be used on a pre-existing
			FloPy model (flopymodel=FloPy.modflow.Modflow model object) or if
			a new FloPy model is created together with CFPy (None); default is
			None, flopy.modflow.Modflow onject or None

		NOTE: if a flopy.modflow.Modflow object is given, all information
		from the flopy model used by CFPy (i.e., spatial discretization, layer
		elevations) are taken from the existing model automatically. With that,
		the user does not have to specify this information in CFPy manually
		"""

		self.network = network
		self.elevations = elevations
		self.flopymodel = flopymodel

		# if a flopy model is given, extract necessary data
		if self.flopymodel is not None:
			self.nrows = self.flopymodel.nrow_ncol_nlay_nper[0]
			self.ncols = self.flopymodel.nrow_ncol_nlay_nper[1]
			self.nlays = self.flopymodel.nrow_ncol_nlay_nper[2]
			self.layer_elevations = self.flopymodel.modelgrid.top_botm

	def export_network(self, path=None):
		
		"""
		Export the validated network to a separate text file given by the path argument.
		The contents of the text file can be directly included in the .nbr-file, which
		is used by CFPy at later stages.

		Arguments
		---------
		path : the path to the directory in which the file should be saved

		NOTE: during execution, the network array is already structured to fit the
		requirements of CFPy: zeros (0.) are replaced with -999. Take into account that
		nodes, which have a height of 0. will be changed as well! If this is not what you
		want, try giving these nodes a value such as, e.g.,  1e-5 or 1e-10.
		"""

		if path is None:
			path = os.getcwd()
			path = os.path.join(path, "CFPy_exported_network.txt")

		# copy network
		network_ = self.network.copy()
		
		# exporting to file
		with open(path, "w") as file:
			file.write("# node network structure for the .nbr-file\n")
			for i in network_:
				for j in i:
					if j == 0.:
						file.write("-999 ")
					else:
						file.write(str(j) + " ")
				file.write("\n")

		return

	def generate_nbr(self, path, nplanes, nrows=None, ncols=None, nlays=None, layer_elevations=None):
		
		"""
		Generate the .nbr-file for CFPy from the given information and the validated
		network and save it in the directory given by path.

		Arguments
		---------
		path : the directory to save the generated file in, str
		nrows : number of rows in the model domain, int
		ncols : number of cols in the model domain, int
		nlays : number of MODFLOW layers in the model domain, int
		nplanes : number of node planes in the model domain, int
		layer_elevations : (nlays + 1) MODFLOW layer elevations
			(model top + n_lays) layer information data; the elevation
			data can either be a single float (specifying uniform layer
			elevations) or arrays of shape (nrows, ncols) specifying
			cell-by-cell layer elevations; in all cases, the information
			sould be given as a list-like, list of ints or numpy array

			examples: [[12.4], [5.3]] for uniform top and bottom elevations
			or [[[5.0, 5.0], [5.0, 5.0]], [[1.0, 1.0], [1.0, 1.0]]] for
			cell-by-cell top and bottom elevations

		Returns
		-------
		"""

		# use pre-existing flopy model data if available
		if self.flopymodel is not None:
			nrows = self.nrows
			ncols = self.ncols
			nlays = self.nlays
			layer_elevations = self.layer_elevations

		if (nrows is None or ncols is None or nlays is None or layer_elevations is None):
			raise ValueError("Not all necessary input data (nrows, ncols, nlays, layer_elevations) are given!")

		if path is not None and os.path.isdir(path):
			path = os.path.join(path, "CFPy_nbr.nbr")
		if path is None:
			path = os.getcwd()
			path = os.path.join(path, "CFPy_nbr.nbr")

		if self.valid == True:
			with open(path, "w+") as f:
				# write nrows, ncols
				f.write(str(nrows) + " " + str(ncols) + "\n")
				# write nlays, nplanes
				f.write(str(nlays) + " " + str(1) + "\n") # the 1 here defines the number of node planes
				# comment line
				f.write("#\n")
				# write layer elevations
				for l in layer_elevations:
					if len(l) == 1:
						f.write("{0:4.2f} ".format(l[0]))
						f.write("\n")
					else:
						for i in l:
							for j in i:
								f.write("{0:4.2f} ".format(j))
							f.write("\n")
					f.write("#" + "\n")
				# write node plane
				for i in self.network:
					for j in i:
						if j == 0.:
							f.write("-999 ")
						else:
							f.write("{0:4.2f} ".format(j))
					f.write("\n")

		else:
			raise ValueError("The network has not been validated yet. Validate the network first!")

class GeneralValidator(Preprocessor):

	"""
	A class for validating the given 2D node network to ensure that further
	computations give the desired results.
	"""

	def __init__(self, network, elevations, flopymodel=None):

		"""
		network : a 2D numpy array of shape (n_rows, n_cols) representing the
			node network, numpy ndarray
		elevations : a 2D numpy array of shape (n_rows, n_cols) representing
			the node elevations; note that it may happen that nodes are added
			to the network during pre-processing - give node elevations at
			each cell, if possible, numpy ndarray
		"""

		Preprocessor.__init__(
			self,
			network=network,
			elevations=elevations,
			flopymodel=flopymodel
			)
		# set the valid attribute to False initially; if the network gets validated
		# later, the value is set to True
		self.valid = False

	def validate_network(self):

		"""
		Validate a given 2D node network. Because diagonally adjecant nodes will not be connected
		to conduits in CFPy, networks can be pre-processed accordingly by adding the
		neccessary nodes to close "gaps".
		If a 3D node network should be validated, use the validate_network method for each
		individual 2D node plane separately.

		Returns
		-------
			network : the validated node network. If elevations are given, the validated
			network nodes are at the defined elevations, numpy ndarray

		NOTE: If multiple isolated conduit branches are present in the network, it may happen
		that the CFPy pre-processing (validate_network method) connects the branches if they
		are close (when the gap between nodes is not more than one cell).
		Always visually inspect the resulting node network afterwards or ensure that such
		behavior does not affect the results!
		"""

		if self.network is None:
			raise ValueError("No network data is given!")
		if self.elevations is None:
			print("Node elevations are not given! The remaining calculations are now carried out.")
		print("\nAlways visually check the validated network for structural correctness! \n\
			i.e., whether branches are correctly isolated or if they got connected during processing.")

		# check whether network and layer_elevations have the same number of rows and columns
		if self.flopymodel is not None:
			if self.network.shape[0] != self.nrows:
				raise ValueError(f"The network has {self.network.shape[0]} rows but the model has {self.nrows} rows!")
			if self.network.shape[1] != self.ncols:
				raise ValueError(f"The network has {self.network.shape[1]} columns but the model has {self.cols} columns!")

		# iterate over the network
		for i in range(self.network.shape[0]):
			# iterate over columns
			for j in range(self.network.shape[1]):
				# check if node is present
				if self.network[i, j] != 1.:
					continue
				else:
					# check that cell is not on boundary
					if (i > 0 and
						j > 0 and
						i < self.network.shape[0] - 1 and
						j < self.network.shape[1] - 1):
						
						# check if node has diagonal neighbor at all:
						if (self.network[i-1, j-1] == 1. or
							self.network[i-1, j+1] == 1. or
							self.network[i+1, j+1] == 1. or
							self.network[i+1, j-1] == 1.):
							
							# check if a node is present at (ROW-1, COL-1) (up-left diagonal neighbor)
							if (self.network[i-1, j-1] == 1.):
								# don't do anything if a node is already present as a upper / sideways
								# neighbor
								if (self.network[i-1, j] == 1. or self.network[i, j-1] == 1.):
									pass
								# randomly add a node either as a upper or sideways neighbor 
								else:
									if (np.random.random() >= 0.5):
										self.network[i-1, j] = 1.
									else:
										self.network[i, j-1] = 1.
							
							# check if a node is present at (ROW-1, COL+1) (up-right diagonal neighbor)
							if (self.network[i-1, j+1] == 1.):
								# don't do anything if a node is already present as a upper / sideways
								# neighbor
								if (self.network[i-1, j] == 1. or self.network[i, j+1] == 1.):
									pass
								# randomly add a node either as a upper or sideways neighbor 
								else:
									if (np.random.random() >= 0.5):
										self.network[i-1, j] = 1.
									else:
										self.network[i, j+1] = 1.
							# check if a node is present at (ROW+1, COL+1) (down-right diagonal neighbor)
							if (self.network[i+1, j+1] == 1.):
								# don't do anything if a node is already present as a upper / sideways
								# neighbor
								if (self.network[i+1, j] == 1. or self.network[i, j+1] == 1.):
									pass
								# randomly add a node either as a upper or sideways neighbor 
								else:
									if (np.random.random() >= 0.5):
										self.network[i+1, j] = 1.
									else:
										self.network[i, j+1] = 1.
							
							# check if a node is present at (ROW+1, COL-1) (down-left diagonal neighbor)	
							if (self.network[i+1, j-1] == 1.):
								# don't do anything if a node is already present as a upper / sideways
								# neighbor
								if (self.network[i+1, j] == 1. or self.network[i, j-1] == 1.):
									pass
								# randomly add a node either as a upper or sideways neighbor 
								else:
									if (np.random.random() >= 0.5):
										self.network[i+1, j] = 1.
									else:
										self.network[i, j-1] = 1.
						
						# if no diagonal neighbor is found: don't do anything and proceed to the next node
						else:
							continue

		# if elevations are given, change node heights to corresponding elevations
		if self.elevations is not None:
			self.network *= self.elevations
		
		# set validity of the network to True
		self.valid = True
		return self.network

class pyKassoValidator(Preprocessor):

	"""
	A class for validating the given 2D node network calculated by pyKasso to ensure that further
	computations give the desired results.

	If this class is used, the network argument should be given as a pyKasso.SKS()-catchment
	instance!
	"""

	def __init__(self, network, elevations, flopymodel=None, sim_num=-1):

		"""
		network : a pyKasso.SKS()-catchment instance
		elevations : a 2D numpy array of shape (n_rows, n_cols) representing
			the node elevations; note that it may happen that nodes are added
			to the network during pre-processing - give node elevations at
			each cell, if possible, numpy ndarray
		sim_num : the simulation number for the pyKasso.SKS()-catchment
			instance, defining the karst network which should be used (default
			is -1, using the last iteration), int
		"""

		Preprocessor.__init__(
			self,
			network=network,
			elevations=elevations,
			flopymodel=flopymodel
			)

		self.catchment = network
		self.sim_num = sim_num
		self.valid = False

	def validate_network(self):

		"""
		Validate a given 2D node network. Because diagonally adjecant nodes will not be connected
		to conduits in CFPy, networks can be pre-processed accordingly by adding the
		neccessary nodes to close "gaps".
		Contrary to the GeneralValidator, the pyKassoValidator also takes the inlet and outlet
		information into account and ensures that nodes are present at these locations.

		Returns
		-------
			network : the validated node network. If elevations are given, the validated
			network nodes are at the defined elevations, numpy ndarray

		NOTE: If multiple isolated conduit branches are present in the network, it may happen
		that the CFPy pre-processing (validate_network method) connects the branches if they
		are close (when the gap between nodes is not more than one cell).
		Always visually inspect the resulting node network afterwards or ensure that such
		behavior does not affect the results!

		NOTE: matrix cells are defined as 0.0 - avoid node elevations of 0.0 in the
		elevations argumemt! Having 0.0s in the elevations data structure can lead
		to unwanted behavior, because these locations would be treated as metrix cells!
		"""

		if self.network is None:
			raise ValueError("No network data (catchment instance) is given!")
		if self.elevations is None:
			print("Node elevations are not given! The remaining calculations are now carried out.")
		print("\nAlways visually check the validated network for structural correctness! \n\
			i.e., whether branches are correctly isolated or if they got connected during processing.")

		# get the network from the pyKasso.SKS-catchment as array
		self.network = np.array(self.network.karst_simulations[self.sim_num].maps['karst'][-1])

		# get inlet locations as cell indices from (x, y)-information
		inlet_locs = []
		for i in self.catchment.inlets:
			loc_x = int(np.floor(i[0] / self.catchment.settings["dx"]))
			loc_y = int(np.floor(i[1] / self.catchment.settings["dy"]))
			inlet_locs.append([loc_y, loc_x])

		# get outlet locations as cell indices from (x, y)-information
		outlet_locs = []
		for i in self.catchment.outlets:
			loc_x = int(np.floor(i[0] / self.catchment.settings["dx"]))
			loc_y = int(np.floor(i[1] / self.catchment.settings["dy"]))
			outlet_locs.append([loc_y, loc_x])

		# check presence of inlet locations in the node network
		# if no node is present at the inlet location, add a node
		for i in inlet_locs:
			if self.network[i[0], i[1]] == 0.:
				self.network[i[0], i[1]] = 1.

		# check presence of outlet locations in the node network
		# if no node is present at the outlet location, add a node
		for i in outlet_locs:
			if self.network[i[0], i[1]] == 0.:
				self.network[i[0], i[1]] = 1.

		# flip the network such that (ROW, COL)-indices are in MODFLOW format
		self.network = np.flip(self.network, 0)

		# check whether network and layer_elevations have the same number of rows and columns
		if self.flopymodel is not None:
			if self.network.shape[0] != self.nrows:
				raise ValueError(f"The network has {self.network.shape[0]} rows but the model has {self.nrows} rows!")
			if self.network.shape[1] != self.ncols:
				raise ValueError(f"The network has {self.network.shape[1]} columns but the model has {self.cols} columns!")

		# iterate over the network
		for i in range(self.network.shape[0]):
			# iterate over columns
			for j in range(self.network.shape[1]):
				# check if node is present
				if self.network[i, j] != 1.:
					continue
				else:
					# check that cell is not on boundary
					if (i > 0 and
						j > 0 and
						i < self.network.shape[0] - 1 and
						j < self.network.shape[1] - 1):
						
						# check if node has diagonal neighbor at all:
						if (self.network[i-1, j-1] == 1. or
							self.network[i-1, j+1] == 1. or
							self.network[i+1, j+1] == 1. or
							self.network[i+1, j-1] == 1.):
							
							# check if a node is present at (ROW-1, COL-1) (up-left diagonal neighbor)
							if (self.network[i-1, j-1] == 1.):
								# don't do anything if a node is already present as a upper / sideways
								# neighbor
								if (self.network[i-1, j] == 1. or self.network[i, j-1] == 1.):
									pass
								# randomly add a node either as a upper or sideways neighbor 
								else:
									if (np.random.random() >= 0.5):
										self.network[i-1, j] = 1.
									else:
										self.network[i, j-1] = 1.
							
							# check if a node is present at (ROW-1, COL+1) (up-right diagonal neighbor)
							if (self.network[i-1, j+1] == 1.):
								# don't do anything if a node is already present as a upper / sideways
								# neighbor
								if (self.network[i-1, j] == 1. or self.network[i, j+1] == 1.):
									pass
								# randomly add a node either as a upper or sideways neighbor 
								else:
									if (np.random.random() >= 0.5):
										self.network[i-1, j] = 1.
									else:
										self.network[i, j+1] = 1.
							# check if a node is present at (ROW+1, COL+1) (down-right diagonal neighbor)
							if (self.network[i+1, j+1] == 1.):
								# don't do anything if a node is already present as a upper / sideways
								# neighbor
								if (self.network[i+1, j] == 1. or self.network[i, j+1] == 1.):
									pass
								# randomly add a node either as a upper or sideways neighbor 
								else:
									if (np.random.random() >= 0.5):
										self.network[i+1, j] = 1.
									else:
										self.network[i, j+1] = 1.
							
							# check if a node is present at (ROW+1, COL-1) (down-left diagonal neighbor)	
							if (self.network[i+1, j-1] == 1.):
								# don't do anything if a node is already present as a upper / sideways
								# neighbor
								if (self.network[i+1, j] == 1. or self.network[i, j-1] == 1.):
									pass
								# randomly add a node either as a upper or sideways neighbor 
								else:
									if (np.random.random() >= 0.5):
										self.network[i+1, j] = 1.
									else:
										self.network[i, j-1] = 1.
						
						# if no diagonal neighbor is found: don't do anything and proceed to the next node
						else:
							continue

		# if elevations are given, change node heights to corresponding elevations
		if self.elevations is not None:
			self.network *= self.elevations

		# set validity to true
		self.valid = True
		return self.network
