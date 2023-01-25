import numpy as np
import pandas as pd
import os

class Postprocessor():
	
	"""
	A class for post-processing CFPy results.
	While MODFLOW-related standard outputs can be read with flopy
	utilities, CFP-specific output data can be read with this
	post-processing capability.
	"""

	def __init__(self, modelname=None, path=None):
		"""
		modelname : the name of the model, str
		path : the path to the listing file, str
		"""

		self.modelname = modelname
		self.path = path

class FileReader(Postprocessor):
	"""
	A class to read the CFP output files.
	It is searched for LIST-files in the active working directory
	and all sub-directories (if modelname is not None) or for a specific
	file (if path is not None).
	If no LIST-file corresponding to the given modelname or path is
	available in the working directory or its sub-directories, an error is raised.

	Example:

	fr = CFPy.postprocessing.FileReader("MyModel")
	node_df, tube_df = fr.read_output(node_num=1, tube_num=2)

	fr = CFPy.postprocessing.FileReader("MyModel.list")
	node_df, tube_df = fr.read_output(node_num=1, tube_num=2)
	"""

	def __init__(self, modelname=None, path=None):
		"""
		modelname : the name of the model, str
		"""

		# check that either modelname or pat are given
		try:
			if path is None and modelname is None:
				msg = ("Either path or modelname need to be given")
				raise NameError(msg)
			if path is not None and modelname is not None:
				msg = ("Either path or modelname need to be given")
				raise NameError(msg)
		except NameError:
			raise

		if modelname is not None:
			# check that a listing file is available corresponding to the
			#	modelname; if such a file is found, change the directory
			#	to that location
			# find all files in the current working directory and all sub-
			# 	directories
			dirs_ = []
			for root, dirs, files in os.walk(os.getcwd()):
				if modelname + ".list" in files:
					# if a list-file with the modelname is available,
					#	append the corresponding directory
					dirs_.append(root)
			try:
				if len(dirs_) == 0:
					msg = ("There is no listing file available corresponding to "
						   "the modelname {}".format(modelname))
					raise NameError(msg)
			except NameError:
				raise

			# change directory to where the correct listing file is located
			os.chdir(dirs_[0])

		# initialize the parent object
		Postprocessor.__init__(
			self,
			modelname=modelname,
			path=path
			)

		self.tube_num = None
		self.node_num = None
		self.list_lines = None

	def read_output(self, node_num=None, tube_num=None):
		"""
		Read the listing file and return DataFrames of node and tube state
		variables with each row representing a single time step.
		
		:: Parameters ::
		node_num : an integer or a list of integers specifying the
			node numbers for which to obtain results for, integer or
			list-like
		tube_num : an integer or a list of integers specifying the
			tube numbers for which to obtain results for, integer or
			list-like

		:: Returns ::
		node_df : a pandas DataFrame containing all node state variables
			for all time steps for a given node; pd.DataFrame
		tube_df : a pandas DataFrame containing all tube state variables
			for all time steps for a given tube; pd.DataFrame

		NOTE: This works fine for CFPv2 because the column names are always
			the same. If the column names change, change the corresponding
			list here.
		"""	

		self.node_num = node_num
		self.tube_num = tube_num
		node_df = None 
		tube_df = None

		# check if either a node or tube number is given
		try:
			if self.node_num is None and self.tube_num is None:
				msg = ("Neither node nor tube numbers are given!")
				raise ValueError(msg)
		except ValueError:
			raise
		
		# read all lines from the listing file; but only if it has not
		#	been read yet
		if self.list_lines is None:
			lines = []
			if self.modelname is not None:
				with open(self.modelname + ".list", "r") as file:
					for line in file:
						lines.append(line)
				self.list_lines = lines
			if self.path is not None:
				with open(self.path, "r") as file:
					for line in file:
						lines.append(line)
				self.list_lines = lines

		# handle node
		# make list with lines relevant for the current node
		if self.node_num is not None:
			node_lines = []
			for num, line in enumerate(self.list_lines):
				# find the location where node data is given
				if "RESULTS OF FLOW CALCULATION" in line:
					# append the line corresponding to the desired node
					node_lines.append(self.list_lines[num + 2 + self.node_num])

			# write relevant lines to a temporary file
			with open("temp_node.csv", "w") as file:
				for line in node_lines:
					# if it is a fixed head node, remove the keyword
					if "FIX" in line:
						line = line.replace("FIX", "")
					file.write(line)

			# assign column names
			node_header = ["Node#", "Node Head [L]", "Matrix Head [L]", "Exchange [L3 T-1]",
				"CADS Flow [L3 T-1]", "PFPS Flow [L3 T-1]", "Direct Recharge [L3 T-1]",
				"Q Well [L3 T-1]", "FHLQ", "Cauchy", "Cauchy LQ", "QLH", "Q Fix [L3 T-1]"]

			# read the temporary file to a DataFrame
			node_df = pd.read_table(
				"temp_node.csv",
				header=None,
				delim_whitespace=True,
				names=node_header
				)

			# remove temporary file
			if os.path.exists("temp_node.csv"):
				os.remove("temp_node.csv")

		# handle tube
		# make list with lines relevant for the current tube
		if self.tube_num is not None:
			tube_lines = []
			for num, line in enumerate(self.list_lines):
				# find the location where tube data is given; use the following string to find the data in the model listing file           2023 01 20 !TR: Modified the search string and added explanation
				if "TUBE    B    E" in line:
					# append the line corresponding to the desired tube
					tube_lines.append(self.list_lines[num + self.tube_num])

			# write relevant lines to a temporary file
			with open("temp_tube.csv", "w") as file:
				for line in tube_lines:
					file.write(line)

			# assign column names
			tube_header = ["Tube", "Beginning Node#", "Ending Node#", "Flow Type",
				"Q [L3 T-1]", "Diam. [L]", "Len. [L]", "Re [-]", "Residence Time [T]"]

			# read the temporary file to a DataFrame
			tube_df = pd.read_table(
				"temp_tube.csv",
				header=None,
				delim_whitespace=True,
				names=tube_header
				)

			# remove temporary file
			#if os.path.exists("temp_tube.csv"):
			#	os.remove("temp_tube.csv")

		return node_df, tube_df