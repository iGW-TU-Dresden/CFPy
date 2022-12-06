Workflow Description for Coupling Stochastic Network Generation (e.g., pyKasso),
CFPy and FloPy to Create a Script-Based Framework for Spatially Distributed
Karst Simulation with MODFLOW-CFP

(1) Definition of Model Domain
	- define a rectangular (in the 2D case) or cuboid (in the 3D case) model
		domain together with a spatial discretization (i.e., number and width
		of rows, number and width of columns, number of layers together with
		layer elevation information)
	- if a non-rectangular model domain should be used, still start with a
		rectangular domain. later, after all input data are generated, the
		model domain can be altered in FloPy by setting cells inactive

(2) pyKasso Network Generation
	- include the domain definition from step (1) in the settings.yaml-file for
		pyKasso
	- define other stochastic network generation options in the settings.yaml
	- generate one (or multiple) karst networks with pyKasso

	NOTE: right now, pyKasso can only generate 2D karst networks. However, 
		multiple (vertically connected) node planes can be given to CFPy. To
		give the 2D network some 3D structure nonetheless, node elevations
		need to be given to CFPy. Those elevations can be uniform or non-
		uniform. Also, multiple 2D node planes (with each different uniform
		or non-uniform node elevations) can be connected vertically. See the
		Jupyter Notebook "pyKasso_CFPy_coupling.ipynb" for more information.

(3) Coupling pyKasso and CFPy
	NOTE: CFPy generally relies on one single input-file, where all neccessary
		information is summarized (domain discretization, number of node
		planes, MODFLOW layer elevations, node network and elevations). this
		file can either be generated automatically with the generate_nbr me-
		thod of the CFPy.preprocessing module (preferred option) or can also
		be generated manually

	- use the pyKassoValidator object in the CFPy.preprocessing module to
		validate the node network generated with pyKasso in step (2). provide
		node elevation information in this step
	- option 1 (preferred option): generate the CFPy input information (the
		.nbr-file) automatically via the generate_nbr method in the CFPy.pre-
		processing module (providing the pyKasso.SKS catchment as well as
		domain discretization and layer elevation information)
	- option 2: or generate the .nbr-file manually (also see the Jupyter
		Notebook "pyKasso_CFPy_coupling.ipynb" for more information) by expor-
		ting the validated network and copying to the .nbr-file

(4) Set Up the MODFLOW-CFP Model with CFPy and FloPy
	- create the model with FloPy (using a MODFLOW-CFP distribution, prefer-
		rably the cfpv2.exe version distributed by TU Dresden)
	- create all input-files for CFP with CFPy based on the .nbr-file gener-
		ated in step (3)
	- create the remaining input-files for MODFLOW with FloPy
	- simulate the model with FloPy (model.run_model method)

(5) Post-Processing of Results
	- use FloPy post-processing methods to process classical MODFLOW results
		such as matrix head information etc.
	- use the methods in the CFPy.postprocessing module to obtain and pro-
		cess CFP-specific results (node- and tube- related data)

NOTE: all steps are shown in the example notebooks!
