"""
Module for plotting
	(1) the nodes and conduit nbr_data structure (in 3D)
	(2) the results (heads) of a calculation together with the conduit nbr_data structure

The methods work in a Jupyter Notebook environment with %matplotlib inline for static
plotting and with %matplotlib widget for interactive and movable visualization
"""

import copy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

class NetworkCreator():

	"""
	Class to create data structures used to plot the conduit network

	Parameters
	----------
	elevs : list-like of MODFLOW layer elevations with shape (n_lays, n_rows,
		n_cols); list-like of floats
	nbr_data : the full nbr_data array from the nbr.nbr() method
	"""

	def __init__(
		self,
		elevs,
		nbr_data
		):
		# assuming that there are two elements (multidimensional lists)
		# in elevs, i.e., the bottom elevation array and the conduit
		# elevation array from nbr.nbr_read()

		# the same holds for the nbr_data variable, which holds all 8
		# results of nbr.nbr(), i.e., node numbers, plane numbers, 
		# node locations, conduit locations, node neighbors, tube
		# numbers, tube pairs and tube neighbors

		# it is necessary to make deep copies of the original lists
		# because otherwise the nbr_data gets changed and the module
		# does not work, i.e., MODFLOW / CFP raises an error

		self.bot_elev = copy.deepcopy(elevs[0])
		self.cond_elev = copy.deepcopy(elevs[1])

		self.node_numbers = copy.deepcopy(nbr_data[0])
		self.plane_numbers = copy.deepcopy(nbr_data[1])
		self.node_loc = copy.deepcopy(nbr_data[2])
		self.cond_loc = copy.deepcopy(nbr_data[3])
		self.node_nbr = copy.deepcopy(nbr_data[4])
		self.tube_numbers = copy.deepcopy(nbr_data[5])
		self.tube_pairs = copy.deepcopy(nbr_data[6])
		self.tube_nbr = copy.deepcopy(nbr_data[7])

	def create_network(self):
		"""
		Create the network data structure for plotting

		Parameters
		----------
		None

		Returns
		-------
		pipes_coords
		node_x : list-like of node x-coordinates (along a row)
		node_y : list-like of node y-coordinates (along a col)
		node_z : list-like of node z-coordinates
		node_numbers : list-like of node numbers
		n_rows : integer number of rows
		n_cols : integer number of cols
		"""		
		# get and copy data
		bot_elev = self.bot_elev
		cond_elev = self.cond_elev
		node_loc = self.node_loc
		node_numbers = self.node_numbers
		plane_numbers = self.plane_numbers
		tube_pairs = self.tube_pairs

		pipes_coords = []
		node_x = []
		node_y = []
		node_z = []

		geoheight = []

		# delete following statements if number of rows and columns gets included
		# in the nbr.nbr_read() or nbr.nbr() methods
		# get number of rows
		n_rows = np.shape(bot_elev)[1]
		# get number of columns
		n_cols = np.shape(bot_elev)[2]

		for plane in range(int(max(plane_numbers)[0])):
			for row in range(n_rows):
				for col in range(n_cols):
					# separate value if it has a leading "c" (e.g., "c25.1")
					#     meaning the node is vertically connected
					geoheight_val = str(cond_elev[plane][row][col])
					if geoheight_val[0] == "c":
					    geoheight_val = geoheight_val[1:]
					if float(geoheight_val) != -999:
					    geoheight.append(geoheight_val)

		# get coordinates
		# create list of x-values (columns)
		node_loc_height = []

		for loc, height in zip(node_loc, geoheight):
			loc[2] = float(height)
			node_loc_height.append(loc)

		for loc in node_loc:
			# subtract 0.5 because locs are 1-based indexed
			# center coordinates to middle of cell
			node_x.append(loc[0] - 0.5)
			node_y.append(loc[1] - 0.5)
			node_z.append(loc[2])
		    
		# get pipe endpoint coordinates
		pipe_x1, pipe_y1, pipe_z1 = [], [], []
		pipe_x2, pipe_y2, pipe_z2 = [], [], []
		
		for pipe_num, con in enumerate(tube_pairs):
			pipe_x1.append(node_x[con[0] - 1])
			pipe_y1.append(node_y[con[0] - 1])
			pipe_z1.append(node_z[con[0] - 1])

			pipe_x2.append(node_x[con[1] - 1])
			pipe_y2.append(node_y[con[1] - 1])
			pipe_z2.append(node_z[con[1] - 1])
    	
		# create list which contains pipe_numbers lists which in turn contain
		# three lists (x-pair, y-pair and z-pair)
		for i in range(len(tube_pairs)):
			pipes_coords.append([])

		for num, pipe in enumerate(pipes_coords):
			pipe.append([pipe_x1[num], pipe_x2[num]])
			pipe.append([pipe_y1[num], pipe_y2[num]])
			pipe.append([pipe_z1[num], pipe_z2[num]])

		return pipes_coords, node_x, node_y, node_z, node_numbers, n_rows, n_cols

class Network(NetworkCreator):

	"""
	Class stores plotting methods to handle the plottable data structures

	Parameters
	----------
	elevs : list-like of MODFLOW layer elevations with shape (n_lays, n_rows,
		n_cols); list-like of floats
	nbr_data : the full nbr_data array from the nbr.nbr() method
	plot_nums : bool specifying whether to plot node and pipe numbers; default
		is True
	rot_x : float specifying the rotation around the x-axis in degrees; default
		is 40
	rot_z : float specifying the rotation around the z-axis in degrees; default
		is 120
	text_shift : float specifying how much to shift the text-boxes away from the
		nodes / tubes (only has an effect if plot_nums is True); default is .1
	dpi : integer number specifying the plotting resolution; default is 100
	kind : string specifying whether MODFLOW layer elevations should be
		interpolated with triangles ("triangular") or with rectangles
		("rectangular"); default is "rectangular"

	NOTE: at the moment, only the plotting utilities for the "raw" network
	visualization and for the combined visualization of network and simulation
	results are included
	"""

	def __init__(
		self,
		elevs,
		nbr_data,
		plot_nums=True,
		rot_x=40,
		rot_z=120,
		text_shift=0.1,
		dpi=100,
		kind="rectangular"
		):

		self.plot_nums = plot_nums
		self.rot_x = rot_x
		self.rot_z = rot_z
		self.text_shift = text_shift
		self.dpi = dpi
		self.kind = kind

		# initialize plottable network data structures
		# here, the elevs and nbr_data arrays are splitted to store the later relevant
		# information
		NetworkCreator.__init__(self, elevs=elevs, nbr_data=nbr_data)
		network_ = self.create_network()
		self.pipes_coords = network_[0]
		self.node_x = network_[1]
		self.node_y = network_[2]
		self.node_z = network_[3]
		self.node_numbers = network_[4]
		self.n_rows = network_[5]
		self.n_cols = network_[6]

	def plot_network(self, plot_nums=None, rot_x=None, rot_z=None,
		text_shift=None, dpi=None, kind=None, alpha=None, outlet_nodes=None,
		inlet_nodes=None, node_size=50):
		"""
		plotting the raw conduit network without simulation results
		
		Parameters
		----------
		plot_nums : bool specifying whether to plot node and pipe numbers; default
		is True
		rot_x : float specifying the rotation around the x-axis in degrees; default
			is 40
		rot_z : float specifying the rotation around the z-axis in degrees; default
			is 120
		text_shift : float specifying how much to shift the text-boxes away from the
			nodes / tubes (only has an effect if plot_nums is True); default is .1
		dpi : integer number specifying the plotting resolution; default is 100
		kind : string specifying whether MODFLOW layer elevations should be
			interpolated with triangles ("triangular") or with rectangles
			("rectangular"); default is "rectangular"
		alpha : the opacity of the MODFLOW layer elevation surfaces where 0 is
			translucent and 1 is fully opaque; float
		outlet_nodes : an optional list with node numbers (0-based indexing)
			representing outlets / spring nodes of the conduit network and if None,
			all network nodes are drawn equally; list of ints
		inlet_nodes : an optional list with node numbers (0-based indexing)
			representing inlet nodes of the conduit network and if None,
			all network nodes are drawn equally; list of ints
		node_size : a float representing the size of the drawn nodes; float

		Returns
		-------
		ax : the matplotlib axis
		xx_center : x-coordinates of the model meshgrid
		yy_center : y-coordinates of the model meshgrid
		"""

		if plot_nums is None:
			plot_nums = self.plot_nums
		if rot_x is None:
			rot_x = self.rot_x
		if rot_z is None:
			rot_z = self.rot_z
		if text_shift is None:
			text_shift = self.text_shift
		if dpi is None:
			dpi = self.dpi
		if kind is None:
			kind = self.kind
		if alpha is None:
			alpha = 0.1

		pipes_coords = self.pipes_coords
		node_x = self.node_x
		node_y = self.node_y
		node_z = self.node_z
		node_numbers = self.node_numbers
		n_rows = self.n_rows
		n_cols = self.n_cols

		# initialize figure
		self.fig = plt.figure(figsize=(10, 10), dpi=dpi)
		self.ax = self.fig.add_subplot(111, projection="3d")
		cnorm = matplotlib.colors.Normalize(vmin=min(node_z), vmax=max(node_z))
		
		# plot pipe locations
		for num, pipe in enumerate(pipes_coords):
			self.ax.plot(
				pipe[0],
				pipe[1],
				pipe[2],
				c="grey",
				linewidth=5,
				zorder=-1,
				alpha=0.6
				)
			im1 = self.ax.plot(
				pipe[0],
				pipe[1],
				pipe[2],
				c=matplotlib.cm.brg(cnorm((pipe[2][0] + pipe[2][1]) / 2)),
				zorder=-1,
				alpha=0.6
				)
			if plot_nums == True:
				props = dict(boxstyle="round", facecolor="white", alpha=0.8, ls="--")
				self.ax.text(
					(pipe[0][0] + pipe[0][1]) / 2 + text_shift,
					(pipe[1][0] + pipe[1][1]) / 2 + text_shift,
					(pipe[2][0] + pipe[2][1]) / 2 + text_shift,
					str(num + 1),
					c=matplotlib.cm.brg(cnorm((pipe[2][0] + pipe[2][1]) / 2)),
					fontweight="bold",
					bbox=props
					)

		# scatter node locations
		for num, vals in enumerate(zip(node_numbers, node_x, node_y, node_z)):
			if outlet_nodes is not None and num in outlet_nodes:
				marker_ = "D"
				s_ = node_size + 10
			elif inlet_nodes is not None and num in inlet_nodes:
				marker_ = "^"
				s_ = node_size + 10
			else:
				marker_ = "o"
				s_ = node_size
			self.ax.scatter(
				vals[1],
				vals[2],
				vals[3],
				zorder=100,
				color=matplotlib.cm.brg(cnorm(vals[3])),
				s=s_,
				edgecolors="k",
				marker=marker_
				)
			if plot_nums == True:
				props = dict(boxstyle="round", facecolor="white", alpha=0.8)
				self.ax.text(
					vals[1] + text_shift,
					vals[2] + text_shift,
					vals[3] + text_shift,
					vals[0],
					c=matplotlib.cm.brg(cnorm(vals[3])),
					fontweight="bold",
					bbox=props
					)
		
		# create lists with coords of cell centers
		x_center = np.linspace(0.5, n_cols + 0.5, n_cols, False)
		y_center = np.linspace(0.5, n_rows + 0.5, n_rows, False)
		xx_center, yy_center = np.meshgrid(x_center, y_center)
		
		for lay in self.bot_elev:
			lay_elev = np.array(lay, dtype=np.float64)
			if kind == "rectangular":
				self.ax.plot_surface(
					X=xx_center,
					Y=yy_center,
					Z=lay_elev,
					ccount=np.shape(yy_center)[1],
					rcount=np.shape(yy_center)[0],
					shade=True,
					alpha=alpha
					)
				self.ax.plot_wireframe(
					X=xx_center,
					Y=yy_center,
					Z=lay_elev,
					ccount=np.shape(yy_center)[1],
					rcount=np.shape(yy_center)[0],
					alpha=alpha
					)
			elif kind == "triangular":
				lay_elev_flat = lay_elev.flatten()
				X, Y = xx_center.flatten(), yy_center.flatten()
				tri = matplotlib.tri.Triangulation(X, Y)

				# plot triangulated surfaces
				self.ax.plot_trisurf(
					X,
					Y,
					lay_elev_flat,
					triangles=tri.triangles,
					alpha=alpha,
					linewidth=1,
					edgecolor="black"
					)
			else:
				raise ValueError("attribute 'kind' has to be either 'rectangular' or 'triangular'!")
		
		# self.ax.invert_zaxis()
		self.ax.invert_yaxis()
		self.ax.grid(True, which="major")
		if n_cols < 20:
			self.ax.set_xticks([i for i in range(1, n_cols + 1)])
		else:
			self.ax.set_xticks([i for i in range(1, n_cols + 1, 5)])
		if n_rows < 20:
			self.ax.set_yticks([i for i in range(1, n_rows + 1)])
		else:
			self.ax.set_yticks([i for i in range(1, n_rows + 1, 5)])
		# self.ax.set_zticks([i for i in range(1, n_lays + 1)])
		self.ax.set_xlabel("Columns")
		self.ax.set_ylabel("Rows")
		self.ax.set_zlabel("z")
		self.ax.set_title("Node and Conduit Network", loc="center", y=0.9)
		
		# set view
		self.ax.view_init(elev=rot_x, azim=rot_z)
		
		plt.tight_layout()

		return self.ax, xx_center, yy_center

	def plot_results(self, heads, time, layer, n_contours=10, plot_nums=None,
		rot_x=None, rot_z=None, text_shift=None, dpi=None, alpha=None,
		outlet_nodes=None, inlet_nodes=None, node_size=50):
		"""
		plotting the conduit network with simulation results as contour lines
		
		Parameters
		----------
		heads : array of head values with shape (n_timesteps, n_lays, n_rows,
			n_cols); numpy ndarray
		time : time step index for which to plot heads for; int
		layer : the layer index for which to plot heads for; int
		plot_nums : bool specifying whether to plot node and pipe numbers; default
			is True
		rot_x : float specifying the rotation around the x-axis in degrees; default
			is 40
		rot_z : float specifying the rotation around the z-axis in degrees; default
			is 120
		text_shift : float specifying how much to shift the text-boxes away from the
			nodes / tubes (only has an effect if plot_nums is True); default is .1
		dpi : integer number specifying the plotting resolution; default is 100
		kind : string specifying whether MODFLOW layer elevations should be
			interpolated with triangles ("triangular") or with rectangles
			("rectangular"); default is "rectangular"
		alpha : the opacity of the MODFLOW layer elevation surfaces where 0 is
			translucent and 1 is fully opaque; float
		outlet_nodes : an optional list with node numbers (0-based indexing)
			representing outlets / spring nodes of the conduit network and if None,
			all network nodes are drawn equally; list of ints
		inlet_nodes : an optional list with node numbers (0-based indexing)
			representing inlet nodes of the conduit network and if None,
			all network nodes are drawn equally; list of ints
		node_size : a float representing the size of the drawn nodes; float
		"""

		if plot_nums is None:
			plot_nums = self.plot_nums
		if rot_x is None:
			rot_x = self.rot_x
		if rot_z is None:
			rot_z = self.rot_z
		if text_shift is None:
			text_shift = self.text_shift
		if dpi is None:
			dpi = self.dpi
		if alpha is None:
			alpha = 0.1

		self.ax, xx_center, yy_center = self.plot_network(
			plot_nums=plot_nums,
			rot_x=rot_x,
			rot_z=rot_z,
			text_shift=text_shift,
			dpi=dpi,
			alpha=alpha,
			outlet_nodes=outlet_nodes,
			inlet_nodes=inlet_nodes,
			node_size=node_size
			)

		contourdata = np.array(heads[time, layer - 1, :, :])
		ct = self.ax.contour(xx_center, yy_center, contourdata, levels=n_contours,
			alpha=0.5, cmap="rainbow", linewidths=2)

		cb = plt.colorbar(ct, shrink=0.3)
		cb.set_label("Head [m]")
		cb.set_alpha(1)
		cb.draw_all()
		
		return self.ax