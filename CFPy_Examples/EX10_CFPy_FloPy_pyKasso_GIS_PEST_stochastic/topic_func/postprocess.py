# postprocess.py v0.2
# Last updated: 2023-12-20

# Import built-in libraries
import os
import re

# Import third party libraries
import flopy
import flopy.export.utils as fpu
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from flopy.utils import ZoneBudget
from matplotlib.collections import LineCollection
from matplotlib.ticker import AutoMinorLocator
from rasterio.transform import from_origin
from scipy.spatial import ConvexHull
from shapely.geometry import (
    LineString,
    MultiLineString,
    MultiPoint,
    MultiPolygon,
    Point,
    Polygon,
)


class plot_model:
    def __init__(self, model, layer, time, hdobj, cbb, settings):
        """Plot the model.

        Parameters
        ----------
        model : flopy.modflow.Modflow
            The FloPy model.
        layer : int
            The layer to plot.
        time : float
            The time to plot. (see: mf.dis.get_totim())
        hdobj : flopy.utils.HeadFile
            The head object.
        cbb : flopy.utils.CellBudgetFile
            The cell budget file.
        settings : dict
            Parameters for the plot and additional settings.
        """

        self.model = model
        self.layer = layer
        self.time = time
        self.hdobj = hdobj
        self.cbb = cbb
        self.settings = settings

    def head_map(self, ax, **kwargs):
        """
        Plot the simulated heads.

        Parameters:
        ----------------
        - model (flopy.modflow.Modflow): The FloPy model.
        - ax (matplotlib.axes._subplots.AxesSubplot): The axes of the plot.
        - layer (int): The layer to plot.
        - time (float): The time to plot.
        - hdobj (flopy.utils.HeadFile): The head object.
        - settings (dict): Parameters for the plot and additional settings.
        - **kwargs: Additional keyword arguments.
            For example:
            - cmap (str): The colormap.
            - cmap_show (bool): Show the colormap.
            - vmin (int): The minimum value of the colorbar.
            - vmax (int): The maximum value of the colorbar.
            - cmap_alpha (float): The transparency of the colormap.
            - cbar_label (str): The label of the colorbar.
            - masked_values (list): The values to mask.
            - cbar_shrink (float): The shrink of the colorbar.
            - cbar_rotation (float): The rotation of the colorbar.
            - cbar_labelpad (float): The labelpad of the colorbar.
            - cbar_show (bool): Show the colorbar.
            - grid_alpha (float): The transparency of the grid.


        Returns:
        ----------------
        None
        """
        # Extract the settings from the settings dictionary
        if self.settings is None:
            self.settings = {}
        else:
            for key, value in self.settings.items():
                if kwargs.get(key) is None:
                    kwargs[key] = value

        # Extract the heads
        head = self.hdobj.get_data(totim=self.time)
        head[head == -999.99] = np.nan

        # Plot the heads
        mapview = flopy.plot.PlotMapView(model=self.model, layer=self.layer, ax=ax)

        quadmesh = mapview.plot_array(
            head,
            ax=ax,
            masked_values=kwargs.get("masked_values")
            if "masked_values" in kwargs
            else [-1.0e20, -2.0e20],
            alpha=kwargs.get("cmap_alpha") if "cmap_alpha" in kwargs else 0.5,
            cmap=kwargs.get("cmap") if "cmap" in kwargs else "viridis_r",
            vmin=kwargs.get("vmin") if "vmin" in kwargs else 0,
            vmax=kwargs.get("vmax") if "vmax" in kwargs else 100,
        )

        # Plot the colorbar
        cbar = plt.colorbar(
            quadmesh,
            shrink=kwargs.get("cbar_shrink") if "cbar_shrink" in kwargs else 0.7,
        )

        cbar.ax.set_ylabel(
            kwargs.get("cbar_label") if "cbar_label" in kwargs else "Head [m]",
            rotation=kwargs.get("cbar_rotation") if "cbar_rotation" in kwargs else 90,
            labelpad=kwargs.get("cbar_labelpad") if "cbar_labelpad" in kwargs else 10,
        )

        if kwargs.get("cbar_show") == False:
            cbar.remove()

        # Plot the grid
        linecollection = mapview.plot_grid(
            linewidths=0.5,
            color="black",
            alpha=kwargs.get("grid_alpha") if "grid_alpha" in kwargs else 0.5,
            ax=ax,
        )

        # Hide the grid lines
        ax.grid(False, which="both")

        # Plot the labels
        ax.set_xlabel("x-coordinate (m)")
        ax.set_ylabel("y-coordinate (m)")

        # Plot the title
        ax.set_title(f"Head map at time {self.time} s and layer {self.layer}")

    def contour_map(self, ax, levels, **kwargs):
        """
        Plot contour lines of the simulated heads / equipotential lines.

        Parameters:
        ----------------
        - model (flopy.modflow.Modflow): The FloPy model.
        - ax (matplotlib.axes._subplots.AxesSubplot): The axes of the plot.
        - layer (int): The layer to plot.
        - time (float): The time to plot.
        - hdobj (flopy.utils.HeadFile): The head object.
        - settings (dict): Parameters for the plot and additional settings.
        - **kwargs: Additional keyword arguments.
            For example:
            - contour_alpha (float): The transparency of the contour.
            - contour_colors (str): The color of the contour ('whitesmoke').
            - contour_linewidths (float): The linewidths of the contour.
            - masked_values (list): The values to mask.
            - contourlabel_format (str): The format of the contour labels.
            - contourlabel_fontsize (float): The fontsize of the contour labels.
            - cmap (str): The colormap.
            - vmin (int): The minimum value of the colorbar.
            - vmax (int): The maximum value of the colorbar.
            - cbar_show (bool): Show the colorbar.



        Returns:
        ----------------
        None
        """

        # Extract the settings from the settings dictionary
        if self.settings is None:
            self.settings = {}
        else:
            for key, value in self.settings.items():
                if kwargs.get(key) is None:
                    kwargs[key] = value

        # Extract the heads
        head = self.hdobj.get_data(totim=self.time)
        head[head == -999.99] = np.nan

        # Plot the contour lines
        mapview = flopy.plot.PlotMapView(model=self.model, layer=self.layer, ax=ax)

        c = mapview.contour_array(
            head,
            masked_values=kwargs.get("masked_values")
            if "masked_values" in kwargs
            else [-1.0e20, -2.0e20],
            linewidths=kwargs.get("contour_linewidths")
            if "contour_linewidths" in kwargs
            else 1.0,
            levels=levels,
            colors=kwargs.get("contour_colors") if "contour_colors" in kwargs else None,
            alpha=kwargs.get("contour_alpha") if "contour_alpha" in kwargs else 0.8,
            ax=ax,
            cmap=kwargs.get("cmap") if "cmap" in kwargs else None,
            vmin=kwargs.get("vmin") if "vmin" in kwargs else None,
            vmax=kwargs.get("vmax") if "vmax" in kwargs else None,
        )

        if not "contour_colors" in kwargs:
            # Plot the colorbar
            cbar = plt.colorbar(
                c, shrink=kwargs.get("cbar_shrink") if "cbar_shrink" in kwargs else 0.7
            )
            if kwargs.get("cbar_show") == False:
                cbar.remove()

        plt.clabel(
            c,
            fmt=kwargs.get("contourlabel_format")
            if "contourlabel_format" in kwargs
            else "%.3fm",
            fontsize=kwargs.get("contourlabel_fontsize")
            if "contourlabel_fontsize" in kwargs
            else 12.0,
        )

    def vector_map(self, ax, vec, **kwargs):
        """Plot the flow direction.

        Parameters
        ----------
        ax : Matplotlib axes object
            The plot axes.
        vec : str
            The type of vector (volumetric or specific).
        """

        # Extract the settings from the settings dictionary
        if self.settings is None:
            self.settings = {}
        else:
            for key, value in self.settings.items():
                if kwargs.get(key) is None:
                    kwargs[key] = value

        # Extract the heads
        head = self.hdobj.get_data(totim=self.time)
        head[head == -999.99] = np.nan

        # Extract the flow direction
        frf = self.cbb.get_data(text="FLOW RIGHT FACE")[0]
        fff = self.cbb.get_data(text="FLOW FRONT FACE")[0]
        if head.shape[0] == 1:
            flf = None
        else:
            flf = self.cbb.get_data(text="FLOW LOWER FACE")[0]

        qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(
            (frf, fff, flf), self.model
        )  # no head array for volumetric discharge

        sqx, sqy, sqz = flopy.utils.postprocessing.get_specific_discharge(
            (frf, fff, flf), self.model, head
        )

        # Plot the flow direction
        mapview = flopy.plot.PlotMapView(model=self.model, layer=self.layer, ax=ax)

        if vec == "volumetric":
            quiver = mapview.plot_vector(
                qx,
                qy,
                istep=kwargs.get("istep") if "istep" in kwargs else 4,
                jstep=kwargs.get("jstep") if "jstep" in kwargs else 4,
                normalize=kwargs.get("normalize") if "normalize" in kwargs else False,
                color=kwargs.get("vector_color")
                if "vector_color" in kwargs
                else "black",
                alpha=kwargs.get("vector_alpha") if "vector_alpha" in kwargs else 0.5,
                scale=kwargs.get("vector_scale") if "vector_scale" in kwargs else 35,
            )

            # ax.set_title("Volumetric discharge vectors (" + r"$L^3/T$" + ")")

        elif vec == "specific":
            quiver = mapview.plot_vector(
                sqx,
                sqy,
                istep=kwargs.get("istep") if "istep" in kwargs else 4,
                jstep=kwargs.get("jstep") if "jstep" in kwargs else 4,
                normalize=kwargs.get("normalize") if "normalize" in kwargs else False,
                color=kwargs.get("vector_color")
                if "vector_color" in kwargs
                else "black",
                alpha=kwargs.get("vector_alpha") if "vector_alpha" in kwargs else 0.5,
                scale=kwargs.get("vector_scale") if "vector_scale" in kwargs else 35,
            )

            # ax.set_title("Specific discharge vectors (" + r"$L/T$" + ")")

    def crosssection_map(self, ax, line, **kwargs):
        """
        Plot the simulated heads.

        Parameters:
        ----------------
        - model (flopy.modflow.Modflow): The FloPy model.
        - ax (matplotlib.axes._subplots.AxesSubplot): The axes of the plot.
        - time (float): The time to plot.
        - line (dict): The line to plot.
        - hdobj (flopy.utils.HeadFile): The head object.
        - cbb (flopy.utils.CellBudgetFile): The cell budget file.
        - settings (dict): Parameters for the plot and additional settings.
        - **kwargs: Additional keyword arguments.
            For example:
            - cmap (str): The colormap.
            - cmap_show (bool): Show the colormap.
            - vmin (int): The minimum value of the colormap.
            - vmax (int): The maximum value of the colormap.
            - cbar_label (str): The label of the colorbar.
            - cmap_alpha (float): The transparency of the colormap.
            - normalize (bool): Normalize the arrows.
            - cbar_rotation (float): The rotation of the colorbar.
            - cbar_labelpad (float): The labelpad of the colorbar.
            - hstep (int): The step for the x-axis.
            - vector_color (str): The color of the arrows.
            - vector_alpha (float): The transparency of the arrows.
            - vector_scale (float): The scale of the arrows.
            - d_vec (str): The type of vector (volumetric or specific).



        Returns:
        ----------------
        None
        """

        # Extract the settings from the settings dictionary
        if self.settings is None:
            self.settings = {}
        else:
            for key, value in self.settings.items():
                if kwargs.get(key) is None:
                    kwargs[key] = value

        # Extract the heads
        head = self.hdobj.get_data(totim=self.time)
        head[head == -999.99] = np.nan

        # Plot the heads
        modelxsect = flopy.plot.PlotCrossSection(model=self.model, line=line, ax=ax)
        # Plot the grid
        linecollection = modelxsect.plot_grid(
            linewidths=0.5, color="black", alpha=0.5, ax=ax
        )

        t = modelxsect.plot_array(
            head,
            masked_values=[-1.0e20, -2.0e20],
            alpha=kwargs.get("cmap_alpha") if "cmap_alpha" in kwargs else 0.5,
            cmap=kwargs.get("cmap") if "cmap" in kwargs else None,
            vmin=kwargs.get("vmin") if "vmin" in kwargs else None,
            vmax=kwargs.get("vmax") if "vmax" in kwargs else None,
            head=head,
            ax=ax,
        )
        cbar = plt.colorbar(t, shrink=1.0, ax=ax)
        cbar.ax.set_ylabel(
            kwargs.get("cbar_label") if "cbar_label" in kwargs else "Head [m]",
            rotation=kwargs.get("cbar_rotation") if "cbar_rotation" in kwargs else 90,
            labelpad=kwargs.get("cbar_labelpad") if "cbar_labelpad" in kwargs else 10,
        )

        if kwargs.get("cmap_show") == False:
            cbar.remove()

        modelxsect.plot_ibound()

        # Extract the flow direction
        frf = self.cbb.get_data(text="FLOW RIGHT FACE")[0]
        fff = self.cbb.get_data(text="FLOW FRONT FACE")[0]
        if head.shape[0] == 1:
            flf = None
        else:
            flf = self.cbb.get_data(text="FLOW LOWER FACE")[0]

        qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(
            (frf, fff, flf), self.model
        )  # no head array for volumetric discharge

        sqx, sqy, sqz = flopy.utils.postprocessing.get_specific_discharge(
            (frf, fff, flf), self.model, head
        )

        dvec = kwargs.get("vec") if "vec" in kwargs else "specific"

        if dvec == "volumetric":
            quiver = modelxsect.plot_vector(
                qx,
                qy,
                qz,
                head=head,
                hstep=kwargs.get("hstep") if "hstep" in kwargs else 3,
                normalize=kwargs.get("normalize") if "normalize" in kwargs else False,
                color=kwargs.get("vector_color")
                if "vector_color" in kwargs
                else "black",
                alpha=kwargs.get("vector_alpha") if "vector_alpha" in kwargs else 0.5,
                scale=kwargs.get("vector_scale") if "vector_scale" in kwargs else 50,
                headwidth=2.0,
                headlength=1.5,
                headaxislength=1.5,
                zorder=10,
                ax=ax,
            )

        elif dvec == "specific":
            quiver = modelxsect.plot_vector(
                sqx,
                sqy,
                sqz,
                head=head,
                hstep=kwargs.get("hstep") if "hstep" in kwargs else 3,
                normalize=kwargs.get("normalize") if "normalize" in kwargs else False,
                color=kwargs.get("vector_color")
                if "vector_color" in kwargs
                else "black",
                alpha=kwargs.get("vector_alpha") if "vector_alpha" in kwargs else 0.5,
                scale=kwargs.get("vector_scale") if "vector_scale" in kwargs else 50,
                headwidth=2.0,
                headlength=1.5,
                headaxislength=1.5,
                zorder=10,
                ax=ax,
            )

        # Hide the grid lines
        ax.grid(False, which="both")
        # Plot the labels
        ax.set_xlabel("Distance (m)")
        ax.set_ylabel("z-coordinate (m)")
        # Plot the title
        # ax.set_title(f"Head map at time {self.time} s")

    def hdobj_to_rst(self, rst_path):
        """Export the heads as a raster.

        Parameters
        ----------
        rst_path : str
            Path to save the raster.

        """

        head = self.hdobj.get_data(totim=self.time)
        head[head == -999.99] = np.nan

        # Get the model grid information
        xll, yll = (
            self.model.modelgrid.xvertices[0, 0],
            self.model.modelgrid.yvertices[0, 0],
        )
        cellsize = self.model.modelgrid.delr[0]

        # Calculate the geotransform for the raster
        transform = from_origin(xll, yll, cellsize, cellsize)

        # Save the head data as a raster
        with rasterio.open(
            rst_path,
            "w",
            driver="GTiff",
            height=head.shape[1],
            width=head.shape[2],
            count=head.shape[0],
            dtype=head.dtype,
            crs=self.model.modelgrid.proj4,
            transform=transform,
        ) as dst:
            dst.write(head, range(1, head.shape[0] + 1))

        print("Raster file was saved at:", rst_path)

    def contour_to_shp(self, shp_path, levels):
        """Export the contour lines as a shapefile.

        Parameters
        ----------
        shp_path : str
            Path to save the shapefile.
        levels : list
            List with the contour levels.
        """

        head = self.hdobj.get_data(totim=self.time)
        head[head == -999.99] = np.nan

        # Export contours to shapefile
        fpu.export_array_contours(
            filename=shp_path,
            a=head[self.layer, :, :],
            modelgrid=self.model.modelgrid,
            levels=levels,
            crs=self.model.modelgrid.proj4,
        )


def bar_plot_budget(data, ax, stper, settings, **kwargs):
    """
    Plot the budget of the model.

    Parameters:
    - data (pd.DataFrame): The budget data.
    - ax (matplotlib.axes.Axes): The axes to plot on.
    - stper (int): The stress periods as indexes. (e.g. 0 or mf.dis.nper-1;
    - settings (dict): Parameters for the plot and additional settings.
    - **kwargs: Additional keyword arguments.
        For example:
        - ylabel (str): The label of the y-axis.
        - y_pos (float): The position of the text on the y-axis (between -1 and 1). Default is 1.

    Returns:
    Bar plot of the budget.
    """

    # Extract the settings from the settings dictionary
    if settings is None:
        settings = {}
    else:
        for key, value in settings.items():
            if kwargs.get(key) is None:
                kwargs[key] = value

    color = ["lightblue" if i >= 0 else "tomato" for i in data.iloc[stper, :]]

    data.iloc[stper, :].plot(kind="bar", color=color, alpha=0.75, ax=ax)

    ylim = ax.get_ylim()

    for x, y in enumerate(data.iloc[stper, :]):
        ax.text(
            x,
            kwargs.get("y_pos") * max([abs(x) for x in ylim]) * 0.9
            if "y_pos" in kwargs
            else y * 0.9,
            kwargs.get("text_format").format(y)
            if "text_format" in kwargs
            else "{:,.0f}".format(y),
            ha="center",
            rotation=0,
            va="center",
            bbox=dict(facecolor="white", alpha=0.0),
        )

    # Show minor ticks for y-axis
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    ax.grid(color=(0.5, 0.5, 0.5, 0.15), linestyle="-", linewidth=0.2, which="minor")
    ax.grid(color=(0.5, 0.5, 0.5, 0.3), linestyle="-", linewidth=0.5, which="major")

    ax.axhline(0, color="black", alpha=0.3)
    ax.set_ylabel(kwargs.get("ylabel") if "ylabel" in kwargs else "")
    ax.set_ylim(-max([abs(x) for x in ylim]), max([abs(x) for x in ylim]))


def PointsOnCircle(r=0.5, N=24):  # r should be smaller or equal then half of cell size
    """Create points on a circle.

    Args:
        r (float, optional): Relative radius of circle. Defaults to 0.5.
        N (int, optional): Number of points on circle. Defaults to 24.

    Returns:
        np.ndarray: Array with points on circle.
    """

    rv = [2 * np.pi / N * n_i for n_i in range(0, N)]
    return np.vstack((r * np.cos(rv), r * np.sin(rv))).T


def partlocs_structured_grid(model, stpd, poc, d8, steps):
    """Create particle locations for structured grid.

    Args:
        model (flopy.modflow.Modflow): The FloPy model.
        stpd (rec.array): Stress period data.
        poc (np.array): Array with points on circle.
        d8 (bool): If True, create 8 points also around the selected grid cell.
        steps (int): Steps between cells. (If all cells are selected, steps = 1)

    Returns:
        localx: List with x-coordinates of particles.
        localy: List with y-coordinates of particles.
        localz: List with z-coordinates of particles.
        partlocs: List with particle locations.
    """

    localx = []
    localy = []
    localz = []
    partlocs = []

    for i in range(0, stpd.shape[0] - 1, steps):
        for lay, row, col, *_ in [stpd[i]]:
            # Your code here
            r = row
            c = col
            for pr_i in poc:
                lx = (model.dis.delc.array[0] / 2 + pr_i[0]) / model.dis.delc.array[0]
                ly = (model.dis.delr.array[0] / 2 + pr_i[0]) / model.dis.delr.array[0]
                lz = 0.5  # This should be between 0 or 1
                partlocs.append((lay, r, c))  # layer (0 --> top layer), row, column
                if d8 == True:
                    partlocs.append((lay, r - 1, c))
                    partlocs.append((lay, r + 1, c))
                    partlocs.append((lay, r, c + 1))
                    partlocs.append((lay, r, c - 1))
                    partlocs.append((lay, r - 1, c - 1))
                    partlocs.append((lay, r + 1, c + 1))
                    partlocs.append((lay, r - 1, c + 1))
                    partlocs.append((lay, r + 1, c - 1))

                localx.append(lx)
                localy.append(ly)
                localz.append(lz)

    if d8 == True:
        localx = localx * 9
        localy = localy * 9
        localz = localz * 9

    return localx, localy, localz, partlocs


def plot_pathlines_col(mppth_path, cmap, vmin, vmax, ax, label, lw):
    """Plot pathlines with cmap.

    Args:
        mppth_path (str): Path to the pathline file.
        cmap (str): Colormap.
        vmin (float): Minimum value of the colormap. (e.g. 0 seconds)
        vmax (float): Maximum value of the colormap. (e.g. 86400 seconds)
        ax (matplotlib.axes.Axes): The axes to plot on.
        label (str): Colorbar label.
        lw (float): Linewidth of the pathlines.
    """

    # Pathlines

    fname = mppth_path
    p = flopy.utils.PathlineFile(fname)

    # Define colormap and normalization
    colormap = plt.cm.get_cmap(cmap)
    normalize = plt.Normalize(vmin / 86400, vmax / 86400)

    lines = []

    for i in range(0, p.get_maxid() + 1):
        pi = p.get_data(partid=i)
        x = pi["x"]
        y = pi["y"]
        time = pi["time"]

        points = np.column_stack((x, y))
        segments = np.concatenate([points[:-1, None], points[1:, None]], axis=1)
        time_normalized = normalize(time[:-1] / 86400)
        colors = colormap(time_normalized)

        # Create LineCollection segments with corresponding colors
        lines.append(LineCollection(segments, colors=colors, lw=lw, alpha=0.6))

    for line in lines:
        ax.add_collection(line)

    # Adding a colorbar
    sm = plt.cm.ScalarMappable(cmap=colormap, norm=normalize)
    sm.set_array([])

    # Add colorbar
    cbar = plt.colorbar(sm, label=label, ax=ax, extend="max")


def save_pathlines_shp(filepath, mppth_path):
    """Save pathlines as shapefile.

    Args:
        filepath (str): Path to save the shapefile.
        mppth_path (str): Path to the pathline file.
    """

    fname = mppth_path
    p = flopy.utils.PathlineFile(fname)

    lines = []

    for i in range(0, p.get_maxid() + 1, 1):
        pi = p.get_data(partid=i)
        x = pi["x"]
        y = pi["y"]
        time = pi["time"]

        points = np.column_stack((x, y))
        segments = np.concatenate([points[:-1, None], points[1:, None]], axis=1)
        if i == 0:
            times = time[1:]
        else:
            times = np.concatenate([times, time[1:]], axis=0)

        lines.append(segments)

    # Convert LineCollection segments to LineString geometries
    geometries = [LineString(segment) for line in lines for segment in line]

    # Create a GeoDataFrame
    gdf = gpd.GeoDataFrame(geometry=geometries)
    gdf["id"] = gdf.index
    gdf["time_s"] = times
    gdf["time_d"] = times / 86400
    gdf["time_y"] = times / (86400 * 365.25)

    # Save as a shapefile
    gdf.to_file(filepath)


def pathline_catchment(mppth_path, shapefile_path):
    """Create a convex hull around the pathlines. This can be used to create a catchment area.

    Args:
        mppth_path (str): Path to the pathline file.
        shapefile_path (str): Path to save the shapefile. If None the shapefile will not be saved.

    returns:
        hull_gdf (gpd.GeoDataFrame): GeoDataFrame with the convex hull.
    """

    fname = mppth_path
    p = flopy.utils.PathlineFile(fname)

    for i in range(0, p.get_maxid() + 1, 1):
        pi = p.get_data(partid=i)
        x = pi["x"]
        y = pi["y"]
        time = pi["time"]

        points = np.column_stack((x, y)).tolist()
        if i == 0:
            all_points = points
        else:
            all_points += points

    points = np.array(all_points)

    # Compute the convex hull
    hull = ConvexHull(points, incremental=True)

    # Plotting the convex hull
    for simplex in hull.simplices:
        plt.plot(points[simplex, 0], points[simplex, 1], "k-", linewidth=1)

    # Compute the convex hull using scipy
    hull = ConvexHull(points)

    # Extract vertices of the convex hull
    hull_points = points[hull.vertices]

    # Create a MultiPoint object from the hull points
    multi_point = MultiPoint(hull_points)

    # Compute the convex hull within Shapely
    convex_hull = multi_point.convex_hull

    # Apply buffer(0) to ensure there are no intersections
    clean_convex_hull = convex_hull.buffer(0)

    # Convert the convex hull to a GeoDataFrame
    hull_gdf = gpd.GeoDataFrame(geometry=[clean_convex_hull])

    # Save the GeoDataFrame as a shapefile
    if shapefile_path != None:
        hull_gdf.to_file(shapefile_path)

    return hull_gdf


def rename_zonbud_col_df(df, aliases, **kwargs):
    """Rename the columns of the ZoneBudget DataFrame
    (e.g. FROM_ZONE_0 --> Inactive)

    Parameters
    ----------
    df : pd.DataFrame
        ZoneBudget DataFrame
    aliases : dict
        Dictionary with the aliases of the zone numbers
    **kwargs : dict
        Additional keyword arguments to rename the column.
        For example:
        new = {"FROM_ZONE_0": "Inactive"}

    Returns
    -------
    pd.DataFrame
        ZoneBudget DataFrame with renamed columns
    names : dict
        Dictionary with the previous and new names
    """

    names = {}

    for key, value in aliases.items():
        names[f"FROM_{value}"] = f"To/From {value}"
        names[f"TO_{value}"] = f"To/From {value}"

    names["ZONE_0"] = "To/From inactive"
    names["FROM_ZONE_0"] = "To/From inactive"
    names["TO_ZONE_0"] = "To/From inactive"
    names["FROM_CONSTANT_HEAD"] = "Isohypse"
    names["TO_CONSTANT_HEAD"] = "Isohypse"
    names["FROM_STORAGE"] = "Storage"
    names["TO_STORAGE"] = "Storage"

    if kwargs.get(key) is None:
        dict_ = kwargs
        for key, value in dict_.items():
            for key_2, value_2 in dict_[key].items():
                names[key_2] = value_2

    df.rename(columns=names, inplace=True)
    return df, names


def bar_plot_zb(df_to, df_from, ax, cmap, stacked, leg_col):
    """Create a bar plot of the zone budget.

    Parameters
    ----------
    df_to : pd.DataFrame
            DataFrame containing the zone budget to values.
    df_from : pd.DataFrame
            DataFrame containing the zone budget from values.
    ax : matplotlib.pyplot.axes
            Axes to plot the bar plot on.
    cmap : str
            Colormap to use for the bar plot.
    stacked : bool
            Whether to stack the bars or not.
    leg_col : int
            Number of columns in the legend.
    """

    if df_to is not None:
        df_to.drop(["totim"], axis=0).plot(
            kind="bar",
            stacked=stacked,
            cmap=cmap,
            ax=ax,
            linewidth=0.5,
            edgecolor="black",
        )

    if df_from is not None:
        df_from.drop(["totim"], axis=0).plot(
            kind="bar",
            stacked=stacked,
            cmap=cmap,
            ax=ax,
            linewidth=0.5,
            edgecolor="black",
            alpha=0.75,
        )

    plt.xticks(rotation=0, ha="center", va="top")  # ha... horizontal align*

    plt.legend(
        loc="upper center",
        frameon=True,
        bbox_to_anchor=(0.5, -0.05),
        ncol=leg_col,
        fontsize=12,
    )

    # Show minor ticks for y-axis
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    ax.grid(
        color=(0.5, 0.5, 0.5, 0.15),
        linestyle="-",
        linewidth=0.2,
        alpha=0.3,
        which="minor",
    )
    ax.grid(
        color=(0.5, 0.5, 0.5, 0.3),
        linestyle="-",
        linewidth=0.5,
        alpha=0.5,
        which="major",
    )

    plt.axhline(0, color="black", alpha=0.3, label=None)
    plt.ylabel("Flux m$^3$/s \n negative = outflow (To), positive = inflow (From)")
    plt.title("Zone budget")


def read_gwb(modelname, budgetkey, timeunit, path, foldernames, start_datetime):
    """Read the listing file of the model and return the dataframes of the global water budget within an dictionary

    Parameters
    ----------
    modelname : flopy.modflow.Modflow
        The FloPy model.
    budgetkey : str
        Key of the global water budget (e.g. "VOLUMETRIC BUDGET FOR ENTIRE MODEL")
    timeunit : str
        Unit of the time
    path : str
        Path to the folder that contains the model results
    foldernames : list
        List of foldernames that contain the model results (listing file)
    start_datetime : datetime or str
        Datetime of the start of the simulation

    Returns
    -------
    mfl: flopy.utils.MfListBudget
        MODFLOW list budget object
    gwb_dict: dict
        Dictionary with the dataframes of the global water budget
    """

    gwb_dict = {}

    for foldername in foldernames:
        os.chdir(os.path.join(path, foldername))

        mfl = flopy.utils.MfListBudget(
            "{}.list".format(modelname),
            budgetkey=budgetkey,
            timeunit=timeunit,
        )

        df_flux, df_vol = mfl.get_dataframes(start_datetime=start_datetime)
        df_flux_new = df_flux.reset_index()
        df_vol_new = df_vol.reset_index()

        groups = df_vol.groupby(lambda x: x.split("_")[-1], axis=1).groups
        df_vol_in = df_vol.loc[:, groups["IN"]]
        df_vol_in.columns = df_vol_in.columns.map(lambda x: x.split("_IN")[0])

        groups = df_flux.groupby(lambda x: x.split("_")[-1], axis=1).groups
        df_flux_in = df_flux.loc[:, groups["IN"]]
        df_flux_in.columns = df_flux_in.columns.map(lambda x: x.split("_IN")[0])

        df_vol_out = df_vol.loc[:, groups["OUT"]]
        df_vol_out.columns = df_vol_out.columns.map(lambda x: x.split("_OUT")[0])

        df_flux_out = df_flux.loc[:, groups["OUT"]]
        df_flux_out.columns = df_flux_out.columns.map(lambda x: x.split("_OUT")[0])

        if budgetkey == "BUDGET OF THE PIPE SYSTEM":  # Water budget of pipe system
            df_vol_in.rename(columns={"CONSTANT_HEAD": "CONSTANT_HEAD2"}, inplace=True)
            gwb_dict[f"df_vol_in_{foldername}"] = df_vol_in

            df_flux_in.rename(columns={"CONSTANT_HEAD": "CONSTANT_HEAD2"}, inplace=True)
            gwb_dict[f"df_flux_in_{foldername}"] = df_flux_in

            gwb_dict[f"df_vol_out_{foldername}"] = df_vol_out[
                [
                    "CONSTANT_HEAD2",
                    "MATRIX_EXCHANGE",
                    "PIPE_RECHARGE",
                    "CAD_STORAGE",
                    "PFP_STORAGE",
                    "CUMULATIVE",
                    "PERCENT_ERROR",
                ]
            ]

            gwb_dict[f"df_flux_out_{foldername}"] = df_flux_out[
                [
                    "CONSTANT_HEAD2",
                    "MATRIX_EXCHANGE",
                    "PIPE_RECHARGE",
                    "CAD_STORAGE",
                    "PFP_STORAGE",
                    "CUMULATIVE",
                    "PERCENT_ERROR",
                ]
            ]

            gwb_dict[f"df_vol_per_{foldername}"] = pd.DataFrame(
                df_vol_out["PERCENT_ERROR"]
            )

            gwb_dict[f"df_flux_per_{foldername}"] = pd.DataFrame(
                df_flux_out["PERCENT_ERROR"]
            )

            try:
                gwb_dict[f"df_vol_out_{foldername}"].drop(
                    ["PERCENT_ERROR"], axis=1, inplace=True
                )

                gwb_dict[f"df_flux_out_{foldername}"].drop(
                    ["PERCENT_ERROR"], axis=1, inplace=True
                )

            except:
                pass

        else:
            gwb_dict[f"df_vol_in_{foldername}"] = df_vol_in

            gwb_dict[f"df_flux_in_{foldername}"] = df_flux_in

            gwb_dict[f"df_vol_out_{foldername}"] = df_vol_out

            gwb_dict[f"df_flux_out_{foldername}"] = df_flux_out

            df_vol_per = df_vol.loc[:, groups["DISCREPANCY"]]
            df_vol_per.columns = df_vol_per.columns.map(lambda x: x.split("_")[0])
            gwb_dict[f"df_vol_per_{foldername}"] = df_vol_per

            df_flux_per = df_flux.loc[:, groups["DISCREPANCY"]]
            df_flux_per.columns = df_flux_per.columns.map(lambda x: x.split("_")[0])
            gwb_dict[f"df_flux_per_{foldername}"] = df_flux_per

        # Needs always a column with _IN and _OUT
        df_vol_delta = (
            gwb_dict[f"df_vol_in_{foldername}"] - gwb_dict[f"df_vol_out_{foldername}"]
        )
        gwb_dict[f"df_vol_delta_{foldername}"] = df_vol_delta

        df_flux_delta = (
            gwb_dict[f"df_flux_in_{foldername}"] - gwb_dict[f"df_flux_out_{foldername}"]
        )
        gwb_dict[f"df_flux_delta_{foldername}"] = df_flux_delta

    return mfl, gwb_dict


def multiple_gwb(model, gwb_dict, foldernames):
    """This function creates a DataFrame of the global waterbudgets as combined results of multiple models

    Parameters
    ----------
    gwb_dict : dict
        Dictionary with the global waterbudgets of the models
    model : flopy.modflow.Modflow
        The FloPy model. This is needed to get the stress periods
    foldernames : list
        List of foldernames that contain the model results (listing file)

    Returns
    -------
    mgwb_dict : dict
        Dictionary with the dataframes of the global waterbudgets of the models
        See the following keys:

    df_vol_in : DataFrame
        DataFrame with volumetric inflow global waterbudgets of the models
    df_vol_out : DataFrame
        DataFrame with volumetric outflow global waterbudgets of the models
    df_flux_in : DataFrame
        DataFrame with specific flux (inflow) global waterbudgets of the models
    df_flux_out : DataFrame
        DataFrame with specific flux (outflow) global waterbudgets of the models
    df_vol_per : DataFrame
        DataFrame with volumetric percentage discrepancy of the global waterbudgets of the models
    df_flux_per : DataFrame
        DataFrame with specific flux percentage discrepancy of the global waterbudgets of the models
    df_vol_delta : DataFrame
        DataFrame with volumetric delta of the global waterbudgets of the models
    df_flux_delta : DataFrame
        DataFrame with specific flux delta of the global waterbudgets of the models
    """

    # Create an empty DataFrame
    df = pd.DataFrame()
    count = 0

    for t in range(model.dis.nper):
        for foldername in foldernames:
            count += 1
            # os.chdir(os.path.join(Model_sim_path, foldername))
            data_vol_in = gwb_dict[f"df_vol_in_{foldername}"].iloc[t, :]
            data_vol_in["Model"] = foldername
            data_vol_in["Stress period"] = t + 1
            data_vol_in["Time"] = model.dis.get_totim()[t]
            data_vol_in = pd.DataFrame(data_vol_in).T
            if count == 1:
                df_vol_in = pd.concat([df, data_vol_in], ignore_index=True, axis=0)
            else:
                df_vol_in = pd.concat(
                    [df_vol_in, data_vol_in], ignore_index=True, axis=0
                )

            data_vol_out = gwb_dict[f"df_vol_out_{foldername}"].iloc[t, :]
            data_vol_out["Model"] = foldername
            data_vol_out["Stress period"] = t + 1
            data_vol_out["Time"] = model.dis.get_totim()[t]
            data_vol_out = pd.DataFrame(data_vol_out).T
            if count == 1:
                df_vol_out = pd.concat([df, data_vol_out], ignore_index=True, axis=0)
            else:
                df_vol_out = pd.concat(
                    [df_vol_out, data_vol_out], ignore_index=True, axis=0
                )

            data_flux_in = gwb_dict[f"df_flux_in_{foldername}"].iloc[t, :]
            data_flux_in["Model"] = foldername
            data_flux_in["Stress period"] = t + 1
            data_flux_in["Time"] = model.dis.get_totim()[t]
            data_flux_in = pd.DataFrame(data_flux_in).T
            if count == 1:
                df_flux_in = pd.concat([df, data_flux_in], ignore_index=True, axis=0)
            else:
                df_flux_in = pd.concat(
                    [df_flux_in, data_flux_in], ignore_index=True, axis=0
                )

            data_flux_out = gwb_dict[f"df_flux_out_{foldername}"].iloc[t, :]
            data_flux_out["Model"] = foldername
            data_flux_out["Stress period"] = t + 1
            data_flux_out["Time"] = model.dis.get_totim()[t]
            data_flux_out = pd.DataFrame(data_flux_out).T
            if count == 1:
                df_flux_out = pd.concat([df, data_flux_out], ignore_index=True, axis=0)
            else:
                df_flux_out = pd.concat(
                    [df_flux_out, data_flux_out], ignore_index=True, axis=0
                )

            data_vol_per = gwb_dict[f"df_vol_per_{foldername}"].iloc[t, :]
            data_vol_per["Model"] = foldername
            data_vol_per["Stress period"] = t + 1
            data_vol_per["Time"] = model.dis.get_totim()[t]
            data_vol_per = pd.DataFrame(data_vol_per).T
            if count == 1:
                df_vol_per = pd.concat([df, data_vol_per], ignore_index=True, axis=0)
            else:
                df_vol_per = pd.concat(
                    [df_vol_per, data_vol_per], ignore_index=True, axis=0
                )

            data_flux_per = gwb_dict[f"df_flux_per_{foldername}"].iloc[t, :]
            data_flux_per["Model"] = foldername
            data_flux_per["Stress period"] = t + 1
            data_flux_per["Time"] = model.dis.get_totim()[t]
            data_flux_per = pd.DataFrame(data_flux_per).T
            if count == 1:
                df_flux_per = pd.concat([df, data_flux_per], ignore_index=True, axis=0)
            else:
                df_flux_per = pd.concat(
                    [df_flux_per, data_flux_per], ignore_index=True, axis=0
                )

            data_vol_delta = gwb_dict[f"df_vol_delta_{foldername}"].iloc[t, :]
            data_vol_delta["Model"] = foldername
            data_vol_delta["Stress period"] = t + 1
            data_vol_delta["Time"] = model.dis.get_totim()[t]
            data_vol_delta = pd.DataFrame(data_vol_delta).T
            if count == 1:
                df_vol_delta = pd.concat(
                    [df, data_vol_delta], ignore_index=True, axis=0
                )
            else:
                df_vol_delta = pd.concat(
                    [df_vol_delta, data_vol_delta], ignore_index=True, axis=0
                )

            data_flux_delta = gwb_dict[f"df_flux_delta_{foldername}"].iloc[t, :]
            data_flux_delta["Model"] = foldername
            data_flux_delta["Stress period"] = t + 1
            data_flux_delta["Time"] = model.dis.get_totim()[t]
            data_flux_delta = pd.DataFrame(data_flux_delta).T
            if count == 1:
                df_flux_delta = pd.concat(
                    [df, data_flux_delta], ignore_index=True, axis=0
                )
            else:
                df_flux_delta = pd.concat(
                    [df_flux_delta, data_flux_delta], ignore_index=True, axis=0
                )

            mgwb_dict = {
                "df_vol_in": df_vol_in,
                "df_vol_out": df_vol_out,
                "df_flux_in": df_flux_in,
                "df_flux_out": df_flux_out,
                "df_vol_per": df_vol_per,
                "df_flux_per": df_flux_per,
                "df_vol_delta": df_vol_delta,
                "df_flux_delta": df_flux_delta,
            }

    return mgwb_dict


def rename_gwb_col_df(df, **kwargs):
    """Rename the columns of the global water budget DataFrame
    (e.g. FROM_ZONE_0 --> Inactive)

    Parameters
    ----------
    df : pd.DataFrame
        Global water budget DataFrame
    **kwargs : dict
        Additional keyword arguments to rename the column.
        For example:
        new = {'RIVER_LEAKAGE': 'River'}

    Returns
    -------
    pd.DataFrame
        ZoneBudget DataFrame with renamed columns
    names : dict
        Dictionary with the previous and new names
    """

    names = {
        "STORAGE": "Storage",
        "CONSTANT_HEAD": "Constant head\n(matrix)",
        "RECHARGE": "Recharge",
        "RIVER_LEAKAGE": "River leakage",
        "WELLS": "Well",
        "ET": "Evapo-\ntranspiration",
        "PIPES": "Pipe\nexchange",
        "CAD_STORAGE": "CAD\nstorage",
        "CONSTANT_HEAD2": "Constant head\n(conduit)",
        "MATRIX_EXCHANGE": "Matrix\nexchange",
        "PFP_STORAGE": "PFP\nstorage",
        "PIPE_RECHARGE": "Pipe\nrecharge",
        "PERCENT": "Matrix\npercent error",
        "PERCENT_ERROR": "Conduit\npercent error",
    }

    try:
        for i in kwargs:
            names.update(kwargs[i])
    except:
        pass

    df.rename(columns=names, inplace=True)

    return df, names


def zb_data_spec(data, zone):
    """Generate a dataframe with the data of a specific zone

    Parameters
    ----------
    data : pd.DataFrame
        Zonebudget dataframe
    zone : str
        Name of the zone

    Returns
    -------
    zone_df : pd.DataFrame
        Dataframe with the data of the specific zone
    """

    df = data.loc[zone]

    try:
        df.drop(["totim"], axis=0, inplace=True)
    except:
        pass

    # Drop rows that are zero

    df = pd.DataFrame(df)
    df = df.loc[(df != 0).any(axis=1)]

    zone_df = df[zone]

    return zone_df


def create_hob_out_file(path_folder, filename):
    """Creates a .hob_out file from a .list file similr to ModelMuse

    Parameters
    ----------
    path_folder : str
        path to the folder that contains the .list file
    filename : str
        name of the .list file

    Returns:
        filename.hob_out : .hob_out file in the same folder as the .list file
    """

    path_list = path_folder
    path_file = os.path.join(path_list, "{}.list".format(filename))

    # Search for a specific line in the file
    search_text = "   NAME              VALUE              VALUE             DIFFERENCE"
    with open(path_file, "r") as f:
        # Read the file line by line
        for line in f:
            # Check if the line contains the search text
            if search_text in line:
                # If the search line is found, read the next lines into a dataframe
                data = []
                for next_line in f:
                    next_line = next_line.strip()
                    # Check if the next line has values separated by commas
                    if "--------------------------" in next_line:
                        print("Skipped")
                    elif " " in next_line:
                        columns = re.split("\s{2,}", next_line.strip())
                        data.append(columns)
                    # If the next line does not have values, stop reading
                    else:
                        break
                # Create a dataframe from the data and print it
                column_names = re.split("\s{1,}", line.strip())
                df = pd.DataFrame(data)  # , columns=column_names)
                # print(df)
                break  # Exit the loop after finding the search line

    df = df.rename(
        columns={0: "NAME", 1: "OBSERVED VALUE", 2: "SIMULATED VALUE", 3: "DIFFERENCE"}
    )

    hob_out = df[["SIMULATED VALUE", "OBSERVED VALUE", "NAME"]]
    hob_out = hob_out.rename(
        columns={
            "NAME": "OBSERVATION NAME",
            "OBSERVED VALUE": "OBSERVED VALUE",
            "SIMULATED VALUE": "SIMULATED EQUIVALENT",
        }
    )

    # write the dataframe to a text file
    path_hob_out = os.path.join(path_list, "{}.hob_out".format(filename))
    hob_out.to_csv(path_hob_out, sep="\t", index=False)

    hob_out["SIMULATED EQUIVALENT"] = hob_out["SIMULATED EQUIVALENT"].astype(float)
    hob_out["OBSERVED VALUE"] = hob_out["OBSERVED VALUE"].astype(float)

    return hob_out


def statistics_mgwb(mgwb_dict):
    """Create a dictionary with the statistics of the global water budget of the (multiple) models

    Parameters
    ----------
    mgwb_dict : dict
        Dictionary with the dataframes of the global water budget of the (multiple) models

    Returns
    -------
    data_dict_temp : dict
        Dictionary with the statistics of the global water budget of the (multiple) models
    """

    data_dict_temp = {}

    for key, df in mgwb_dict.items():
        # Remove not necessary columns of that new dataframe
        # Create a new dataframe with the mean of all values based on the column 'Time' (simulation time steps) for all models
        ### MEAN
        data_dict_temp[key + "_mean_t"] = df.groupby("Time").mean()
        try:
            data_dict_temp[key + "_mean_t"].drop("Stress period", axis=1, inplace=True)
        except:
            pass
        # Create a new dataframe with the mean of all values based on the column 'Stress period' (stress period index) for all models
        data_dict_temp[key + "_mean_st"] = df.groupby("Stress period").mean()
        try:
            data_dict_temp[key + "_mean_st"].drop("Time", axis=1, inplace=True)
        except:
            pass
        # Same procedure for the other statistics
        ### STD (standard deviation)
        data_dict_temp[key + "_std_t"] = df.groupby("Time").std()
        try:
            data_dict_temp[key + "_std_t"].drop("Stress period", axis=1, inplace=True)
        except:
            pass

        data_dict_temp[key + "_std_st"] = df.groupby("Stress period").std()
        try:
            data_dict_temp[key + "_std_st"].drop("Time", axis=1, inplace=True)
        except:
            pass

        ### MIN (minimum)
        data_dict_temp[key + "_min_t"] = df.groupby("Time").min()
        try:
            data_dict_temp[key + "_min_t"].drop("Stress period", axis=1, inplace=True)
        except:
            pass

        data_dict_temp[key + "_min_st"] = df.groupby("Stress period").min()
        try:
            data_dict_temp[key + "_min_st"].drop("Time", axis=1, inplace=True)
        except:
            pass

        ### MAX (maximum)
        data_dict_temp[key + "_max_t"] = df.groupby("Time").max()
        try:
            data_dict_temp[key + "_max_t"].drop("Stress period", axis=1, inplace=True)
        except:
            pass

        data_dict_temp[key + "_max_st"] = df.groupby("Stress period").max()
        try:
            data_dict_temp[key + "_max_st"].drop("Time", axis=1, inplace=True)
        except:
            pass

        ### MEDIAN
        data_dict_temp[key + "_median_t"] = df.groupby("Time").median()
        try:
            data_dict_temp[key + "_median_t"].drop(
                "Stress period", axis=1, inplace=True
            )
        except:
            pass

        data_dict_temp[key + "_median_st"] = df.groupby("Stress period").median()
        try:
            data_dict_temp[key + "_median_st"].drop("Time", axis=1, inplace=True)
        except:
            pass

        # ### 75th percentile
        # data_dict_temp[key + "_75th_t"] = df.groupby("Time").quantile(0.75)
        # try:
        #     data_dict_temp[key + "_75th_t"].drop("Stress period", axis=1, inplace=True)
        # except:
        #     pass

        # data_dict_temp[key + "_75th_st"] = df.groupby("Stress period").quantile(0.75)
        # try:
        #     data_dict_temp[key + "_75th_st"].drop("Time", axis=1, inplace=True)
        # except:
        #     pass

        # ### 25th percentile
        # data_dict_temp[key + "_25th_t"] = df.groupby("Time").quantile(0.25)

        # try:
        #     data_dict_temp[key + "_25th_t"].drop("Stress period", axis=1, inplace=True)
        # except:
        #     pass

        # data_dict_temp[key + "_25th_st"] = df.groupby("Stress period").quantile(0.25)
        # try:
        #     data_dict_temp[key + "_25th_st"].drop("Time", axis=1, inplace=True)
        # except:
        #     pass

    return data_dict_temp
