# Stochastic modeling example <!-- omit in toc -->

- [Objective](#objective)
- [Software and Dependencies](#software-and-dependencies)

# Objective

The primary objective of this study is to model the dynamic behavior of a karst aquifer system, which transfers water to an alluvial aquifer through subsurface pathways, including a subsurface spring and the karst rock matrix. The investigation focuses on the interconnected region between the karst and alluvial systems to assess system behavior and quantify boundary inflow. Synthetic observation wells are strategically placed to generate a time series of virtual potentiometric head, serving as a reference dataset for simulated data (`01_Syn_model.ipynb`).

Stochastic modeling techniques (`02_Stoch_model.ipynb`) are used to replicate this reference system, considering various realizations of the karst network. Models meeting specific confidence criteria are selected for further analysis such as evaluating the impact of a virtual Managed Aquifer Recharge (MAR) application based on the chosen models.

The examples presented here aim to illustrate diverse aspects of geospatial information utilization, such as using shapefiles (no real world system) to delineate spatial features and create the grid/model domain for Flopy and pyKasso's network generation. Post-processing routines are also provided, which include error analysis, global and local water budget analysis, network probability assessment, and using multiple models to investigate various injection scenarios. Furthermore, the example demonstrate the application of PEST for automatic calibration, providing insights into sensitivity and uncertainty analysis.

>:white_check_mark: **Done:** Synthetic model has been completed, and the stochastic model is nearing completion.

>:construction: **Work in progress:** Currently implementing Zonebudget and comparing simulated and virtual observations of the reference system. Additionally, working on PEST calibration and injection scenario analysis.

>:warning: **Warning:** The current scripts are not entirely finished yet.

# Software and Dependencies

For proper execution of the scripts, follow the instruction for CFPy and pyKasso provided [here][CFPy]. Furthermore, supplementary packages are required, which can be installed using `pip`:

`pip install flopy geopandas mpmath matplotlib notebook numpy pandas pyproj pyshp pyyaml openpyxl rasterio scipy scikit-fmm`

[Links]: #
[CFPy]: https://github.com/marcusgenzel/CFPy/tree/dev "Link to CFPy github"


