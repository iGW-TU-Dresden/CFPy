# str_grid_processing.py v0.1
# Last updated: 2023-12-21

import flopy
from flopy.utils import GridIntersect

from shapely.geometry import (
    LineString,
    MultiLineString,
    MultiPoint,
    MultiPolygon,
    Point,
    Polygon,
)


def grid_ix(data_dict, ix, s_key, single_features):
    """Retrieve the grid intersection of a shapefile with the model grid.

    Parameters
    ----------
    data_dict : dict
        Dictionary with the data.
    ix : flopy.utils.GridIntersect
        Grid intersection object.
    s_key : str
        The starting key of the data_dict to search for the shapefile. (e.g. "D-Sp")
    single_features : bool
        If True, only  that have single features of the geodataframes will be intersected with the grid.
        If False, all features of the shapefile will be intersected with the grid.

    Returns
    -------
    dict
        Dictionary with the grid intersection of the shapefile.
    """

    data_dict_temp = {}

    for key, shp in data_dict.items():
        if key.startswith(s_key):
            try:
                # print(key, shp.geometry.geom_type.unique())

                if single_features == True:
                    if len(shp.geometry) == 1:
                        print(key, shp.geometry.geom_type.unique())
                        if shp.geometry.geom_type.unique() == ["Point"]:
                            for nr in range(len(shp.geometry)):
                                data_dict_temp[key + "-ix"] = ix.intersect(
                                    Point(shp.geometry[nr])
                                )

                        elif shp.geometry.geom_type.unique() == ["LineString"]:
                            for nr in range(len(shp.geometry)):
                                data_dict_temp[key + "-ix"] = ix.intersect(
                                    LineString(shp.geometry[nr])
                                )

                        elif shp.geometry.geom_type.unique() == ["Polygon"]:
                            for nr in range(len(shp.geometry)):
                                data_dict_temp[key + "-ix"] = ix.intersect(
                                    Polygon(shp.geometry[nr])
                                )
                else:
                    print(key, shp.geometry.geom_type.unique())
                    if shp.geometry.geom_type.unique() == ["Point"]:
                        for nr in range(len(shp.geometry)):
                            data_dict_temp[key + "-ix"] = ix.intersect(
                                Point(shp.geometry[nr])
                            )

                    elif shp.geometry.geom_type.unique() == ["LineString"]:
                        for nr in range(len(shp.geometry)):
                            data_dict_temp[key + "-ix"] = ix.intersect(
                                LineString(shp.geometry[nr])
                            )

                    elif shp.geometry.geom_type.unique() == ["Polygon"]:
                        for nr in range(len(shp.geometry)):
                            data_dict_temp[key + "-ix"] = ix.intersect(
                                Polygon(shp.geometry[nr])
                            )

            except:
                print(
                    f"Error in {key}, probably no shapefile/geometry. Is this maybe a raster file or an alredy existing intersection dataset?"
                )
    return data_dict_temp
