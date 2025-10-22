"""Geometry preprocessing utilities for GMSHFlow."""

import geopandas as gpd
import shapely
import shapely.ops
import topojson as tp


def merge_many_multilinestring_into_one_linestring(gdf):
    '''
    Merges multiple MultiLineString geometries in a GeoDataFrame into single LineString geometries.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        A GeoDataFrame containing MultiLineString geometries.

    Raises
    ------
    AssertionError
        If the input GeoDataFrame does not contain any
         MultiLineString geometries or if the merging process fails.
       
    Returns
    -------
    gdf : geopandas.GeoDataFrame
        A GeoDataFrame with merged LineString geometries.

    '''

    assert any(gdf.geom_type == 'MultiLineString'), 'The geodataframe must have multilinestrings'
    # explode the multilinestrings
    gdf2 = gdf.explode().reset_index()
    # group the linestrings by the index
    gdf.geometry = gdf2.groupby('index')['geometry'].apply(
        lambda x: shapely.ops.linemerge(list(x), directed=True))
    #check that there are no Multilinestrings anymore and, that the number of lines is the same as the original
    assert all(gdf.geom_type == 'LineString'), 'The geodataframe must have only linestrings, then the algorithm didnt work'

    return gdf


def simplify_keeping_topology(gdf, cs, plot=False):
    '''
    Simplifies the geometries in a GeoDataFrame while keeping the topology of the original geometries.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        A GeoDataFrame containing the geometries to be simplified.
    cs : float
        The cell size to be used for simplifying the geometries.
    plot : bool, optional
        If True, the simplified geometries are plotted. The default is False.
    
    Returns
    -------
    gdf : geopandas.GeoDataFrame
        A GeoDataFrame with the simplified geometries.
    
    '''
    topo = tp.Topology(gdf.to_crs(gdf.crs), prequantize=False)
    simple = topo.toposimplify(cs/2).to_gdf()
    if plot:
        simple.plot()
    gdf.geometry = simple.geometry
    return gdf