"""Geometry preprocessing utilities for GMSHFlow."""

from typing import Union
import geopandas as gpd
import shapely
import shapely.ops
import topojson as tp


def merge_many_multilinestring_into_one_linestring(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Merges multiple MultiLineString geometries in a GeoDataFrame into single LineString geometries.

    Args:
        gdf: A GeoDataFrame containing MultiLineString geometries.

    Returns:
        A GeoDataFrame with merged LineString geometries.

    Raises:
        AssertionError: If the input GeoDataFrame does not contain any
            MultiLineString geometries or if the merging process fails.

    Example:
        >>> gdf_with_multilines = load_multilinestring_data()
        >>> gdf_merged = merge_many_multilinestring_into_one_linestring(gdf_with_multilines)
        >>> assert all(gdf_merged.geom_type == 'LineString')
    """

    assert any(gdf.geom_type == 'MultiLineString'), 'The geodataframe must have multilinestrings'
    # explode the multilinestrings
    gdf2 = gdf.explode().reset_index()
    # group the linestrings by the index
    gdf.geometry = gdf2.groupby('index')['geometry'].apply(
        lambda x: shapely.ops.linemerge(list(x), directed=True))
    #check that there are no Multilinestrings anymore and, that the number of lines is the same as the original
    assert all(gdf.geom_type == 'LineString'), 'The geodataframe must have only linestrings, then the algorithm didnt work'

    return gdf


def simplify_keeping_topology(gdf: gpd.GeoDataFrame, cs: float, plot: bool = False) -> gpd.GeoDataFrame:
    """Simplifies geometries in a GeoDataFrame while preserving topology.

    Uses TopJSON topology-preserving simplification to reduce geometry complexity
    while maintaining spatial relationships between features.

    Args:
        gdf: A GeoDataFrame containing the geometries to be simplified.
        cs: The cell size used for simplification tolerance (tolerance = cs/2).
        plot: If True, plots the simplified geometries for visualization.

    Returns:
        A GeoDataFrame with simplified geometries that preserve topology.

    Example:
        >>> gdf_complex = load_complex_polygons()
        >>> gdf_simple = simplify_keeping_topology(gdf_complex, cs=100.0)
        >>> # Geometries are simplified but topology is preserved
    """
    topo = tp.Topology(gdf.to_crs(gdf.crs), prequantize=False)
    simple = topo.toposimplify(cs/2).to_gdf()
    if plot:
        simple.plot()
    gdf.geometry = simple.geometry
    return gdf