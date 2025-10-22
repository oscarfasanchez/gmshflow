"""Sample geometries for testing."""

import geopandas as gpd
from shapely.geometry import LineString, Point, Polygon


def create_test_polygon():
    """Create a simple test polygon."""
    coords = [(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)]
    poly = Polygon(coords)
    gdf = gpd.GeoDataFrame({'cs': [5.0]}, geometry=[poly], crs='EPSG:4326')
    return gdf

def create_test_linestring():
    """Create a simple test linestring."""
    coords = [(0, 0), (5, 5), (10, 0)]
    line = LineString(coords)
    gdf = gpd.GeoDataFrame({'cs': [2.0]}, geometry=[line], crs='EPSG:4326')
    return gdf

def create_test_points():
    """Create test points."""
    points = [Point(1, 1), Point(5, 5), Point(9, 9)]
    gdf = gpd.GeoDataFrame({'cs': [1.0, 1.5, 2.0]}, geometry=points, crs='EPSG:4326')
    return gdf

def create_domain_polygon():
    """Create a domain polygon for mesh testing."""
    coords = [(0, 0), (100, 0), (100, 100), (0, 100), (0, 0)]
    poly = Polygon(coords)
    gdf = gpd.GeoDataFrame({'name': ['domain']}, geometry=[poly], crs='EPSG:4326')
    return gdf
