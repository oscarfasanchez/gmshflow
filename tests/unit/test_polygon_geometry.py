"""Test polygon geometry handler unit tests with validation focus."""

import sys
from pathlib import Path
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

def test_polygon_handler_validation():
    """Test PolygonGeometryHandler input validation."""
    from gmshflow.geometry.polygon import PolyGeometryHandler as PolygonGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_domain_polygon, create_test_points

    handler = PolygonGeometryHandler()

    # Test with valid Polygon data
    poly_gdf = create_domain_polygon()
    handler.set_gpd_poly(poly_gdf)
    assert len(handler.gdf_poly) > 0

    # Test with invalid geometry types
    point_gdf = create_test_points()
    with pytest.raises(ValueError, match="All geometries must be Polygon type"):
        handler.set_gpd_poly(point_gdf)

def test_polygon_handler_missing_cs_column():
    """Test PolygonGeometryHandler when cs column is missing."""
    from gmshflow.geometry.polygon import PolyGeometryHandler as PolygonGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_domain_polygon

    handler = PolygonGeometryHandler()
    poly_gdf = create_domain_polygon()

    # Remove cs column and don't set cs_poly
    poly_gdf = poly_gdf.drop('cs', axis=1, errors='ignore')

    with pytest.raises(ValueError, match="Either 'cs' column must exist in GeoDataFrame or cs_poly must be set"):
        handler.set_gpd_poly(poly_gdf)

def test_polygon_handler_with_cs_poly_parameter():
    """Test PolygonGeometryHandler using cs_poly parameter."""
    from gmshflow.geometry.polygon import PolyGeometryHandler as PolygonGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_domain_polygon

    handler = PolygonGeometryHandler(cs_poly=10.0)
    poly_gdf = create_domain_polygon()
    poly_gdf = poly_gdf.drop('cs', axis=1, errors='ignore')  # Remove cs column
    handler.set_gpd_poly(poly_gdf)
    assert 'cs' in handler.gdf_poly.columns
    assert all(handler.gdf_poly['cs'] == 10.0)

def test_polygon_handler_keep_topology_branch():
    """Test PolygonGeometryHandler with keep_topology=True branch."""
    from gmshflow.geometry.polygon import PolyGeometryHandler as PolygonGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_domain_polygon

    handler = PolygonGeometryHandler(cs_poly=15.0)
    poly_gdf = create_domain_polygon()

    with patch('builtins.print'):  # Suppress warning prints
        with patch('gmshflow.geometry.polygon.simplify_keeping_topology') as mock_simplify:
            mock_simplify.return_value = poly_gdf.copy()
            handler.set_gpd_poly(poly_gdf, keep_topology=True)
            mock_simplify.assert_called_once()

def test_polygon_handler_per_feature_simplify():
    """Test PolygonGeometryHandler with per-feature simplification."""
    from gmshflow.geometry.polygon import PolyGeometryHandler as PolygonGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))

    handler = PolygonGeometryHandler()
    # Create multiple polygons for per-feature testing
    import geopandas as gpd
    from shapely.geometry import Polygon
    polygons = [
        Polygon([(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)]),
        Polygon([(20, 0), (30, 0), (30, 10), (20, 10), (20, 0)])
    ]
    poly_gdf = gpd.GeoDataFrame({'cs': [5.0, 8.0]}, geometry=polygons, crs='EPSG:4326')

    with patch('builtins.print'):  # Suppress warning prints
        handler.set_gpd_poly(poly_gdf)

    assert len(handler.gdf_poly) >= 1  # At least one polygon
    assert 'cs' in handler.gdf_poly.columns

def test_polygon_create_loop_validation():
    """Test polygon loop creation validation without GMSH."""
    from gmshflow.geometry.polygon import PolyGeometryHandler as PolygonGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_domain_polygon

    handler = PolygonGeometryHandler()
    poly_gdf = create_domain_polygon()

    with patch('builtins.print'):
        handler.set_gpd_poly(poly_gdf)

    # Test that geometry is properly prepared for loop creation
    assert hasattr(handler, 'gdf_poly')
    assert len(handler.gdf_poly) == 1
    assert handler.gdf_poly.geometry[0].geom_type == 'Polygon'

    # Test that coordinates can be extracted properly
    coords = list(handler.gdf_poly.geometry[0].exterior.coords)
    assert len(coords) == 5  # Polygon should have 5 coordinates (including closed loop)
    assert coords[0] == coords[-1]  # First and last coordinates should be the same (closed)

    # The actual GMSH loop creation is tested in integration tests

def test_polygon_buffer_thickness_validation():
    """Test buffer thickness validation in surface grid creation."""
    from gmshflow.geometry.polygon import PolyGeometryHandler as PolygonGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_domain_polygon

    handler = PolygonGeometryHandler()
    poly_gdf = create_domain_polygon()

    with patch('builtins.print'):
        handler.set_gpd_poly(poly_gdf)

    # Should raise error for invalid thickness
    with pytest.raises(ValueError, match="Buffer thickness must be 1 or 2 cells"):
        handler.create_surfacegrid_from_buffer_poly(cs_thick=3)

    with pytest.raises(ValueError, match="Buffer thickness must be 1 or 2 cells"):
        handler.create_surfacegrid_from_buffer_poly(cs_thick=0.5)

def test_polygon_surface_grid_validation():
    """Test surface grid validation without GMSH."""
    from gmshflow.geometry.polygon import PolyGeometryHandler as PolygonGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_domain_polygon

    handler = PolygonGeometryHandler(cs_poly=10.0)
    poly_gdf = create_domain_polygon()

    with patch('builtins.print'):
        handler.set_gpd_poly(poly_gdf)

    # Test that geometry is properly prepared
    assert hasattr(handler, 'gdf_poly')
    assert len(handler.gdf_poly) == 1
    assert handler.cs_poly == 10.0

    # Test input validation - thickness must be between 1 and 2
    with pytest.raises(ValueError, match="Buffer thickness must be 1 or 2 cells"):
        handler.create_surfacegrid_from_buffer_poly(cs_thick=3)

    with pytest.raises(ValueError, match="Buffer thickness must be 1 or 2 cells"):
        handler.create_surfacegrid_from_buffer_poly(cs_thick=0.5)
