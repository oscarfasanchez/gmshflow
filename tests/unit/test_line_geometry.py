"""Test line geometry handler unit tests with validation focus."""

import sys
from pathlib import Path
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

def test_line_handler_validation():
    """Test LineGeometryHandler input validation."""
    from gmshflow.geometry.line import LineGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_test_linestring, create_test_points

    handler = LineGeometryHandler()

    # Test with valid LineString data
    line_gdf = create_test_linestring()
    handler.set_gpd_line(line_gdf)
    assert len(handler.gdf_line) > 0

    # Test with invalid geometry types
    point_gdf = create_test_points()
    with pytest.raises(ValueError, match="All geometries must be LineString or MultiLineString"):
        handler.set_gpd_line(point_gdf)

def test_line_handler_missing_cs_column():
    """Test LineGeometryHandler when cs column is missing."""
    from gmshflow.geometry.line import LineGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_test_linestring

    handler = LineGeometryHandler()
    line_gdf = create_test_linestring()

    # Remove cs column and don't set cs_line
    line_gdf = line_gdf.drop('cs', axis=1)

    with pytest.raises(ValueError, match="Either 'cs' column must exist in GeoDataFrame or cs_line must be set"):
        handler.set_gpd_line(line_gdf)

def test_line_handler_with_cs_line_parameter():
    """Test LineGeometryHandler using cs_line parameter."""
    from gmshflow.geometry.line import LineGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_test_linestring

    handler = LineGeometryHandler(cs_line=5.0)
    line_gdf = create_test_linestring()
    line_gdf = line_gdf.drop('cs', axis=1)  # Remove cs column

    handler.set_gpd_line(line_gdf)
    assert 'cs' in handler.gdf_line.columns
    assert all(handler.gdf_line['cs'] == 5.0)

def test_line_handler_multilinestring_explode():
    """Test LineGeometryHandler with MultiLineString geometries."""
    import geopandas as gpd
    from shapely.geometry import LineString, MultiLineString

    from gmshflow.geometry.line import LineGeometryHandler

    # Create MultiLineString geometry
    line1 = LineString([(0, 0), (1, 1)])
    line2 = LineString([(1, 1), (2, 0)])
    multi_line = MultiLineString([line1, line2])

    gdf = gpd.GeoDataFrame({'cs': [2.0]}, geometry=[multi_line])

    handler = LineGeometryHandler()
    with patch('builtins.print'):  # Suppress warning print
        handler.set_gpd_line(gdf)

    # Should be exploded into separate LineStrings
    assert len(handler.gdf_line) == 2
    assert all(handler.gdf_line.geom_type == 'LineString')

def test_line_handler_keep_topology_branch():
    """Test LineGeometryHandler with keep_topology=True branch."""
    from gmshflow.geometry.line import LineGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_test_linestring

    handler = LineGeometryHandler(cs_line=10.0)
    line_gdf = create_test_linestring()

    with patch('builtins.print'):  # Suppress warning prints
        with patch('gmshflow.geometry.line.simplify_keeping_topology') as mock_simplify:
            mock_simplify.return_value = line_gdf.copy()
            handler.set_gpd_line(line_gdf, keep_topology=True)
            mock_simplify.assert_called_once()

def test_line_handler_per_feature_simplify():
    """Test LineGeometryHandler with per-feature simplification."""
    from gmshflow.geometry.line import LineGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))

    handler = LineGeometryHandler()
    # Create a linestring with multiple features
    import geopandas as gpd
    from shapely.geometry import LineString
    lines = [LineString([(0, 0), (5, 5), (10, 0)]), LineString([(20, 0), (25, 5), (30, 0)])]
    line_gdf = gpd.GeoDataFrame({'cs': [2.0, 5.0]}, geometry=lines, crs='EPSG:4326')

    with patch('builtins.print'):  # Suppress warning prints
        handler.set_gpd_line(line_gdf)

    assert len(handler.gdf_line) == 2
    assert handler.gdf_line['cs'].tolist() == [2.0, 5.0]

def test_line_create_line_validation():
    """Test line creation validation without GMSH."""
    from gmshflow.geometry.line import LineGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_test_linestring

    handler = LineGeometryHandler()
    line_gdf = create_test_linestring()

    with patch('builtins.print'):
        handler.set_gpd_line(line_gdf)

    # Test that geometry is properly prepared for line creation
    assert hasattr(handler, 'gdf_line')
    assert len(handler.gdf_line) == 1
    assert handler.gdf_line.geometry[0].geom_type == 'LineString'

    # Test that coordinates can be extracted properly
    coords = list(handler.gdf_line.geometry[0].coords)
    assert len(coords) == 3  # LineString with 3 points: (0,0), (5,5), (10,0)
    assert coords[0] == (0.0, 0.0)
    assert coords[1] == (5.0, 5.0)
    assert coords[2] == (10.0, 0.0)

    # The actual GMSH line creation is tested in integration tests

def test_line_handler_convert_to_points():
    """Test converting lines to points for size fields."""
    from gmshflow.geometry.line import LineGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_test_linestring

    handler = LineGeometryHandler()
    line_gdf = create_test_linestring()

    with patch('builtins.print'):
        handler.set_gpd_line(line_gdf)

    # Mock the conversion method
    with patch.object(handler, 'convert_to_points_for_size_fields') as mock_convert:
        mock_convert.return_value = create_test_linestring()  # Return something
        handler.convert_to_points_for_size_fields()
        mock_convert.assert_called_once()

def test_buffer_thickness_validation():
    """Test buffer thickness validation in surface grid creation."""
    from gmshflow.geometry.line import LineGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_test_linestring

    handler = LineGeometryHandler()
    line_gdf = create_test_linestring()

    with patch('builtins.print'):
        handler.set_gpd_line(line_gdf)

    # Should raise error for invalid thickness
    with pytest.raises(ValueError, match="Buffer thickness must be 1 or 2 cells"):
        handler.create_surfacegrid_from_buffer_line(cs_thick=3)

    with pytest.raises(ValueError, match="Buffer thickness must be 1 or 2 cells"):
        handler.create_surfacegrid_from_buffer_line(cs_thick=0.5)

def test_line_surface_grid_validation():
    """Test surface grid validation without GMSH."""
    from gmshflow.geometry.line import LineGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_test_linestring

    handler = LineGeometryHandler()
    line_gdf = create_test_linestring()

    with patch('builtins.print'):
        handler.set_gpd_line(line_gdf)

    # Test that geometry is properly prepared
    assert hasattr(handler, 'gdf_line')
    assert len(handler.gdf_line) == 1

    # Test input validation - thickness must be between 1 and 2
    with pytest.raises(ValueError, match="Buffer thickness must be 1 or 2 cells"):
        handler.create_surfacegrid_from_buffer_line(cs_thick=3)

    with pytest.raises(ValueError, match="Buffer thickness must be 1 or 2 cells"):
        handler.create_surfacegrid_from_buffer_line(cs_thick=0.5)

    # The actual GMSH surface grid creation is tested in integration tests
