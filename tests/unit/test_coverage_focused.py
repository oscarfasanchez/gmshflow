"""Focused unit tests for improving coverage on key code branches."""

import sys
from pathlib import Path
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

def test_line_handler_basic_validation():
    """Test basic LineGeometryHandler validation without complex mocking."""
    import geopandas as gpd
    from shapely.geometry import Point

    from gmshflow.geometry.line import LineGeometryHandler

    handler = LineGeometryHandler()

    # Test with invalid geometry type
    point_gdf = gpd.GeoDataFrame({'cs': [1.0]}, geometry=[Point(0, 0)])
    with pytest.raises(ValueError, match="All geometries must be LineString or MultiLineString"):
        handler.set_gpd_line(point_gdf)

def test_line_handler_missing_cs_validation():
    """Test LineGeometryHandler cs column validation."""
    import geopandas as gpd
    from shapely.geometry import LineString

    from gmshflow.geometry.line import LineGeometryHandler

    handler = LineGeometryHandler()

    # Create line without cs column and no cs_line parameter
    line_gdf = gpd.GeoDataFrame({'id': [1]}, geometry=[LineString([(0, 0), (1, 1)])])

    with pytest.raises(ValueError, match="Either 'cs' column must exist in GeoDataFrame or cs_line must be set"):
        handler.set_gpd_line(line_gdf)

def test_line_handler_cs_line_parameter():
    """Test LineGeometryHandler with cs_line parameter."""
    import geopandas as gpd
    from shapely.geometry import LineString

    from gmshflow.geometry.line import LineGeometryHandler

    handler = LineGeometryHandler(cs_line=5.0)
    line_gdf = gpd.GeoDataFrame({'id': [1]}, geometry=[LineString([(0, 0), (1, 1)])])

    with patch('builtins.print'):  # Suppress warnings
        handler.set_gpd_line(line_gdf)
        assert 'cs' in handler.gdf_line.columns
        assert handler.gdf_line['cs'].iloc[0] == 5.0

def test_polygon_handler_basic_validation():
    """Test basic PolyGeometryHandler validation."""
    import geopandas as gpd
    from shapely.geometry import Point

    from gmshflow.geometry.polygon import PolyGeometryHandler

    handler = PolyGeometryHandler()

    # Test with invalid geometry type
    point_gdf = gpd.GeoDataFrame({'cs': [1.0]}, geometry=[Point(0, 0)])
    with pytest.raises(ValueError, match="All geometries must be Polygon type"):
        handler.set_gpd_poly(point_gdf)

def test_polygon_handler_missing_cs_validation():
    """Test PolyGeometryHandler cs column validation."""
    import geopandas as gpd
    from shapely.geometry import Polygon

    from gmshflow.geometry.polygon import PolyGeometryHandler

    handler = PolyGeometryHandler()

    # Create polygon without cs column and no cs_poly parameter
    poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    poly_gdf = gpd.GeoDataFrame({'id': [1]}, geometry=[poly])

    with pytest.raises(ValueError, match="Either 'cs' column must exist in GeoDataFrame or cs_poly must be set"):
        handler.set_gpd_poly(poly_gdf)

def test_polygon_handler_cs_poly_parameter():
    """Test PolyGeometryHandler with cs_poly parameter."""
    import geopandas as gpd
    from shapely.geometry import Polygon

    from gmshflow.geometry.polygon import PolyGeometryHandler

    handler = PolyGeometryHandler(cs_poly=10.0)
    poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    poly_gdf = gpd.GeoDataFrame({'id': [1]}, geometry=[poly])

    handler.set_gpd_poly(poly_gdf)
    assert 'cs' in handler.gdf_poly.columns
    assert handler.gdf_poly['cs'].iloc[0] == 10.0

def test_point_handler_basic_validation():
    """Test basic PointGeometryHandler validation."""
    import geopandas as gpd
    from shapely.geometry import Polygon

    from gmshflow.geometry.point import PointGeometryHandler

    handler = PointGeometryHandler()

    # Test with invalid geometry type
    poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    poly_gdf = gpd.GeoDataFrame({'cs': [1.0]}, geometry=[poly])
    with pytest.raises(ValueError, match="All geometries must be Point type"):
        handler.set_gdf_point(poly_gdf)

def test_point_handler_missing_cs_validation():
    """Test PointGeometryHandler cs column validation."""
    import geopandas as gpd
    from shapely.geometry import Point

    from gmshflow.geometry.point import PointGeometryHandler

    handler = PointGeometryHandler()

    # Create point without cs column and no cs_point parameter
    point_gdf = gpd.GeoDataFrame({'id': [1]}, geometry=[Point(0, 0)])

    with pytest.raises(ValueError, match="Either 'cs' column must exist in GeoDataFrame or cs_point must be set"):
        handler.set_gdf_point(point_gdf)

def test_point_handler_cs_point_parameter():
    """Test PointGeometryHandler with cs_point parameter."""
    import geopandas as gpd
    from shapely.geometry import Point

    from gmshflow.geometry.point import PointGeometryHandler

    handler = PointGeometryHandler(cs_point=7.5)
    point_gdf = gpd.GeoDataFrame({'id': [1]}, geometry=[Point(0, 0)])

    handler.set_gdf_point(point_gdf)
    assert 'cs' in handler.gdf_point.columns
    assert handler.gdf_point['cs'].iloc[0] == 7.5

def test_domain_initialization():
    """Test GmshMeshDomain initialization."""
    import geopandas as gpd
    from shapely.geometry import Polygon

    from gmshflow.core.domain import GmshMeshDomain

    poly = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    domain_gdf = gpd.GeoDataFrame({'cs': [5.0]}, geometry=[poly])

    domain = GmshMeshDomain("test_domain", domain_gdf, cs_dom=2.5)

    assert domain.name == "test_domain"
    assert domain.cs_dom == 2.5
    assert len(domain.gdf_dom) == 1
    assert domain.gdf_list == []

def test_domain_add_polygon_geometry_type_validation():
    """Test domain polygon geometry type validation."""
    import geopandas as gpd
    from shapely.geometry import Point, Polygon

    from gmshflow.core.domain import GmshMeshDomain

    poly = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    domain_gdf = gpd.GeoDataFrame({'cs': [5.0]}, geometry=[poly])
    domain = GmshMeshDomain("test", domain_gdf, cs_dom=10.0)

    # Test invalid input type
    with pytest.raises(TypeError, match="Expected GeoDataFrame"):
        domain.add_domain_polygon_geometry("not a geodataframe")

    # Test invalid geometry types
    point_gdf = gpd.GeoDataFrame({'cs': [1.0]}, geometry=[Point(0, 0)])
    with pytest.raises(ValueError, match="All geometries must be Polygon type"):
        domain.add_domain_polygon_geometry(point_gdf)

def test_domain_add_valid_polygon_geometry():
    """Test adding valid polygon geometry to domain."""
    import geopandas as gpd
    from shapely.geometry import Polygon

    from gmshflow.core.domain import GmshMeshDomain

    poly = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    domain_gdf = gpd.GeoDataFrame({'cs': [5.0]}, geometry=[poly])
    domain = GmshMeshDomain("test", domain_gdf, cs_dom=10.0)

    # Add valid polygon
    additional_poly = gpd.GeoDataFrame({'cs': [3.0]}, geometry=[poly])
    domain.add_domain_polygon_geometry(additional_poly)

    assert len(domain.gdf_list) == 1
    assert domain.gdf_list[0] is additional_poly

def test_model_initialization():
    """Test GmshModel initialization without context manager."""
    from gmshflow.core.model import GmshModel

    model = GmshModel("test_model")
    assert model.name == "test_model"
    assert not model._initialized
    # Model should not be initialized yet
    assert hasattr(model, '_initialized')

def test_model_error_messages():
    """Test GmshModel error messages for uninitialized operations."""
    from gmshflow.core.model import GmshModel

    model = GmshModel("test_model")

    # Test operations without context manager
    with pytest.raises(RuntimeError, match="is not initialized"):
        model.synchronize()

    with pytest.raises(RuntimeError, match="is not initialized"):
        model.generate_mesh(dimension=2)

    with pytest.raises(RuntimeError, match="is not initialized"):
        model.write("test.msh")

    with pytest.raises(RuntimeError, match="is not initialized"):
        model.get_triangular_quality()

def test_model_finalize_safety():
    """Test that model finalize is safe to call multiple times."""
    from gmshflow.core.model import GmshModel

    model = GmshModel("test_model")

    # Should not raise error when finalizing uninitialized model
    model.finalize()
    model.finalize()  # Should be safe to call multiple times

def test_buffer_thickness_validation():
    """Test buffer thickness validation across geometry handlers."""
    import geopandas as gpd
    from shapely.geometry import LineString, Polygon

    from gmshflow.geometry.line import LineGeometryHandler
    from gmshflow.geometry.polygon import PolyGeometryHandler

    # Test line handler
    line_handler = LineGeometryHandler()
    line_gdf = gpd.GeoDataFrame({'cs': [2.0]}, geometry=[LineString([(0, 0), (1, 1)])])

    with patch('builtins.print'):
        line_handler.set_gpd_line(line_gdf)

    with pytest.raises(ValueError, match="Buffer thickness must be 1 or 2 cells"):
        line_handler.create_surfacegrid_from_buffer_line(cs_thick=3)

    with pytest.raises(ValueError, match="Buffer thickness must be 1 or 2 cells"):
        line_handler.create_surfacegrid_from_buffer_line(cs_thick=0.5)

    # Test polygon handler
    poly_handler = PolyGeometryHandler()
    poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    poly_gdf = gpd.GeoDataFrame({'cs': [2.0]}, geometry=[poly])

    poly_handler.set_gpd_poly(poly_gdf)

    with pytest.raises(ValueError, match="Buffer thickness must be 1 or 2 cells"):
        poly_handler.create_surfacegrid_from_buffer_poly(cs_thick=3)

    with pytest.raises(ValueError, match="Buffer thickness must be 1 or 2 cells"):
        poly_handler.create_surfacegrid_from_buffer_poly(cs_thick=0.5)

def test_preprocessing_functions_exist():
    """Test that preprocessing functions exist."""
    from gmshflow.utils import preprocessing

    # Test that the preprocessing module has expected functions
    assert hasattr(preprocessing, 'simplify_keeping_topology')

    # Test simplify_keeping_topology can be imported
    from gmshflow.utils.preprocessing import simplify_keeping_topology
    assert callable(simplify_keeping_topology)

def test_domain_prepare_mesh_branches():
    """Test different mesh area branches in prepare_mesh_domain."""
    import geopandas as gpd
    from shapely.geometry import Polygon

    from gmshflow.core.domain import GmshMeshDomain

    poly = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    domain_gdf = gpd.GeoDataFrame({'cs': [5.0]}, geometry=[poly])
    domain = GmshMeshDomain("test", domain_gdf, cs_dom=10.0)

    # Test different mesh_area values
    domain.prepare_mesh_domain(mesh_area=0)  # convex hull
    assert hasattr(domain, 'shp_dom')

    domain.prepare_mesh_domain(mesh_area=1)  # oriented envelope
    assert hasattr(domain, 'shp_dom')

    domain.prepare_mesh_domain(mesh_area=2)  # bounding box
    assert hasattr(domain, 'shp_dom')

def test_line_multilinestring_explode():
    """Test LineGeometryHandler MultiLineString explosion."""
    import geopandas as gpd
    from shapely.geometry import LineString, MultiLineString

    from gmshflow.geometry.line import LineGeometryHandler

    # Create MultiLineString
    line1 = LineString([(0, 0), (1, 1)])
    line2 = LineString([(1, 1), (2, 0)])
    multi_line = MultiLineString([line1, line2])

    gdf = gpd.GeoDataFrame({'cs': [2.0]}, geometry=[multi_line])

    handler = LineGeometryHandler()
    with patch('builtins.print'):
        handler.set_gpd_line(gdf)

    # Should be exploded into separate LineStrings
    assert len(handler.gdf_line) == 2
    assert all(handler.gdf_line.geom_type == 'LineString')

def test_geometry_simplification_paths():
    """Test different geometry simplification code paths."""
    import geopandas as gpd
    from shapely.geometry import LineString

    from gmshflow.geometry.line import LineGeometryHandler

    # Test with cs_line parameter (should use global simplification)
    handler = LineGeometryHandler(cs_line=5.0)
    line_gdf = gpd.GeoDataFrame({'id': [1]}, geometry=[LineString([(0, 0), (5, 5)])])

    with patch('builtins.print'):
        handler.set_gpd_line(line_gdf)
        assert 'cs' in handler.gdf_line.columns
        assert handler.gdf_line['cs'].iloc[0] == 5.0

    # Test without cs_line (should use per-feature simplification)
    handler2 = LineGeometryHandler()
    line_gdf2 = gpd.GeoDataFrame({'cs': [3.0]}, geometry=[LineString([(0, 0), (3, 3)])])

    with patch('builtins.print'):
        handler2.set_gpd_line(line_gdf2)
        assert handler2.gdf_line['cs'].iloc[0] == 3.0
