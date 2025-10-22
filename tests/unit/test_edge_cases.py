"""Test point geometry handler and model edge cases."""

import pytest
import sys
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

def test_point_handler_validation():
    """Test PointGeometryHandler input validation."""
    from gmshflow.geometry.point import PointGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_test_points, create_domain_polygon
    
    handler = PointGeometryHandler()
    
    # Test with valid Point data
    point_gdf = create_test_points()
    handler.set_gdf_point(point_gdf)
    assert len(handler.gdf_point) > 0
    
    # Test with invalid geometry types  
    poly_gdf = create_domain_polygon()
    with pytest.raises(ValueError, match="All geometries must be Point type"):
        handler.set_gdf_point(poly_gdf)

def test_point_handler_missing_cs_column():
    """Test PointGeometryHandler when cs column is missing."""
    from gmshflow.geometry.point import PointGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_test_points
    
    handler = PointGeometryHandler()
    point_gdf = create_test_points()
    
    # Remove cs column and don't set cs_point
    point_gdf = point_gdf.drop('cs', axis=1)
    
    with pytest.raises(ValueError, match="Either 'cs' column must exist in GeoDataFrame or cs_point must be set"):
        handler.set_gdf_point(point_gdf)

def test_point_handler_with_cs_point_parameter():
    """Test PointGeometryHandler using cs_point parameter."""
    from gmshflow.geometry.point import PointGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_test_points
    
    handler = PointGeometryHandler(cs_point=7.5)
    point_gdf = create_test_points()
    point_gdf = point_gdf.drop('cs', axis=1)  # Remove cs column
    
    handler.set_gdf_point(point_gdf)
    assert 'cs' in handler.gdf_point.columns
    assert all(handler.gdf_point['cs'] == 7.5)

def test_point_create_point_validation():
    """Test point creation validation without GMSH."""
    from gmshflow.geometry.point import PointGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_test_points
    
    handler = PointGeometryHandler()
    point_gdf = create_test_points()
    
    handler.set_gdf_point(point_gdf)
    
    # Test that geometry is properly prepared for point creation
    assert hasattr(handler, 'gdf_point')
    assert len(handler.gdf_point) == len(point_gdf)
    assert all(geom.geom_type == 'Point' for geom in handler.gdf_point.geometry)
    
    # Test that all points have valid coordinates
    for geom in handler.gdf_point.geometry:
        coords = geom.coords[0]
        assert len(coords) == 2  # x, y coordinates
        assert all(isinstance(coord, (int, float)) for coord in coords)
    
    # The actual GMSH point creation is tested in integration tests

def test_point_create_point_with_coord_df():
    """Test creating points with coordinate DataFrame output."""
    from gmshflow.geometry.point import PointGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_test_points
    
    handler = PointGeometryHandler()
    point_gdf = create_test_points()
    
    with patch('gmsh.model.geo.addPoint', side_effect=[1, 2, 3, 4]):
        handler.set_gdf_point(point_gdf)
        point_ids = handler.create_point_from_point(df_coord=True)
        
        assert isinstance(point_ids, list)
        assert hasattr(handler, 'gdf_coord')
        assert 'id_gmsh' in handler.gdf_coord.columns
        assert len(handler.gdf_coord) == len(point_gdf)

def test_model_error_cases():
    """Test GmshModel error handling cases."""
    from gmshflow.core.model import GmshModel
    
    model = GmshModel("test_model")
    
    # Test operations without GMSH context (should raise errors)
    with pytest.raises(RuntimeError, match="not initialized"):
        model.synchronize()
    
    with pytest.raises(RuntimeError, match="is not initialized"):
        model.generate_mesh(dimension=2)
    
    with pytest.raises(RuntimeError, match="is not initialized"):
        model.write("test.msh")

def test_model_finalize_without_init():
    """Test model finalization without initialization."""
    from gmshflow.core.model import GmshModel
    
    model = GmshModel("test_model")
    
    # Should not raise error when finalizing uninitialized model
    model.finalize()
    
    # Should be safe to call multiple times
    model.finalize()
    model.finalize()

@patch('gmsh.model.mesh.getNodes')
@patch('gmsh.model.mesh.getElements')
def test_model_get_triangular_quality_empty_mesh(mock_get_elements, mock_get_nodes):
    """Test get_triangular_quality with empty mesh."""
    from gmshflow.core.model import GmshModel
    
    # Mock empty mesh
    mock_get_nodes.return_value = ([], [])
    mock_get_elements.return_value = ([], [], [])
    
    model = GmshModel("test_model")
    model._initialized = True  # Bypass initialization check (correct attribute)
    
    result = model.get_triangular_quality()
    
    # Should return empty DataFrame with correct columns
    assert len(result) == 0
    expected_columns = ['minSICN', 'minDetJac', 'maxDetJac', 'minSJ', 'minSIGE', 'gamma', 
                       'innerRadius', 'outerRadius', 'minIsotropy', 'angleShape', 'minEdge', 'maxEdge']
    assert list(result.columns) == expected_columns

def test_preprocessing_edge_cases():
    """Test preprocessing utility edge cases."""
    from gmshflow.utils.preprocessing import merge_many_multilinestring_into_one_linestring
    import geopandas as gpd
    from shapely.geometry import LineString
    
    # Test with regular LineString (should raise error)
    line_gdf = gpd.GeoDataFrame({'cs': [1.0]}, geometry=[LineString([(0, 0), (1, 1)])])
    
    with pytest.raises(ValueError, match="MultiLineString geometries"):
        merge_many_multilinestring_into_one_linestring(line_gdf)

def test_preprocessing_simplify_edge_cases():
    """Test simplify_keeping_topology edge cases."""
    from gmshflow.utils.preprocessing import simplify_keeping_topology
    import geopandas as gpd
    from shapely.geometry import Polygon
    
    # Test with very small tolerance (should warn about keeping topology)
    poly_gdf = gpd.GeoDataFrame({'cs': [1.0]}, geometry=[Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])], crs='EPSG:4326')

    with patch('builtins.print'):  # Suppress warning print
        result = simplify_keeping_topology(poly_gdf, 0.001)
        assert len(result) == 1
        assert result.geom_type[0] == 'Polygon'

def test_model_context_manager_validation():
    """Test context manager basic validation without GMSH."""
    from gmshflow.core.model import GmshModel
    
    # Test model creation
    model = GmshModel("test_model")
    assert model.name == "test_model"
    assert not model._initialized
    
    # Test double initialization check
    model._initialized = True
    with pytest.raises(RuntimeError, match="already initialized"):
        model.__enter__()

def test_domain_voronoi_export_edge_cases():
    """Test Voronoi export with edge cases."""
    from gmshflow.core.domain import GmshMeshDomain
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_domain_polygon
    
    domain_gdf = create_domain_polygon()
    domain = GmshMeshDomain("test", domain_gdf, cs_dom=10.0)
    
    # Test with empty surface list
    with pytest.raises(ValueError, match="No valid surfaces provided"):
        domain.export_to_voronoi("/tmp", "test", [], int_org_dom=True)

def test_geometry_simplification_branches():
    """Test different geometry simplification code paths."""
    from gmshflow.geometry.line import LineGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_test_linestring
    
    # Test with no cs_line set and per-feature cs values
    handler = LineGeometryHandler()  # No cs_line set
    line_gdf = create_test_linestring()
    
    with patch('builtins.print'):  # Suppress warnings
        with patch('gmshflow.geometry.line.simplify_keeping_topology') as mock_simplify:
            mock_simplify.return_value = line_gdf.copy()
            
            # This should trigger the per-feature simplification branch
            handler.set_gpd_line(line_gdf, keep_topology=True)  # keep_topology=True with per-feature cs
            
            # Should not call global simplify_keeping_topology since cs_line is None
            # but should simplify per feature
            assert len(handler.gdf_line) > 0

def test_size_field_warning_branches():
    """Test size field creation warning branches.""" 
    from gmshflow.geometry.line import LineGeometryHandler
    import geopandas as gpd
    from shapely.geometry import LineString
    
    # Create linestring with cs column but no cs_line parameter
    line = LineString([(0, 0), (5, 5)])
    line_gdf = gpd.GeoDataFrame({'cs': [2.0, 8.0]}, geometry=[line, line])
    
    handler = LineGeometryHandler()  # No cs_line set
    
    with patch('builtins.print') as mock_print:
        with patch('gmshflow.geometry.line.simplify_keeping_topology') as mock_simplify:
            mock_simplify.return_value = line_gdf.copy()
            
            # This should trigger the warning about keeping topology with biggest cell size
            handler.set_gpd_line(line_gdf, keep_topology=True)
            
            # Check that warning was printed
            warning_calls = [call for call in mock_print.call_args_list 
                           if 'Warning: the topology will be kept' in str(call)]
            assert len(warning_calls) > 0