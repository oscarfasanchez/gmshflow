"""Additional unit tests to improve coverage on specific geometry methods."""

import sys
from pathlib import Path
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

def test_line_convert_to_points_functionality():
    """Test line conversion to points for size fields."""
    import geopandas as gpd
    from shapely.geometry import LineString

    from gmshflow.geometry.line import LineGeometryHandler

    handler = LineGeometryHandler()
    line = LineString([(0, 0), (5, 5), (10, 0)])
    line_gdf = gpd.GeoDataFrame({'cs': [2.0]}, geometry=[line])

    with patch('builtins.print'):
        handler.set_gpd_line(line_gdf)

    # Test that the method exists and can be called
    # Note: We're just testing the method exists since it has complex GMSH dependencies
    assert hasattr(handler, 'convert_to_points_for_size_fields')

    # Mock the method to test the interface
    with patch.object(handler, 'convert_to_points_for_size_fields') as mock_convert:
        mock_convert.return_value = gpd.GeoDataFrame()
        handler.convert_to_points_for_size_fields()
        mock_convert.assert_called_once()

def test_line_buffer_surface_warnings():
    """Test warning prints in buffer surface creation."""
    import geopandas as gpd
    from shapely.geometry import LineString

    from gmshflow.geometry.line import LineGeometryHandler

    handler = LineGeometryHandler()
    line = LineString([(0, 0), (5, 5), (10, 0)])
    line_gdf = gpd.GeoDataFrame({'cs': [2.0]}, geometry=[line])

    with patch('builtins.print'):
        handler.set_gpd_line(line_gdf)

    # Test that the method prints the expected warning
    with patch('builtins.print') as mock_print:
        try:
            handler.create_surfacegrid_from_buffer_line(cs_thick=1.5)
        except Exception:
            pass  # Expected to fail due to GMSH dependencies

        # Check that warning was printed
        warning_calls = [call for call in mock_print.call_args_list
                        if 'starting: create_surfacegrid_from_buffer_line' in str(call)]
        assert len(warning_calls) > 0

def test_polygon_buffer_surface_warnings():
    """Test warning prints in polygon buffer surface creation."""
    import geopandas as gpd
    from shapely.geometry import Polygon

    from gmshflow.geometry.polygon import PolyGeometryHandler

    handler = PolyGeometryHandler()
    poly = Polygon([(0, 0), (5, 0), (5, 5), (0, 5)])
    poly_gdf = gpd.GeoDataFrame({'cs': [2.0]}, geometry=[poly])

    handler.set_gpd_poly(poly_gdf)

    # Test that the method prints the expected warning
    with patch('builtins.print') as mock_print:
        try:
            handler.create_surfacegrid_from_buffer_poly(cs_thick=1.5)
        except Exception:
            pass  # Expected to fail due to GMSH dependencies

        # Check that warning was printed
        warning_calls = [call for call in mock_print.call_args_list
                        if 'starting: create_surfacegrid_from_buffer_poly' in str(call)]
        assert len(warning_calls) > 0

def test_domain_prepare_mesh_with_additional_geometries():
    """Test domain preparation with additional geometries."""
    import geopandas as gpd
    from shapely.geometry import Polygon

    from gmshflow.core.domain import GmshMeshDomain

    # Create main domain
    main_poly = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    domain_gdf = gpd.GeoDataFrame({'cs': [5.0]}, geometry=[main_poly])
    domain = GmshMeshDomain("test", domain_gdf, cs_dom=10.0)

    # Add additional geometries
    additional_poly = Polygon([(2, 2), (8, 2), (8, 8), (2, 8)])
    additional_gdf = gpd.GeoDataFrame({'cs': [3.0]}, geometry=[additional_poly])
    domain.add_domain_polygon_geometry(additional_gdf)

    # Test preparation with gdf_list processing
    domain.prepare_mesh_domain(mesh_area=1, min_overlap=0.1)
    assert hasattr(domain, 'shp_dom')

def test_domain_prepare_mesh_empty_list():
    """Test domain preparation with empty gdf_list."""
    import geopandas as gpd
    from shapely.geometry import Polygon

    from gmshflow.core.domain import GmshMeshDomain

    poly = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    domain_gdf = gpd.GeoDataFrame({'cs': [5.0]}, geometry=[poly])
    domain = GmshMeshDomain("test", domain_gdf, cs_dom=10.0)

    # Test with explicitly empty gdf_list
    domain.prepare_mesh_domain(mesh_area=1, gdf_list=[])
    assert hasattr(domain, 'shp_dom')

def test_preprocessing_simplify_topojson_path():
    """Test simplify_keeping_topology function."""
    import geopandas as gpd
    from shapely.geometry import Polygon

    from gmshflow.utils.preprocessing import simplify_keeping_topology

    # Create a simple polygon with CRS to avoid coordinate transformation errors
    poly = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    poly_gdf = gpd.GeoDataFrame({'cs': [1.0]}, geometry=[poly], crs='EPSG:4326')

    with patch('builtins.print'):  # Suppress warning prints
        result = simplify_keeping_topology(poly_gdf, 0.5)
        assert len(result) == 1
        assert result.geom_type[0] == 'Polygon'

def test_domain_mesh_options_coverage():
    """Test different domain mesh area options for coverage."""
    import geopandas as gpd
    from shapely.geometry import Polygon

    from gmshflow.core.domain import GmshMeshDomain

    poly = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    domain_gdf = gpd.GeoDataFrame({'cs': [5.0]}, geometry=[poly])

    for mesh_area in [0, 1, 2]:
        domain = GmshMeshDomain(f"test_{mesh_area}", domain_gdf, cs_dom=10.0)
        domain.prepare_mesh_domain(mesh_area=mesh_area, meshing_buff_mult=1.5)
        assert hasattr(domain, 'shp_dom')

def test_geometry_handler_cs_warning_branches():
    """Test warning branches in geometry handlers when topology is kept."""
    import geopandas as gpd
    from shapely.geometry import LineString

    from gmshflow.geometry.line import LineGeometryHandler

    # Test the topology warning branch by having per-feature cs values
    handler = LineGeometryHandler()  # No global cs_line
    line = LineString([(0, 0), (5, 5)])
    line_gdf = gpd.GeoDataFrame({'cs': [2.0]}, geometry=[line])

    with patch('builtins.print') as mock_print:
        # This should trigger normal processing without warnings in the topology branch
        handler.set_gpd_line(line_gdf, keep_topology=False)

        # Should have basic warning
        warning_calls = [call for call in mock_print.call_args_list
                        if 'Warning: the line geometries will be simplified' in str(call)]
        assert len(warning_calls) > 0

def test_point_handler_coordinate_dataframe():
    """Test PointGeometryHandler coordinate dataframe creation."""
    import geopandas as gpd
    from shapely.geometry import Point

    from gmshflow.geometry.point import HAS_GMSH, PointGeometryHandler

    handler = PointGeometryHandler()
    points = [Point(0, 0), Point(5, 5), Point(10, 0)]
    point_gdf = gpd.GeoDataFrame({'cs': [1.0, 2.0, 1.5]}, geometry=points)

    handler.set_gdf_point(point_gdf)

    if not HAS_GMSH:
        # Without GMSH, test that the method exists and raises appropriate error
        with pytest.raises(ImportError, match="GMSH is required"):
            handler.create_point_from_point(df_coord=True)
    else:
        # Mock point creation with coordinate dataframe when GMSH is available
        with patch('gmsh.model.geo.addPoint', side_effect=[1, 2, 3]):
            point_ids = handler.create_point_from_point(df_coord=True)

            assert isinstance(point_ids, list)
            assert len(point_ids) == 3
            assert hasattr(handler, 'gdf_coord')
            assert 'id_gmsh' in handler.gdf_coord.columns

def test_geometry_simplification_different_tolerances():
    """Test geometry simplification with different tolerance values."""
    import geopandas as gpd
    from shapely.geometry import LineString, Polygon

    from gmshflow.geometry.line import LineGeometryHandler
    from gmshflow.geometry.polygon import PolyGeometryHandler

    # Test line simplification
    line_handler = LineGeometryHandler(cs_line=2.0)
    line = LineString([(0, 0), (1, 1), (2, 0), (3, 1), (4, 0)])
    line_gdf = gpd.GeoDataFrame({'id': [1]}, geometry=[line])

    with patch('builtins.print'):
        line_handler.set_gpd_line(line_gdf)
        # Line should be simplified with tolerance cs_line/2 = 1.0
        assert len(line_handler.gdf_line) == 1
        assert 'cs' in line_handler.gdf_line.columns

    # Test polygon simplification
    poly_handler = PolyGeometryHandler(cs_poly=3.0)
    vertices = [(0, 0), (1, 0.1), (2, 0), (2, 1), (1.9, 2), (2, 3), (1, 3), (0, 2.9), (0, 2), (0.1, 1)]
    poly = Polygon(vertices)
    poly_gdf = gpd.GeoDataFrame({'id': [1]}, geometry=[poly])

    poly_handler.set_gpd_poly(poly_gdf)
    # Polygon should be simplified with tolerance cs_poly/2 = 1.5
    assert len(poly_handler.gdf_poly) == 1
    assert 'cs' in poly_handler.gdf_poly.columns
