"""Integration tests for polygon surface grid creation with real GMSH."""

import sys
from pathlib import Path
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

@pytest.mark.gmsh
def test_polygon_surface_grid_creation_with_gmsh():
    """Test surface grid creation from buffer polygons (requires GMSH)."""
    gmsh = pytest.importorskip("gmsh")
    from gmshflow.geometry.polygon import PolyGeometryHandler as PolygonGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_domain_polygon

    # Initialize GMSH for this test
    gmsh.initialize()
    try:
        gmsh.model.add("polygon_surface_test")

        handler = PolygonGeometryHandler(cs_poly=10.0)
        poly_gdf = create_domain_polygon()

        with patch('builtins.print'):  # Suppress print output
            handler.set_gpd_poly(poly_gdf)

            # Test surface grid creation with real GMSH
            result = handler.create_surfacegrid_from_buffer_poly(cs_thick=1.5)

            # Verify result structure
            assert isinstance(result, tuple)
            assert len(result) == 5  # Should return (gdf_poly, c_ind_buf_pos, c_ind_buf_neg, ind_s_buff, ind_s_mid)

            gdf_poly, c_ind_buf_pos, c_ind_buf_neg, ind_s_buff, ind_s_mid = result

            # Verify returned data types
            assert isinstance(c_ind_buf_pos, list)
            assert isinstance(c_ind_buf_neg, list)
            assert isinstance(ind_s_buff, list)
            assert isinstance(ind_s_mid, list)

    finally:
        gmsh.finalize()

@pytest.mark.gmsh
def test_polygon_surface_grid_buffer_validation():
    """Test surface grid buffer thickness validation with GMSH."""
    pytest.importorskip("gmsh")
    from gmshflow.geometry.polygon import PolyGeometryHandler as PolygonGeometryHandler
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_domain_polygon

    handler = PolygonGeometryHandler(cs_poly=10.0)
    poly_gdf = create_domain_polygon()

    with patch('builtins.print'):
        handler.set_gpd_poly(poly_gdf)

    # Test validation - these should raise errors even with GMSH available
    with pytest.raises(ValueError, match="Buffer thickness must be 1 or 2 cells"):
        handler.create_surfacegrid_from_buffer_poly(cs_thick=3)

    with pytest.raises(ValueError, match="Buffer thickness must be 1 or 2 cells"):
        handler.create_surfacegrid_from_buffer_poly(cs_thick=0.5)
