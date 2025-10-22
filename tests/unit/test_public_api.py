"""Unit tests for public API imports."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

def test_public_api_imports():
    """Test that public API can be imported successfully."""
    from gmshflow import (
        GmshMeshDomain,
        GmshModel,
        LineGeometryHandler,
        PointGeometryHandler,
        PolyGeometryHandler,
        merge_many_multilinestring_into_one_linestring,
        simplify_keeping_topology,
    )

    # Basic instantiation tests (no GMSH calls)
    assert GmshModel is not None
    assert GmshMeshDomain is not None
    assert PolyGeometryHandler is not None
    assert LineGeometryHandler is not None
    assert PointGeometryHandler is not None

    # Function availability
    assert callable(merge_many_multilinestring_into_one_linestring)
    assert callable(simplify_keeping_topology)

def test_class_instantiation():
    """Test basic class instantiation without GMSH calls."""
    from gmshflow import (
        GmshMeshDomain,
        LineGeometryHandler,
        PointGeometryHandler,
        PolyGeometryHandler,
    )
    sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
    from sample_geometries import create_domain_polygon

    # These should work without GMSH
    domain_gdf = create_domain_polygon()
    domain = GmshMeshDomain("test_domain", domain_gdf, 10.0)
    assert domain.name == "test_domain"
    assert domain.cs_dom == 10.0

    poly_handler = PolyGeometryHandler(cs_poly=5.0)
    assert poly_handler.cs_poly == 5.0

    line_handler = LineGeometryHandler(cs_line=2.0)
    assert line_handler.cs_line == 2.0

    point_handler = PointGeometryHandler(cs_point=1.0)
    assert point_handler.cs_point == 1.0
