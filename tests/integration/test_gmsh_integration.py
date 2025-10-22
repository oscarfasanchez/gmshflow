"""Integration tests with real GMSH."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).parent.parent))

from conftest import gmsh_required


@gmsh_required
class TestGmshModelIntegration:
    """Test GmshModel with real GMSH."""

    def test_context_manager_real_gmsh(self):
        """Test context manager with actual GMSH."""
        import gmsh

        from gmshflow.core.model import GmshModel

        # Ensure GMSH is not initialized
        try:
            gmsh.finalize()
        except:
            pass

        with GmshModel("test") as model:
            assert model._initialized
            model.synchronize()  # Real GMSH call

        # GMSH should be finalized
        with pytest.raises(Exception):  # GMSH throws when not initialized
            gmsh.model.add("should_fail")

    def test_simple_mesh_generation(self):
        """Test basic mesh generation."""
        import gmsh

        from gmshflow.core.model import GmshModel

        with GmshModel("square") as model:
            # Create simple square geometry
            p1 = gmsh.model.geo.addPoint(0, 0, 0, 1.0)
            p2 = gmsh.model.geo.addPoint(1, 0, 0, 1.0)
            p3 = gmsh.model.geo.addPoint(1, 1, 0, 1.0)
            p4 = gmsh.model.geo.addPoint(0, 1, 0, 1.0)

            l1 = gmsh.model.geo.addLine(p1, p2)
            l2 = gmsh.model.geo.addLine(p2, p3)
            l3 = gmsh.model.geo.addLine(p3, p4)
            l4 = gmsh.model.geo.addLine(p4, p1)

            loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
            surface = gmsh.model.geo.addPlaneSurface([loop])

            model.synchronize()
            model.generate_mesh(2)

            # Check mesh was created
            quality_df = model.get_triangular_quality()
            assert len(quality_df) > 0
            assert quality_df['minSICN'].min() > 0

@gmsh_required
class TestGeometryHandlerIntegration:
    """Test geometry handlers with real GMSH."""

    def test_polygon_handler_creates_gmsh_entities(self):
        """Test polygon handler with GMSH."""
        from gmshflow import GmshModel, PolyGeometryHandler
        sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
        from sample_geometries import create_test_polygon

        gdf = create_test_polygon()
        handler = PolyGeometryHandler()
        handler.set_gpd_poly(gdf)

        with GmshModel("test") as model:
            loop_ids = handler.create_loop_from_poly()
            model.synchronize()

            # Verify GMSH entities were created
            import gmsh
            curves = gmsh.model.getEntities(1)  # Lines
            assert len(curves) > 0
