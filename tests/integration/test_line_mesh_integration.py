"""Test line geometry integration with GMSH mesh generation."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).parent.parent))

from conftest import gmsh_required


@gmsh_required
class TestLineGeometryMesh:
    """Test line geometry integration with mesh generation."""

    def test_mesh_with_embedded_lines(self):
        """Test mesh generation with embedded line geometries."""
        import gmsh

        from gmshflow import GmshMeshDomain, GmshModel, LineGeometryHandler
        sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
        from sample_geometries import create_domain_polygon, create_test_linestring

        # Create test data
        domain_gdf = create_domain_polygon()
        line_gdf = create_test_linestring()

        with GmshModel("embedded_lines_test") as model:
            # Set up domain
            mesh_domain = GmshMeshDomain("test_domain", domain_gdf, cs_dom=25.0)
            mesh_domain.prepare_mesh_domain(mesh_area=1, min_overlap=0.5)

            mesh_domain.create_domain_loop_from_poly()
            ind_s_dom = mesh_domain.create_domain_surface()

            # Add line geometry
            line_handler = LineGeometryHandler()
            line_handler.set_gpd_line(line_gdf)
            line_ids = line_handler.create_line_from_line()

            # Embed lines in surface
            mesh_domain.add_embedded_lines(id_line_list=line_ids, surface_id=ind_s_dom)

            # Generate mesh
            model.synchronize()
            model.generate_mesh(dimension=2)

            # Verify mesh was created
            elements = gmsh.model.mesh.getElements(dim=2)
            assert len(elements[1][0]) > 0, "No triangular elements generated"

            # Check that lines are embedded - verify 1D elements exist
            line_elements = gmsh.model.mesh.getElements(dim=1)
            assert len(line_elements[1]) > 0, "No line elements found - embedding failed"

            print(f"Generated {len(elements[1][0])} triangular elements")
            print(f"Embedded {len(line_elements[1][0])} line elements")

    def test_line_size_field_integration(self):
        """Test line geometry with size field mesh generation."""
        import gmsh

        from gmshflow import GmshMeshDomain, GmshModel, LineGeometryHandler
        sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
        from sample_geometries import create_domain_polygon, create_test_linestring

        # Create test data with fine cell size on line
        domain_gdf = create_domain_polygon()
        line_gdf = create_test_linestring()
        line_gdf['cs'] = 2.0  # Fine mesh size along lines

        with GmshModel("line_field_test") as model:
            # Set up domain with coarse mesh
            mesh_domain = GmshMeshDomain("test_domain", domain_gdf, cs_dom=20.0)
            mesh_domain.prepare_mesh_domain(mesh_area=1, min_overlap=0.5)

            mesh_domain.create_domain_loop_from_poly()
            mesh_domain.create_domain_surface()

            # Create line and convert to points for size field
            line_handler = LineGeometryHandler()
            line_handler.set_gpd_line(line_gdf)

            # Convert lines to size control points
            df_points = line_handler.convert_to_points_for_size_fields()
            assert len(df_points) > 0, "No control points generated from lines"
            assert 'id_gmsh' in df_points.columns, "Missing GMSH point IDs"

            # Create exponential field from line points
            ind_min_field, ind_exp_field = mesh_domain.create_exponential_field(df_points)
            mesh_domain.set_field()
            mesh_domain.set_mesh_size_from_geometries(use_boundaries=False, use_points=False)

            # Generate mesh
            model.synchronize()
            model.generate_mesh(dimension=2)

            # Verify mesh adapts to line size field
            elements = gmsh.model.mesh.getElements(dim=2)
            assert len(elements[1][0]) > 0, "No elements generated"

            # Check mesh refinement near lines by analyzing edge lengths
            nodes = gmsh.model.mesh.getNodes()
            node_coords = nodes[1].reshape(-1, 3)

            # Get element connectivity
            elem_nodes = elements[2][0].reshape(-1, 3)  # Triangular elements

            # Calculate some edge lengths
            edge_lengths = []
            for i in range(min(100, len(elem_nodes))):  # Sample first 100 elements
                elem = elem_nodes[i]
                n1, n2, n3 = elem - 1  # Convert to 0-based indexing
                if n1 < len(node_coords) and n2 < len(node_coords):
                    p1, p2 = node_coords[n1][:2], node_coords[n2][:2]
                    length = ((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)**0.5
                    edge_lengths.append(length)

            if edge_lengths:
                min_edge = min(edge_lengths)
                max_edge = max(edge_lengths)
                print(f"Edge length range: {min_edge:.2f} - {max_edge:.2f}")

                # Should have some fine elements due to line size field
                fine_edges = sum(1 for e in edge_lengths if e < 5.0)
                assert fine_edges > 0, "No fine mesh refinement detected near lines"

@gmsh_required
class TestFieldGeneration:
    """Test mesh field generation capabilities."""

    def test_exponential_field_creation(self):
        """Test exponential size field creation and application."""
        from gmshflow import GmshMeshDomain, GmshModel, PointGeometryHandler
        sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
        from sample_geometries import create_domain_polygon, create_test_points

        domain_gdf = create_domain_polygon()
        points_gdf = create_test_points()

        with GmshModel("exp_field_test") as model:
            mesh_domain = GmshMeshDomain("test_domain", domain_gdf, cs_dom=30.0)
            mesh_domain.prepare_mesh_domain(mesh_area=1, min_overlap=0.5)

            mesh_domain.create_domain_loop_from_poly()
            mesh_domain.create_domain_surface()

            # Create control points for field
            point_handler = PointGeometryHandler()
            point_handler.set_gdf_point(points_gdf)
            point_handler.create_point_from_point(df_coord=True)

            # Create exponential field
            df_coords = point_handler.gdf_coord
            ind_min_field, ind_exp_fields = mesh_domain.create_exponential_field(df_coords, fac=1.2)

            # Verify field creation
            assert ind_min_field is not None, "Minimum field not created"
            assert len(ind_exp_fields) > 0, "No exponential fields created"

            # Apply field and generate mesh
            mesh_domain.set_field()
            mesh_domain.set_mesh_size_from_geometries(use_boundaries=False, use_points=False)

            model.synchronize()
            model.generate_mesh(dimension=2)

            # Verify mesh quality with field
            quality_df = model.get_triangular_quality()
            assert len(quality_df) > 0, "No mesh quality data"

            # Check for reasonable size variation
            min_edge = quality_df['minEdge'].min()
            max_edge = quality_df['maxEdge'].max()
            size_ratio = max_edge / min_edge if min_edge > 0 else float('inf')

            print(f"Size field mesh: {len(quality_df)} elements")
            print(f"Edge size ratio: {size_ratio:.2f}")

            # Should have size variation due to field
            assert size_ratio > 1.5, "Insufficient size variation from field"
            assert size_ratio < 50.0, "Excessive size variation - field too aggressive"

    def test_linear_threshold_field(self):
        """Test linear threshold field creation."""
        import gmsh

        from gmshflow import GmshMeshDomain, GmshModel, PointGeometryHandler
        sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
        from sample_geometries import create_domain_polygon, create_test_points

        domain_gdf = create_domain_polygon()
        points_gdf = create_test_points()

        with GmshModel("linear_field_test") as model:
            mesh_domain = GmshMeshDomain("test_domain", domain_gdf, cs_dom=25.0)
            mesh_domain.prepare_mesh_domain(mesh_area=1, min_overlap=0.5)

            mesh_domain.create_domain_loop_from_poly()
            mesh_domain.create_domain_surface()

            # Create control points
            point_handler = PointGeometryHandler()
            point_handler.set_gdf_point(points_gdf)
            point_handler.create_point_from_point(df_coord=True)

            # Create linear threshold field
            df_coords = point_handler.gdf_coord
            mesh_domain.create_linear_threshold_field(df_coords, fac=1.3)

            # Apply field and generate mesh
            mesh_domain.set_field()
            mesh_domain.set_mesh_size_from_geometries(use_boundaries=False, use_points=False)

            model.synchronize()
            model.generate_mesh(dimension=2)

            # Verify mesh generation with threshold field
            elements = gmsh.model.mesh.getElements(dim=2)
            assert len(elements[1][0]) > 0, "No elements generated with threshold field"

            quality_df = model.get_triangular_quality()
            good_quality = (quality_df['minSICN'] > 0.2).sum()
            total_elements = len(quality_df)

            print(f"Threshold field mesh: {total_elements} elements")
            print(f"Good quality ratio: {good_quality/total_elements:.1%}")

            assert good_quality / total_elements > 0.7, "Too many poor quality elements with threshold field"
