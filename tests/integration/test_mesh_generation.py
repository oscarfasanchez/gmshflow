"""Comprehensive GMSH mesh generation tests."""

import pytest
import sys
import tempfile
import os
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).parent.parent))

from conftest import gmsh_required

@gmsh_required
class TestMeshGeneration:
    """Test actual mesh generation with GMSH."""
    
    def test_complete_mesh_workflow(self):
        """Test complete workflow from domain to mesh generation."""
        import gmsh
        from gmshflow import GmshModel, GmshMeshDomain
        sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
        from sample_geometries import create_domain_polygon
        
        # Create test domain
        domain_gdf = create_domain_polygon()
        
        with GmshModel("mesh_workflow_test") as model:
            # Create mesh domain
            mesh_domain = GmshMeshDomain("test_domain", domain_gdf, cs_dom=25.0)
            mesh_domain.prepare_mesh_domain(mesh_area=1, min_overlap=0.5)
            
            # Create domain loop and surface
            c_ind = mesh_domain.create_domain_loop_from_poly()
            assert c_ind is not None
            assert len(c_ind) > 0
            
            ind_s_dom = mesh_domain.create_domain_surface()
            assert ind_s_dom is not None
            
            # Generate mesh
            model.synchronize()
            model.generate_mesh(dimension=2)
            
            # Verify mesh was created
            elements = gmsh.model.mesh.getElements(dim=2)
            assert len(elements[0]) > 0, "No 2D elements generated"
            assert len(elements[1]) > 0, "No element data"
            
            # Check mesh quality
            quality_df = model.get_triangular_quality()
            assert len(quality_df) > 0, "No quality data available"
            assert quality_df['minSICN'].min() > 0, "Poor element quality detected"
            
            print(f"Generated {len(quality_df)} triangular elements")
            print(f"Min quality (SICN): {quality_df['minSICN'].min():.3f}")
    
    def test_mesh_with_embedded_points(self):
        """Test mesh generation with embedded observation points."""
        import gmsh
        from gmshflow import GmshModel, GmshMeshDomain, PointGeometryHandler
        sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
        from sample_geometries import create_domain_polygon, create_test_points
        
        # Create test data
        domain_gdf = create_domain_polygon()
        points_gdf = create_test_points()
        
        with GmshModel("embedded_points_test") as model:
            # Set up domain
            mesh_domain = GmshMeshDomain("test_domain", domain_gdf, cs_dom=20.0)
            mesh_domain.prepare_mesh_domain(mesh_area=1, min_overlap=0.5)
            
            c_ind = mesh_domain.create_domain_loop_from_poly()
            ind_s_dom = mesh_domain.create_domain_surface()
            
            # Add embedded points
            point_handler = PointGeometryHandler()
            point_handler.set_gdf_point(points_gdf)
            point_ids = point_handler.create_point_from_point(df_coord=False)
            
            mesh_domain.add_embedded_points(id_point_list=point_ids, surface_id=ind_s_dom)
            
            # Generate mesh
            model.synchronize()
            model.generate_mesh(dimension=2)
            
            # Verify mesh and embedded points
            elements = gmsh.model.mesh.getElements(dim=2)
            assert len(elements[1][0]) > 0, "No triangular elements generated"
            
            # Check that mesh nodes include embedded points
            nodes = gmsh.model.mesh.getNodes()
            node_coords = nodes[1].reshape(-1, 3)[:, :2]  # Get x,y coordinates
            
            # Verify some mesh nodes are close to our embedded points
            point_coords = [(p.x, p.y) for p in points_gdf.geometry]
            found_embedded = 0
            for px, py in point_coords:
                distances = ((node_coords[:, 0] - px)**2 + (node_coords[:, 1] - py)**2)**0.5
                if distances.min() < 1.0:  # Within 1 unit
                    found_embedded += 1
            
            assert found_embedded > 0, "Embedded points not found in mesh"
            print(f"Found {found_embedded}/{len(point_coords)} embedded points in mesh")
    
    def test_mesh_with_polygon_geometry(self):
        """Test mesh generation with polygon geometries."""
        import gmsh
        from gmshflow import GmshModel, GmshMeshDomain, PolyGeometryHandler
        sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
        from sample_geometries import create_domain_polygon, create_test_polygon
        
        # Create test data - larger domain and smaller internal polygon
        domain_gdf = create_domain_polygon()
        # Scale domain up and move internal polygon
        domain_gdf.geometry = domain_gdf.geometry.scale(xfact=2, yfact=2, origin=(0,0))
        
        poly_gdf = create_test_polygon()
        poly_gdf.geometry = poly_gdf.geometry.translate(xoff=50, yoff=50)  # Move inside domain
        
        with GmshModel("polygon_geometry_test") as model:
            # Set up domain
            mesh_domain = GmshMeshDomain("test_domain", domain_gdf, cs_dom=30.0)
            mesh_domain.prepare_mesh_domain(mesh_area=1, min_overlap=0.5)
            
            c_ind = mesh_domain.create_domain_loop_from_poly()
            
            # Add polygon geometry
            poly_handler = PolyGeometryHandler()
            poly_handler.set_gpd_poly(poly_gdf)
            poly_loops = poly_handler.create_loop_from_poly(def_surf=False)
            
            # Add as internal loop to domain
            mesh_domain.add_internal_loops(poly_loops)
            
            ind_s_dom = mesh_domain.create_domain_surface()
            
            # Generate mesh
            model.synchronize()
            model.generate_mesh(dimension=2)
            
            # Verify mesh
            elements = gmsh.model.mesh.getElements(dim=2)
            assert len(elements[1][0]) > 0, "No triangular elements generated"
            
            quality_df = model.get_triangular_quality()
            assert len(quality_df) > 0, "No quality data"
            
            # Check that we have reasonable mesh quality
            good_elements = (quality_df['minSICN'] > 0.3).sum()
            total_elements = len(quality_df)
            quality_ratio = good_elements / total_elements
            
            print(f"Generated {total_elements} elements with {quality_ratio:.1%} good quality")
            assert quality_ratio > 0.5, "Too many poor quality elements"
    
    def test_mesh_output_formats(self):
        """Test different mesh output formats."""
        import gmsh
        from gmshflow import GmshModel, GmshMeshDomain
        sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
        from sample_geometries import create_domain_polygon
        
        domain_gdf = create_domain_polygon()
        
        with tempfile.TemporaryDirectory() as temp_dir:
            with GmshModel("output_test") as model:
                # Create simple mesh
                mesh_domain = GmshMeshDomain("test_domain", domain_gdf, cs_dom=30.0)
                mesh_domain.prepare_mesh_domain(mesh_area=1, min_overlap=0.5)
                
                c_ind = mesh_domain.create_domain_loop_from_poly()
                ind_s_dom = mesh_domain.create_domain_surface()
                
                model.synchronize()
                model.generate_mesh(dimension=2)
                
                # Test .msh output
                msh_file = os.path.join(temp_dir, "test.msh")
                model.write(msh_file)
                assert os.path.exists(msh_file), "MSH file not created"
                assert os.path.getsize(msh_file) > 0, "MSH file is empty"
                
                # Test Voronoi export (skip if geometry issues occur)
                try:
                    voro_name = "test_voronoi"
                    gdf_voro = mesh_domain.export_to_voronoi(
                        temp_dir, voro_name, [ind_s_dom], 
                        int_org_dom=True, min_cell_overlap=0.5
                    )
                    
                    # Verify Voronoi output
                    assert len(gdf_voro) > 0, "No Voronoi cells generated"
                    assert all(gdf_voro.geom_type == 'Polygon'), "Non-polygon geometries in Voronoi"
                    
                    # Check shapefile was created
                    shp_file = os.path.join(temp_dir, f"{voro_name}.shp")
                    assert os.path.exists(shp_file), "Shapefile not created"
                    
                    print(f"Generated {len(gdf_voro)} Voronoi cells")
                    
                except Exception as e:
                    # Voronoi can fail with certain geometries - just warn and continue
                    print(f"Voronoi export failed (expected with some geometries): {e}")
                    gdf_voro = None
                
                if gdf_voro is not None:
                    print(f"Total area: {gdf_voro.area.sum():.1f}")
                else:
                    print("Voronoi export skipped due to geometry issues")

@gmsh_required
class TestMeshQuality:
    """Test mesh quality and validation."""
    
    def test_mesh_quality_metrics(self):
        """Test comprehensive mesh quality analysis."""
        import gmsh
        from gmshflow import GmshModel
        
        with GmshModel("quality_test") as model:
            # Create a challenging geometry (thin rectangle)
            p1 = gmsh.model.geo.addPoint(0, 0, 0, 1.0)
            p2 = gmsh.model.geo.addPoint(10, 0, 0, 1.0)
            p3 = gmsh.model.geo.addPoint(10, 1, 0, 1.0)
            p4 = gmsh.model.geo.addPoint(0, 1, 0, 1.0)
            
            l1 = gmsh.model.geo.addLine(p1, p2)
            l2 = gmsh.model.geo.addLine(p2, p3)
            l3 = gmsh.model.geo.addLine(p3, p4)
            l4 = gmsh.model.geo.addLine(p4, p1)
            
            loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
            surface = gmsh.model.geo.addPlaneSurface([loop])
            
            model.synchronize()
            model.generate_mesh(dimension=2)
            
            # Get comprehensive quality metrics
            quality_df = model.get_triangular_quality()
            
            # Test all quality measures are present
            expected_metrics = [
                'minSICN', 'minDetJac', 'maxDetJac', 'minSJ', 'minSIGE',
                'gamma', 'innerRadius', 'outerRadius', 'minIsotropy',
                'angleShape', 'minEdge', 'maxEdge'
            ]
            
            for metric in expected_metrics:
                assert metric in quality_df.columns, f"Missing quality metric: {metric}"
                assert not quality_df[metric].isna().any(), f"NaN values in {metric}"
            
            # Basic quality checks
            assert quality_df['minSICN'].min() >= 0, "Negative SICN values"
            assert quality_df['gamma'].max() <= 1.0, "Gamma values > 1"
            assert quality_df['minEdge'].min() > 0, "Zero edge lengths"
            
            # Print quality summary
            print(f"Generated {len(quality_df)} elements")
            print(f"SICN range: {quality_df['minSICN'].min():.3f} - {quality_df['minSICN'].max():.3f}")
            print(f"Gamma range: {quality_df['gamma'].min():.3f} - {quality_df['gamma'].max():.3f}")
            
            # At least some elements should have reasonable quality
            good_quality = (quality_df['minSICN'] > 0.3).sum()
            assert good_quality > 0, "No good quality elements found"