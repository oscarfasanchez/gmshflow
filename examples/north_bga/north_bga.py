#import libraries for doing the mesh
import geopandas as gpd
import pandas as pd
import os
# import gmshflow from the path: ./01_Python/01_python
import sys
sys.path.append(os.path.abspath(os.path.join('..','..' ,'01_python')))
import gmshflow

# set the path to the shapefiles
wdshp = os.path.join('.', 'data')

# Load the shapefiles

#load domain
dom = gpd.read_file(os.path.join(wdshp, 'Domin_Mod.shp'))
#load rivers
riv = gpd.read_file(os.path.join(wdshp, 'Rios.shp'))
#load creeks
crks = gpd.read_file(os.path.join(wdshp, 'queb_corr.shp'))
#load observations and filter pizometer observations with the column DEPTH_MEA > 0
obs = gpd.read_file(os.path.join(wdshp, 'INV_PAS_V5_DEM.shp'))
obs = obs[obs['DEPTH_MEA'] > 0]
#load Faults
faults = gpd.read_file(os.path.join(wdshp, 'Fallas_V2_EditPau.shp'))
#load Gallery
gal = gpd.read_file(os.path.join(wdshp, 'Alineamiento_Galeria.shp'))
#load barrier
bar = gpd.read_file(os.path.join(wdshp, 'barrier.shp'))

# Create the mesh

# Define cell_size in each area
cs_dom = 250.0
cs_riv = 250.0
cs_crk = 150.0
cs_obs = 10.0
cs_fault = 80.0
cs_barrier = 10
cs_gal = 3

# Add cell size column
riv['cs'] = cs_riv
crks['cs'] = cs_crk
obs['cs'] = cs_obs
faults['cs'] = cs_fault
gal['cs'] = cs_gal
bar['cs'] = cs_barrier

# set cs in crks=20 if its less than 500 meters from gal
crks.loc[crks.distance(gal.geometry[0]) < 300, 'cs'] = 50

# Initialize GmshModel and GmshMeshDomain
gmsh_model = gmshflow.GmshModel("north_bga_gmsh")
mesh_domain = gmshflow.GmshMeshDomain("domain", dom, cs_dom)

# Prepare mesh domain
mesh_domain.prepare_mesh_domain(mesh_area=1, min_overlap=0.5)

# Create domain loop
c_ind = mesh_domain.create_domain_loop_from_poly()

# Add river polygon geometry
riv_handler = gmshflow.PolyGeometryHandler()
riv_handler.set_gpd_poly(riv)


# Add creek line geometry
crk_handler = gmshflow.LineGeometryHandler()
crk_handler.set_gpd_line(crks)
crk_handler.create_line_from_line()

#add barrier line geometry

bar_handler = gmshflow.LineGeometryHandler()
bar_handler.set_gpd_line(bar)
bar_handler.create_surfacegrid_from_buffer_line(cs_thick=1)

# Add fault line geometry
fault_handler = gmshflow.LineGeometryHandler()
fault_handler.set_gpd_line(faults)
fault_handler.create_line_from_line()

# Add gallery line geometry, and create a surface
gal_handler = gmshflow.LineGeometryHandler()
gal_handler.set_gpd_line(gal)
gal_handler.create_surfacegrid_from_buffer_line(cs_thick=2)

# # add observation points
# obs_handler = gmshflow.PointGeometryHandler()
# obs_handler.set_gdf_point(obs)
# obs_handler.create_point_from_point(df_coord=True)
# gdf_obs_coord = obs_handler.gdf_coord

# # Convert to points for threshold fields
gdf_riv_coord = riv_handler.convert_to_points_for_size_fields()
gdf_crk_coord = crk_handler.convert_to_points_for_size_fields()
gdf_fault_coord = fault_handler.convert_to_points_for_size_fields()
gdf_gal_coord = gal_handler.convert_to_points_for_size_fields()


# # Create exponential field
# gdf_coords = pd.concat([gdf_obs_coord, gdf_riv_coord, gdf_crk_coord, gdf_fault_coord, gdf_gal_coord])
df_coords = pd.concat([gdf_riv_coord, gdf_gal_coord, gdf_crk_coord, gdf_fault_coord])
# mesh_domain.create_exponential_field(gdf_coords, fac=1.1)
mesh_domain.create_exponential_field(df_coords, fac=1.1)
# # Set exponential field as background mesh
mesh_domain.set_field()

# # Deactivate the mesh size from the geometry points
mesh_domain.set_mesh_size_from_geometries(use_boundaries=True, use_points=False)

# Create domain surface

#need to add internal loops,then we can create the domain surface
mesh_domain.add_internal_loops(gal_handler.c_ind_buf+bar_handler.c_ind_buf)

ind_s_dom = mesh_domain.create_domain_surface()

# Generate mesh
print('Generating the mesh...')
gmsh_model.generate_mesh()
print('Mesh ready')

#get some quality stats of the triangular mesh
df_quality = gmsh_model.get_triangular_quality()
#get gamma( ratio of the inscribed to circumcribed sphere/circle radius)
df_quality.gamma.plot.hist(density=True)
print('min gamma is equal to: ', df_quality.gamma.min())

# Run GUI 
gmsh_model.run_gui()

# Export to Voronoi
surf_tags = [ind_s_dom]+gal_handler.ind_s_buff  # Add other surface tags if needed
voro = mesh_domain.export_to_voronoi(wdshp, "gdf_voro", surf_tags)

# Write mesh to file
gmsh_model.write("north_bga_gmsh.msh")

# Finalize Gmsh
gmsh_model.finalize()









