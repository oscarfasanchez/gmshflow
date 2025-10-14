#import libraries for doing the mesh
import geopandas as gpd
import numpy as np
import os
import gmshflow
import flopy
from flopy.discretization import VertexGrid
from pathlib import Path

# set the path to the shapefiles
# set the path to the shapefiles relative to this script
SCRIPT_DIR = Path(__file__).resolve().parent
wdshp = SCRIPT_DIR / "data"
if not wdshp.exists():
    raise FileNotFoundError(f"Data folder not found: {wdshp}")

# Load the shapefiles
#load domain
dom = gpd.read_file(os.path.join(wdshp, 'Domin_Mod.shp'))
#load observations and filter pizometer observations with the column DEPTH_MEA > 0
obs = gpd.read_file(os.path.join(wdshp, 'INV_PAS_V5_DEM.shp'))
obs = obs[obs['DEPTH_MEA'] > 0]
#get a geoseries with the unique points
obs_geo = obs.remove_repeated_points(tolerance=1.0).normalize().drop_duplicates()
# make a geodataframe with the unique points
obs = gpd.GeoDataFrame(geometry=obs_geo)

# Create the mesh

# Define cell_size in each area
cs_dom = 250.0
cs_obs = 30.0
# Add cell size column
obs['cs'] = cs_obs


# Initialize GmshModel and GmshMeshDomain
gmsh_model = gmshflow.GmshModel("north_bga_gmsh")
mesh_domain = gmshflow.GmshMeshDomain("domain", dom, cs_dom)

# Prepare mesh domain
mesh_domain.prepare_mesh_domain(mesh_area=1, min_overlap=0.5)

# Create domain loop
c_ind = mesh_domain.create_domain_loop_from_poly()

# add observation points
obs_handler = gmshflow.PointGeometryHandler()
obs_handler.set_gdf_point(obs)
obs_ind_list = obs_handler.create_point_from_point(df_coord=False)

# Create domain surface
ind_s_dom = mesh_domain.create_domain_surface()

#now we can embed elements into the surface
mesh_domain.add_embedded_points(id_point_list=obs_ind_list, surface_id=ind_s_dom)

# Generate mesh
print('Generating the mesh...')
gmsh_model.generate_mesh()
print('Mesh ready')

# Export to Voronoi
surf_tags = [ind_s_dom]
shp_mesh_name = "gdf_voro_simple"
gdf_voro = mesh_domain.export_to_voronoi(wdshp, shp_mesh_name, surf_tags)

#savefig of the plot of geodataframe of the mesh
ax = gdf_voro.plot()
fig = ax.get_figure()
fig.savefig(os.path.join(wdshp, f'{shp_mesh_name}.png'))

# Finalize Gmsh
gmsh_model.finalize()

#convert shapefile to cvfd to be imported in modflow6
verts, iverts = flopy.utils.cvfdutil.shapefile_to_cvfd(
                os.path.join(wdshp, f'{shp_mesh_name}.shp'), verbose=True, duplicate_decimals=3, skip_hanging_node_check=True)
gridprops = flopy.utils.cvfdutil.get_disv_gridprops(verts, iverts, xcyc=None)

vertices = gridprops["vertices"]
cell2d = gridprops["cell2d"]
nlay = 1
ncpl = gridprops["ncpl"]
idomain = np.ones((nlay, ncpl), dtype=int)

modelgrid = VertexGrid(
    vertices=vertices,
    cell2d=cell2d,
    idomain=idomain,
    nlay=nlay,
    ncpl=ncpl,
)











