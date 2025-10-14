# gmshflow

`gmshflow` is a Python wrapper for GMSH that allows users to create Voronoi meshes using GeoPandas and GMSH. It is designed to facilitate the creation of meshes for MODFLOW6 groundwater flow models.

## Features

- Create meshes for MODFLOW6 models using shapefiles or other formats compatible with  GeoPandas.
- Generate Voronoi meshes.
- Simplify and buffer geometries for better mesh quality.
- Support for various geometric features such as domains, rivers, barriers, faults, and pits.
- Export meshes to shapefiles and MODFLOW6-compatible formats.

## Generalities
This wrapper inherits the basic structure of the GMSH API.
The most important principles to take into account for this wrapper are:
Entities are built in a bottom-up manner using its Boundary Representation (BRep): 
first points, then lines(made of points), curves/loops (made of lines) and surfaces(made of curves)
Where each entity possess a unique tag

by default the mesh size is defined from points and from boundaries.
The cell sizes are assigned on points, the points inherit the cell size to the lines, and so on.
Then the mesh inherits the cell size from the points that were used to define the boundaries(loop) of the surfaces.

Other entities can print their cell size directly into the surface mesh by explicitly embedding the entity into
the surface, this will force the mesh to the shape of the entity and will fix the maximum cell size.
The definition of the cell size in the surfaces comes from the interpolation between the aforementioned features.

The interpolation can cause over-refinement, so it's recommended to deactivate the interpolation and use fields with the target mesh size.
This package has functions to create those fields easily.
In any case unless explicitly stated, the maximum cell size, will be calculated using the minimum values set by all the previous methods.




## Installation

Recommended on Windows: use Conda for geospatial dependencies, then perform an editable install.

1) Create/activate the Conda environment (once):

```powershell
# From the repo root
conda env create -f src/env.yml
conda activate gmshflow
```

2) Editable install so `import gmshflow` works from anywhere in this env:

```powershell
# From the repo root
pip install -e .
```

If `conda` or `python` are not recognized in PowerShell, use the "Anaconda Prompt" (installed by Anaconda/Miniconda) or ensure they are on your PATH. You can also run inside VS Code after selecting the Conda interpreter.

Notes:
- Heavy deps (geopandas, shapely, pyproj, etc.) are managed by Conda via `src/env.yml`.
- `pyproject.toml` intentionally omits dependencies to avoid pip trying to build geospatial wheels on Windows.
- After the editable install, remove any `sys.path.append(...)` hacks from examples.

## Usage

### Example

Here is an example of how to use `gmshflow` to create a mesh with a domain and observation points:

```python
import os
import geopandas as gpd
from gmshflow import GmshModel, GmshMeshDomain, PointGeometryHandler

# Set the path to the shapefiles
wdshp = os.path.join('.', 'data')

# Load the shapefiles
# Load domain
dom = gpd.read_file(os.path.join(wdshp, 'Domin_Mod.shp'))
# Load observations and filter pizometer observations with the column DEPTH_MEA > 0
obs = gpd.read_file(os.path.join(wdshp, 'INV_PAS_V5_DEM.shp'))
obs = obs[obs['DEPTH_MEA'] > 0]
# Get a geoseries with the unique points
obs_geo = obs.remove_repeated_points(tolerance=1.0).normalize().drop_duplicates()
# Make a geodataframe with the unique points
obs = gpd.GeoDataFrame(geometry=obs_geo)

# Define cell_size in each area
cs_dom = 250.0
cs_obs = 30.0
# Add cell size column
obs['cs'] = cs_obs

# Initialize GmshModel and GmshMeshDomain
gmsh_model = GmshModel("north_bga_gmsh")
mesh_domain = GmshMeshDomain("domain", dom, cs_dom)

# Prepare mesh domain
mesh_domain.prepare_mesh_domain(mesh_area=1, min_overlap=0.5)

# Create domain loop
c_ind = mesh_domain.create_domain_loop_from_poly()

# Add observation points
obs_handler = PointGeometryHandler()
obs_handler.set_gdf_point(obs)
obs_ind_list = obs_handler.create_point_from_point(df_coord=False)

# Create domain surface
ind_s_dom = mesh_domain.create_domain_surface()

# Embed observation points into the surface
mesh_domain.add_embedded_points(id_point_list=obs_ind_list, surface_id=ind_s_dom)

# Generate mesh
print('Generating the mesh...')
gmsh_model.generate_mesh()
print('Mesh ready')

# Export to Voronoi
surf_tags = [ind_s_dom]
shp_mesh_name = "gdf_voro_simple"
gdf_voro = mesh_domain.export_to_voronoi(wdshp, shp_mesh_name, surf_tags)

# Finalize Gmsh
gmsh_model.finalize()
```

### Quick verification

After installation, you can quickly verify the module is importable:

```powershell
python -c "import gmshflow; print(getattr(gmshflow, '__name__', 'ok'))"
```


## License

This project is licensed under the MIT License.
