# gmshflow

`gmshflow` is a Python wrapper for GMSH that allows users to create Voronoi meshes using GeoPandas and GMSH. It is designed to facilitate the creation of meshes for MODFLOW6 groundwater flow models.

## Features

- Create meshes for MODFLOW6 models using shapefiles in GeoPandas format.
- Generate Voronoi and triangular meshes.
- Simplify and buffer geometries for better mesh quality.
- Support for various geometric features such as domains, rivers, barriers, faults, and pits.
- Export meshes to shapefiles and MODFLOW6-compatible formats.

## Generalities
This wrapper inherits the basic structure of GMSH API.
The most important principles to take into account for this wrapper are:
Entities are built in a bottom-up manner using its Boundary Representation (BRep): 
first points, then lines(made of points), curves/loops (made of lines) and surfaces(made of curves)
Where each entity possess a unique tag

by default the mesh size is defined from points and from boundaries
The cell sizes are assigned on points, the points inherit the cell size to the lines, and so on.
Then the mesh inherits the cell size from the points that were used to define the boundaries(loop) of the surfaces.

Other entities can print their cell size directly into the surface mesh by explicitly embedding the entity into
the surface, this will force the mesh to the shape of the entity and will fix the maximum cell size.
The definition of the cell size in the surfaces comes from the interpolation between the aforementioned features.

The interpolation can cause over refinement, then its recommended to use fields with the target mesh size, this package
has functions to create those fields easily. In any case unless explicitly stated, the maximum cell size, will be calculated using the minimum values set by all the previous methods.




## Installation

To install the required packages, use the following command:

```bash
pip install ...???...
```
nah, Probably would be more like copy paste, we'll see

## Usage

### Example

Here is an example of how to use `gmshflow` to create a mesh:

```python
import os
import geopandas as gpd
from gmshflow import GmshModel, GmshMeshDomain, PolyGeometryHandler, LineGeometryHandler

# Import all shapefiles
print('hello world, some documentation in coming')
```

## License

This project is licensed under the XXX License.
