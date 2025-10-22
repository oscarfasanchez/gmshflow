"""GMSHFlow: A wrapper for GMSH to create meshes for MODFLOW6 groundwater flow models.

This package provides an object-oriented interface to convert GIS data (shapefiles
in GeoPandas format) into triangular/Voronoi meshes suitable for groundwater modeling.
"""

from .core.domain import GmshMeshDomain
from .core.model import GmshModel
from .geometry.line import LineGeometryHandler
from .geometry.point import PointGeometryHandler
from .geometry.polygon import PolyGeometryHandler
from .utils.preprocessing import (
    merge_many_multilinestring_into_one_linestring,
    simplify_keeping_topology,
)

__version__ = "0.1.0"
__author__ = "GMSHFlow Contributors"

__all__ = [
    "GmshModel",
    "GmshMeshDomain",
    "PolyGeometryHandler",
    "LineGeometryHandler",
    "PointGeometryHandler",
    "merge_many_multilinestring_into_one_linestring",
    "simplify_keeping_topology"
]
