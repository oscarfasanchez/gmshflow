"""Point geometry handler for GMSHFlow."""

from typing import List, Optional
import pandas as pd
import geopandas as gpd
import gmsh


class PointGeometryHandler:
    """Handler for point geometries in GMSH mesh generation.

    Manages Point geometries for mesh generation, primarily used for
    mesh size control and embedding constraint points in the mesh.

    Args:
        cs_point: Default cell size for point geometries. If None, cell sizes
            must be specified in the 'cs' column of input GeoDataFrames.

    Attributes:
        gdf_point: GeoDataFrame containing point geometries.
        cs_point: Default cell size for points.
        gdf_coord: Coordinate information for points.

    Example:
        >>> handler = PointGeometryHandler(cs_point=10.0)
        >>> handler.set_gdf_point(well_locations_gdf)
        >>> point_indices = handler.create_point_from_point()
    """
    def __init__(self, cs_point: Optional[float] = None) -> None:
        self.gdf_point = None
        self.cs_point = cs_point
        self.gdf_coord = None

    def set_gdf_point(self, gdf_point: gpd.GeoDataFrame):
        '''
        This function sets the geodataframe point for meshing.'''
        # check it is a point dataframe
        assert all(gdf_point.geom_type == 'Point'), 'All geometries must be of type Point'
        #check that it has a column called cs or cs_point is not None
        assert 'cs' in gdf_point.columns or self.cs_point is not None, 'The geodataframe must have a cell size column or cs_point must be defined'
        self.gdf_point = gdf_point
        #simplify the geometries
        if self.cs_point is not None:
            self.gdf_point['cs'] = self.cs_point
        
    def create_point_from_point(self, df_coord: bool = False) -> List[int]:
        """Create GMSH points from point geometries.

        Converts point geometries to GMSH point entities that can be used
        for mesh size control or embedding in the mesh.

        Args:
            df_coord: If True, creates coordinate DataFrame with GMSH IDs.

        Returns:
            List of GMSH point indices created.

        Example:
            >>> point_ids = handler.create_point_from_point(df_coord=True)
            >>> print(f"Created {len(point_ids)} GMSH points")
        """
        # This function creates a gmsh point from a point geometry
        assert self.gdf_point is not None, 'The point geometry is not defined'
        # lets get the points of the polygons, then the lines, and finally the area
        p_ind = []
        for i in self.gdf_point.index:
            point = self.gdf_point.loc[i]
            ind = gmsh.model.geo.addPoint(point.geometry.x, point.geometry.y, 0, self.gdf_point.loc[i, 'cs'])# TODO check
            self.gdf_point.loc[i, 'p_ind'] = ind
            p_ind.append(ind)
        if df_coord:
            self.gdf_coord = self.gdf_point
            self.gdf_coord = pd.concat(
                [self.gdf_coord, self.gdf_point.get_coordinates()], axis=1) 
            self.gdf_coord['id_gmsh'] = self.gdf_coord['p_ind'].astype(int)
        return p_ind