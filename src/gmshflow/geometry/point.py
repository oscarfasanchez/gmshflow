"""Point geometry handler for GMSHFlow."""

import pandas as pd
import geopandas as gpd
import gmsh


class PointGeometryHandler:
    '''
    Class to handle point geometries for meshing in GMSH.
    
    Parameters
    ----------
    cs_point : float, optional
        Cell size of the point geometry. The default is None.
    '''
    def __init__(self, cs_point=None):
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
        
    def create_point_from_point(self, df_coord=False):
        '''
        This function creates a point from a point geometry.
        
        Parameters
        ----------
        p_ind : list
            List of the index of the points of the point geometry.
            '''
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