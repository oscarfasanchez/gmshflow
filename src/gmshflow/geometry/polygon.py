"""Polygon geometry handler for GMSHFlow."""

from typing import List, Optional, Tuple
import geopandas as gpd
import gmsh
from ..utils.preprocessing import simplify_keeping_topology


class PolyGeometryHandler:
    """Handler for polygon geometries in GMSH mesh generation.

    This class manages polygon geometries for mesh generation, including
    simplification, GMSH geometry creation, and specialized meshing strategies
    like buffer-based surface grids.

    Args:
        cs_poly: Default cell size for polygon geometries. If None, cell sizes
            must be specified in the 'cs' column of input GeoDataFrames.

    Attributes:
        gdf_poly: GeoDataFrame containing polygon geometries.
        cs_poly: Default cell size for polygons.
        gdf_coord: Coordinate points derived from polygon boundaries.
        s_ind: List of GMSH surface indices created from polygons.

    Example:
        >>> handler = PolyGeometryHandler(cs_poly=50.0)
        >>> handler.set_gpd_poly(polygon_gdf)
        >>> loops = handler.create_loop_from_poly(def_surf=True)
    """
    def __init__(self, cs_poly: Optional[float] = None) -> None:
        self.gdf_poly = None
        self.cs_poly = cs_poly
        self.gdf_coord = None
        #TODO include the following attributes inside self.gdf_poly
        self.s_ind = []

    def set_gpd_poly(self, gdf_poly: gpd.GeoDataFrame, keep_topology: bool = False) -> None:
        """Set the GeoDataFrame containing polygon geometries for meshing.

        Configures the polygon geometries and applies simplification to improve
        mesh quality. Validates input geometry types and cell size specifications.

        Args:
            gdf_poly: GeoDataFrame containing polygon geometries. Must have either
                a 'cs' column specifying cell sizes or cs_poly must be set.
            keep_topology: If True, uses topology-preserving simplification across
                all polygons. If False, simplifies each polygon independently.

        Raises:
            AssertionError: If geometries are not all Polygon type or if no cell
                size information is available.

        Example:
            >>> poly_gdf = gpd.read_file("barriers.shp")
            >>> handler.set_gpd_poly(poly_gdf, keep_topology=True)
        """
        # check it is a polygon dataframe
        assert all(gdf_poly.geom_type == 'Polygon'), 'All geometries must be of type Polygon'
        #check that it has a column called cs or cs_poly is not None
        assert 'cs' in gdf_poly.columns or self.cs_poly is not None, 'The geodataframe must have a cell size column or cs_poly must be defined'
        self.gdf_poly = gdf_poly
        #simplify the geometries
        if self.cs_poly is not None:
            if keep_topology:
                self.gdf_poly = simplify_keeping_topology(self.gdf_poly, self.cs_poly)
            else:    
                self.gdf_poly.geometry = self.gdf_poly.geometry.simplify(self.cs_poly/2)
            self.gdf_poly['cs'] = self.cs_poly
        else:
            if keep_topology:
                #get the biggest cell size to simplify the geometries
                print('Warning: the topology will be kept, but the geometries will be simplified'
                      'to the biggest cell size, even if some elements have smaller cell sizes'
                        ', this wont affect the final mesh size necessarily')
                self.cs_poly = self.cs_poly.max()
                self.gdf_poly = simplify_keeping_topology(self.gdf_poly, self.cs_poly)
            else:
                self.gdf_poly.geometry = self.gdf_poly.apply(lambda x: x.geometry.simplify(x.cs/2), axis=1)

    def create_loop_from_poly(self, def_surf=False):
        '''
        This function creates a loop from a polygon geometry.
        
        Parameters
        ----------
        def_surf : bool, optional
            If True, the function defines the surface of the polygon. The default is False.
        
        Returns
        -------
        c_ind : list
            List of the index of the curve loops of the polygon geometry.
            '''
        # This function creates a loop from a polygon geometry
        assert self.gdf_poly is not None, 'The polygon geometry is not defined'
        # lets get the points of the polygons, then the lines, and finally the area
        for i in self.gdf_poly.index:
            poly = self.gdf_poly.loc[i, 'geometry']
            poly_xy = list(zip(poly.exterior.xy[0], poly.exterior.xy[1]))
            p_ind = []
            for xy in poly_xy[:-1]:  # to avoid repeated points
                ind = gmsh.model.geo.addPoint(xy[0], xy[1], 0, self.gdf_poly.loc[i, 'cs'])# TODO check
                p_ind.append(ind)
            
            l_ind = []
            for j in range(len(p_ind)):
                ind = gmsh.model.geo.addLine(p_ind[j-1], p_ind[j])
                l_ind.append(ind)

            c_ind = []
            ind = gmsh.model.geo.addCurveLoop(l_ind)
            c_ind.append(ind)
               
            if def_surf:
                ind_s = gmsh.model.geo.addPlaneSurface(c_ind)
                self.s_ind.append(ind_s)
                self.gdf_poly.loc[i, 's_ind'] = ind_s
            # add the indices to a new column in the dataframe
            # TODO decide the final data storage system I will use
            self.gdf_poly.loc[i, 'c_ind'] = c_ind
            
            # self.gdf_poly.loc[i,'l_ind'] = l_ind
            # self.gdf_poly.loc[i,'p_ind'] = p_ind
        return self.gdf_poly['c_ind'].tolist()
        #TODO add a return?

    def create_surfacegrid_from_buffer_poly(self, cs_thick=1, simpl_fac=1.5, def_surf=True):
        '''
        This function creates a surface grid from a buffer polygon.
        this could be desired to have voronoi meshes that follow more accurately
        the original geometry(because normally voronoi elements will have their
        centers on the original gometry instead of their boundaries),
        but it is not recommended for complex geometries.

        Parameters
        ----------
        cs_thick : int, optional
            Offset of the buffer polygon. The default is 1.
        simpl_fac : float, optional
            Simplification factor for the buffer polygon. The default is 1.5.
        def_surf : bool, optional   
            If True, the function defines the surface of the buffer polygon. The default is True.
        
        Returns
        -------
        gdf_poly : geopandas.GeoDataFrame
            Geodataframe with the buffer polygon geometries.
        c_ind_buf_pos : list
            List of the index of the curve loops of the positive buffer polygon.
        c_ind_buf_neg : list
            List of the index of the curve loops of the negative buffer polygon.
        ind_s_buff : list
            List of the index of the buffer polygon surface.
        ind_s_mid : list
            List of the index of the buffer polygon mid surface.

        '''
        # TODO check changes in buffer to avoid small angles,
        # and probably different meshing algorithm?
        print('starting: create_surfacegrid_from_buffer_poly')
        assert cs_thick >= 1 and cs_thick <= 2, 'this function was designed only for offset of 1 or 2 cells, if used for more it results in bad quality meshes'
        c_ind_buf_pos = []
        c_ind_buf_neg = []
        ind_s_buff = []
        ind_s_mid = []
        #the lack of space to mesh requires further simplification
        self.gdf_poly.geometry = self.gdf_poly.geometry.simplify(self.cs_poly * simpl_fac)
        for i in self.gdf_poly.index:
            off_pos = self.gdf_poly.loc[i, "geometry"].buffer(cs_thick * self.cs_poly / 2, quad_segs=1, join_style=2, mitre_limit=5.0).boundary
            off_neg = self.gdf_poly.loc[i, "geometry"].buffer(-cs_thick * self.cs_poly / 2, quad_segs=1, join_style=2, mitre_limit=5.0).boundary

            # get nodes of buffer to assign them to the mesh iteratively
            off_pos_xy = list(zip(off_pos.xy[0], off_pos.xy[1]))
            off_neg_xy = list(zip(off_neg.xy[0], off_neg.xy[1]))

            p_ind_pos = []
            for xy in off_pos_xy[:-1]:  # to avoid repeated points
                ind = gmsh.model.geo.addPoint(xy[0], xy[1], 0, self.cs_poly)
                p_ind_pos.append(ind)

            p_ind_neg = []
            for xy in off_neg_xy[:-1]:  # to avoid repeated points
                ind = gmsh.model.geo.addPoint(xy[0], xy[1], 0, self.cs_poly)
                p_ind_neg.append(ind)

            # define lines
            l_ind_pos = []
            for j in range(len(p_ind_pos)):
                ind = gmsh.model.geo.addLine(p_ind_pos[j-1], p_ind_pos[j])
                l_ind_pos.append(ind)

            l_ind_neg = []
            for i in range(len(p_ind_neg)):
                ind = gmsh.model.geo.addLine(p_ind_neg[i-1], p_ind_neg[i])
                l_ind_neg.append(ind)

            # define loops
            loop_list_neg = l_ind_neg
            loop_list_pos = l_ind_pos

            c_ind_buf_neg.append(gmsh.model.geo.addCurveLoop(loop_list_neg))
            c_ind_buf_pos.append(gmsh.model.geo.addCurveLoop(loop_list_pos))
            self.c_ind.append([c_ind_buf_neg, c_ind_buf_pos])
            print('Loop list created')
            if def_surf:
                # define surfaces
                ind_s_buff.append(gmsh.model.geo.addPlaneSurface([c_ind_buf_pos[-1], c_ind_buf_neg[-1]]))
                ind_s_mid.append(gmsh.model.geo.addPlaneSurface([c_ind_buf_neg[-1]]))
                self.s_ind.append([ind_s_mid, ind_s_buff])

                gmsh.model.geo.synchronize()
                # set quads for buffer area
                gmsh.model.mesh.setRecombine(2, ind_s_buff[-1])
                gmsh.model.mesh.setAlgorithm(2, ind_s_buff[-1], 8)
                gmsh.model.mesh.setAlgorithm(2, ind_s_mid[-1], 6)
                self.gdf_poly.loc[self.gdf_poly.shape[0], ['layer', 'geometry']] = [f'off_pos{i}', off_pos]
                self.gdf_poly.loc[self.gdf_poly.shape[0], ['layer', 'geometry']] = [f'off_neg{i}', off_neg]
        print('finished: create_surfacegrid_from_buffer_poly')
        return self.gdf_poly, c_ind_buf_pos, c_ind_buf_neg, ind_s_buff, ind_s_mid

    def convert_to_points_for_size_fields(self):
        '''
        This function converts the polygon geometry to points for cell size fields.
        
        Returns
        -------
        gdf_coord : geopandas.GeoDataFrame
            Geodataframe with the point geometry.
        '''
        # Implementation of conversion to points for threshold fields
        assert 'cs' in self.gdf_poly.columns, 'gdf must have cell size column(cs)'
        gdf2 = self.gdf_poly.copy()
        # apply to have different point densities in different segments of drn
        gdf2['geometry'] = self.gdf_poly.apply(lambda x: x.geometry.simplify(x.cs), axis=1)  # TODO cs_bar/2??
        gdf2['geometry'] = gdf2.apply(lambda x: x.geometry.segmentize(x.cs), axis=1)
        # get xy values in a new dataframe
        gdf_coord = gdf2.get_coordinates().reset_index(names='ind')
        gdf_coord['cs'] = gdf_coord.apply(  # get cs from previous gdf
            lambda x: gdf2.loc[x.ind, 'cs'], axis=1)  # assign cell size
        # to filter repeated values
        idx = gdf_coord.groupby(['x', 'y'])['cs'].idxmin()
        gdf_coord = gdf_coord.loc[idx]
        # add points to gmsh
        gdf_coord['id_gmsh'] = gdf_coord.apply(
            lambda x: gmsh.model.geo.addPoint(x.x, x.y, 0, x.cs), axis=1)
        self.gf_coord = gdf_coord
        return gdf_coord