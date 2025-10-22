"""Line geometry handler for GMSHFlow."""

from statistics import mean
from typing import Optional

import geopandas as gpd
import numpy as np
from shapely.geometry import Point

try:
    import gmsh
    HAS_GMSH = True
except ImportError:
    gmsh = None
    HAS_GMSH = False

from ..utils.preprocessing import (
    merge_many_multilinestring_into_one_linestring,
    simplify_keeping_topology,
)


class LineGeometryHandler:
    """Handler for line geometries in GMSH mesh generation.

    Manages LineString and MultiLineString geometries for mesh generation,
    including creation of lines, buffer-based surface grids, and conversion
    to point fields for mesh size control.

    Args:
        cs_line: Default cell size for line geometries. If None, cell sizes
            must be specified in the 'cs' column of input GeoDataFrames.

    Attributes:
        gdf_line: GeoDataFrame containing line geometries.
        cs_line: Default cell size for lines.
        gdf_coord: Coordinate points derived from line geometries.
        ind_s_buff: List of buffer surface indices.
        c_ind_buf: List of buffer curve loop indices.
        l_ind_list: List of GMSH line indices.

    Example:
        >>> handler = LineGeometryHandler(cs_line=25.0)
        >>> handler.set_gpd_line(river_lines_gdf)
        >>> lines = handler.create_line_from_line()
    """
    def __init__(self, cs_line: Optional[float] = None) -> None:
        self.gdf_line = None
        self.cs_line = cs_line
        self.gdf_coord = None
        self.ind_s_buff = []
        self.c_ind_buf = []
        self.l_ind_list = []

    def set_gpd_line(self, gdf_line: gpd.GeoDataFrame, keep_topology=False):
        '''
        This function sets the geodataframe line for meshing.

        Parameters
        ----------
        gdf_line : geopandas.GeoDataFrame
            Geodataframe with the line geometry.
        keep_topology : bool, optional
            If True, the topology of the line geometry is kept, but individual cell sizes wont we used
            for the simplifications. The default is False.

        '''
        # check it is a line dataframe
        valid_types = ['LineString', 'MultiLineString']
        if not all(gdf_line.geom_type.isin(valid_types)):
            invalid_types = gdf_line.geom_type[~gdf_line.geom_type.isin(valid_types)].unique()
            raise ValueError(f'All geometries must be LineString or MultiLineString. Found: {invalid_types}')

        # check that it has a column called cs or cs_line is not None
        if 'cs' not in gdf_line.columns and self.cs_line is None:
            available_cols = list(gdf_line.columns)
            raise ValueError(
                f"Either 'cs' column must exist in GeoDataFrame or cs_line must be set. "
                f"Available columns: {available_cols}"
            )

        # TODO simplify keeping topology through the package topojson
        print('Warning: the line geometries will be simplified without keeping topology, check the results')
        if any(gdf_line.geom_type == 'MultiLineString'):
            self.gdf_line = gdf_line.explode().reset_index()
        else:
            self.gdf_line = gdf_line
        #simplify the geometries
        if self.cs_line is not None:
            if keep_topology:
                self.gdf_line = simplify_keeping_topology(self.gdf_line, self.cs_line)
            else:
                self.gdf_line.geometry = self.gdf_line.geometry.simplify(self.cs_line/2)
            self.gdf_line['cs'] = self.cs_line
        else:
            if keep_topology:
                #get the biggest cell size to simplify the geometries
                print('Warning: the topology will be kept, but the geometries will be simplified'
                      ' to the biggest cell size, even if some elements have smaller cell sizes'
                        ', this wont affect the final mesh size necessarily')
                self.cs_line = self.gdf_line['cs'].max()
                self.gdf_line = simplify_keeping_topology(self.gdf_line, self.cs_line)
            else:
                self.gdf_line.geometry = self.gdf_line.apply(lambda x: x.geometry.simplify(x.cs/2), axis=1)

    def create_line_from_line(self):
        '''
        This function creates a line from a line geometry.

        Returns
        -------
        l_ind_list : list
            List of the indices of the lines of the line geometry.'''
        # This function creates a gmsh point and line from a line geometry
        if self.gdf_line is None:
            raise RuntimeError('Line geometry is not defined. Call set_gpd_line() first.')
        # lets get the points of the polygons, then the lines, and finally the area
        for i in self.gdf_line.index:
            line = self.gdf_line.loc[i, 'geometry']
            line_xy = list(zip(line.xy[0], line.xy[1]))
            p_ind = []
            for xy in line_xy[:]:  # to avoid repeated points
                ind = gmsh.model.geo.addPoint(xy[0], xy[1], 0, self.gdf_line.loc[i, 'cs'])# TODO check
                p_ind.append(ind)

            l_ind = []
            for j in range(len(p_ind)-1):
                ind = gmsh.model.geo.addLine(p_ind[j], p_ind[j+1])
                l_ind.append(ind)
            self.l_ind_list += l_ind
        # self.gdf _line['p_ind'] = p_ind

        return self.l_ind_list

    def create_surfacegrid_from_buffer_line(self, cs_thick: int = 1, simpl_fac=1.5, def_surf=True):
        '''
        This function creates a surface grid from a buffer line.
        this could be used to create a clean sharp boundary on the voronoi mesh with the line shape
        (cs_thick = 1) or to create a quasi-rectangular mesh around the line shape (cs_thick = 2)
        Parameters
        ----------
        cs_thick : int, optional
            Offset of the buffer line. The default is 1.
        simpl_fac : float, optional
            Simplification factor for the buffer line. The default is 1.5.
        def_surf : bool, optional
            If True, the function defines the surface of the buffer line. The default is True.
        Internal Returns
        -------
        gdf_line : geopandas.GeoDataFrame
            Geodataframe with the buffer line geometries.
        ind_s_buff : list
            List of the index of the buffer line surface.

        Returns
        -------
        c_ind_buf : list
            List of the index of the curve loops of the buffer line.

        '''
        # lets simplify the barrier to ease the quad meshing
        # Failed attempt to smooth angles, better just enforce 1 or 2 to avoid
        # misuse of it with poor meshes
        print('starting: create_surfacegrid_from_buffer_line')
        print('Warning, this function expects that different feature lines are not intersecting')
        # to get decent quality mesh, the offset should be at least 1 cell and considerable simplification
        if not (1 <= cs_thick <= 2):
            raise ValueError(
                f'Buffer thickness must be 1 or 2 cells. Got: {cs_thick}. '
                'Higher values result in poor mesh quality.'
            )
        ind_s_buff = []

        #the lack of space to mesh requires further simplification
        self.gdf_line.geometry = self.gdf_line.apply(lambda x: x.geometry.simplify(x.cs * simpl_fac), axis=1)
        #  two ways of segmentize, in the original line, or at the buffers,
        # is better in the original to make it symmetrical
        self.gdf_line.geometry = self.gdf_line.apply(lambda x: x.geometry.segmentize(x.cs), axis=1)
        for i in self.gdf_line.index:
            cs_line = self.gdf_line.loc[i, 'cs']
            off_pos = self.gdf_line.loc[i, "geometry"].offset_curve(
                cs_thick * cs_line / 2, quad_segs=1, join_style=2, mitre_limit=5.0)
            off_neg = self.gdf_line.loc[i, "geometry"].offset_curve(
                -cs_thick * cs_line / 2, quad_segs=1, join_style=2, mitre_limit=5.0)

            #save the offset geometry on the original line gdf
            self.gdf_line.loc[self.gdf_line.shape[0], ['layer', 'geometry', 'cs']] = [f'off_pos_{i}', off_pos, self.gdf_line.loc[i, 'cs']]
            self.gdf_line.loc[self.gdf_line.shape[0], ['layer', 'geometry', 'cs']] = [f'off_neg_{i}', off_neg, self.gdf_line.loc[i, 'cs']]
            if any((self.gdf_line.geom_type == 'MultiLineString') & (self.gdf_line.layer.str.contains(f'_{i}'))):
                # workaround for geos bug that creates multilinestring sometimes
                # Make the multilinestring a linesting
                # this split it in two linestrings, but how to merge it into one linestring?
                # self.gdf_line.loc[self.gdf_line.geom_type=='MultiLineString'].explode()
                self.gdf_line = merge_many_multilinestring_into_one_linestring(self.gdf_line)
                # reassign the merged linestring to the original line
                off_pos = self.gdf_line.loc[self.gdf_line.layer == f'off_pos_{i}', 'geometry'].values[0]
                off_neg = self.gdf_line.loc[self.gdf_line.layer == f'off_neg_{i}', 'geometry'].values[0]


            # get nodes of buffer to assign them to the mesh iteratively
            off_pos_xy = list(zip(off_pos.xy[0], off_pos.xy[1]))
            off_neg_xy = list(zip(off_neg.xy[0], off_neg.xy[1]))

            p_ind_pos = []
            for xy in off_pos_xy:  # to avoid repeated points?
                ind = gmsh.model.geo.addPoint(xy[0], xy[1], 0, cs_line)
                p_ind_pos.append(ind)

            p_ind_neg = []
            for xy in off_neg_xy:
                ind = gmsh.model.geo.addPoint(xy[0], xy[1], 0, cs_line)
                p_ind_neg.append(ind)
            # define lines
            l_ind_pos = []
            for j in range(len(p_ind_pos)-1):
                ind = gmsh.model.geo.addLine(p_ind_pos[j], p_ind_pos[j+1])
                l_ind_pos.append(ind)

            l_ind_neg = []
            for j in range(len(p_ind_neg)-1):
                ind = gmsh.model.geo.addLine(p_ind_neg[j], p_ind_neg[j+1])
                l_ind_neg.append(ind)

            l_strt = gmsh.model.geo.addLine(p_ind_pos[0], p_ind_neg[0])
            l_end = gmsh.model.geo.addLine(p_ind_neg[-1], p_ind_pos[-1])
            # define loops, carefully, in the right order
            loop_list = l_ind_neg+ [l_end]+ list(np.array(l_ind_pos[::-1])*-1)+ [l_strt]
            c_ind_buf = gmsh.model.geo.addCurveLoop(loop_list)
            # define surfaces
            ind_s_buff = gmsh.model.geo.addPlaneSurface([c_ind_buf])
            # TODO check if I should restart the variable when the function is called, or create another function to clean
            self.ind_s_buff.append(ind_s_buff)# should I standarize the names using inheritance? s_ind?
            self.c_ind_buf.append(c_ind_buf)

            #to improve the quad mesh quality, the  tranfinite curves are set
            # The side of the rectanule perpendicular to the original line is divided according to the buffer size
            gmsh.model.geo.mesh.setTransfiniteCurve(l_strt, 2+(cs_thick-1)) # transfinite and recombination are incompatible apparently
            gmsh.model.geo.mesh.setTransfiniteCurve(l_end, 2+(cs_thick-1))

            gmsh.model.geo.synchronize()
            # the parallel sides ere divided according to the original line size and the cell size
            for j in range(len(l_ind_neg)):
                # off_neg_xy[i]
                # TODO include both sides on the computation
                disp = Point(off_neg_xy[j]).distance(Point(off_neg_xy[j+1]))
                disn = Point(off_pos_xy[j]).distance(Point(off_pos_xy[j+1]))

                cs_line = self.gdf_line.loc[i, 'cs']
                div = round(mean([disp, disn])/cs_line)
                gmsh.model.geo.mesh.setTransfiniteCurve(l_ind_neg[j], div+1)#2#div+1
                gmsh.model.geo.mesh.setTransfiniteCurve(l_ind_pos[j], div+1)#2#div+1
            #lets setup the Surface
            #set it as a transfinite surface
            gmsh.model.geo.mesh.setTransfiniteSurface(ind_s_buff, "Left", [p_ind_neg[0], p_ind_neg[-1], p_ind_pos[-1], p_ind_pos[0]])

            # to make quad grids
            #TODO write here the algorithm names
            gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 0)# 0, 2,  3
            gmsh.model.geo.mesh.setRecombine(2, ind_s_buff)
            gmsh.model.geo.synchronize()

            gmsh.model.mesh.setAlgorithm(2, ind_s_buff, 8) #frontal delunay quads

        return self.c_ind_buf

    def convert_to_points_for_size_fields(self):
        '''
        This function converts the line geometry to points for cell size fields.

        Returns
        -------
        gdf_coord : geopandas.GeoDataFrame
            Geodataframe with the point geometry.
        '''
        # TODO fix name, and include this function inside create line and poly etc
        # Implementation of conversion to points for threshold fields
        if 'cs' not in self.gdf_line.columns:
            available_cols = list(self.gdf_line.columns)
            raise ValueError(
                f"GeoDataFrame must have 'cs' column for cell sizes. Available columns: {available_cols}"
            )
        gdf2 = self.gdf_line.copy()
        # apply to have different point densities in different segments of drn
        gdf2['geometry'] = self.gdf_line.apply(lambda x: x.geometry.simplify(x.cs), axis=1)
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
            lambda x: gmsh.model.geo.addPoint(x.x, x.y, 0, x.cs), axis=1).astype(int)
        self.gf_coord = gdf_coord
        return gdf_coord
