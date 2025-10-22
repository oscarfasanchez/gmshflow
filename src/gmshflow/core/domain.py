"""Domain mesh creation and management for GMSHFlow."""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import gmsh
import shapely
import shapely.ops
from shapely.geometry import Point


class GmshMeshDomain:
    ''' 
    Class to create a domain for meshing the domain in GMSH.
    
    Parameters
    ----------
    name : str
        Name of the domain.
    gdf_dom : geopandas.GeoDataFrame
        Geodataframe with the domain geometry.
    cs_dom : float, optional
        Cell size of the domain. The default is None, because in some case the cell size is defined
        in a field called 'cs'.
    
    '''
    #define the constructor using the domain, the domain is a 
    #geopandas dataframe that contains a column with the cell size called 'cs' 
    def __init__(self, name, gdf_dom: gpd.GeoDataFrame, cs_dom: float = None):
        self.name = name
        self.cs_dom = cs_dom
        self.gdf_list = []  # list of geopandas dataframes with the geometries to add to the domain shape
        self.gdf_dom = gdf_dom
        self.c_ind = None  # index of the curve loop of the outer domain
        self.shp_dom = None  # final shape of the domain
        self.ind_min_field = None  # index of the minimum field
        self.loop_list_surface_dom = []  # list of internal loops to define the surface of the domain
        self.ind_s_dom = None  # index of the domain surface
        self.ind_embed_lines = []
        self.ind_embed_points = []

    def add_domain_polygon_geometry(self, gdf_dom_geom):
        '''
        This function adds extra geometries to the domain and adds them to the list of gdfs.
        The function returns a list of geopandas dataframes with the domain geometries.
        check that gdf_dom is a geopandas dataframe polygon of one elementand that it has a column called cs

        Parameters
        ----------
        gdf_dom_geom : Geodataframe
            Geometry of the domain, only one feature

        Raises
        ------
        TypeError
            in case a geodataframe is not provided.
        ValueError
            in case a polygon is not provided.

        Returns
        -------
        None.

        '''
        
        if not isinstance(gdf_dom_geom, gpd.GeoDataFrame):
            raise TypeError('The geometry must be a geopandas dataframe')
        if not all(gdf_dom_geom.geom_type == 'Polygon'):
            raise ValueError('All geometries in the domain must be of type Polygon')
       
        self.gdf_list.append(gdf_dom_geom)

    def prepare_mesh_domain(self, mesh_area=1,
                                  gdf_list=[], min_overlap=1,
                                  meshing_buff_mult=1):
        '''
        This function prepares the domain for meshing by either creating a buffer around the domain, simplifying
        the domain, or both. The buffer is also affected for additional geopandas geometries that are provided to
        extend the domain. The function returns a geopandas dataframe with the domain geometry.
        
        Parameters
        ----------
        mesh_area : int, optional
            0=convex hull, 1=oriented envelope, 2=bounding box. The default is 1.
        gdf_list : list, optional
            List of geopandas dataframes with the geometries to add to the domain shape. The default is [].
        min_overlap : float, optional
            Minimum overlap between the domain and the additional geometries. The default is 1.
        meshing_buff_mult : float, optional
            Multiplier for the meshing buffer. The default is 1.
        '''
               
        shp_dom = self.gdf_dom.geometry[0].simplify(self.cs_dom/2)
        if self.gdf_list != []:
            for gdf in gdf_list:
                # lets simplify the domain and delete features not intersecting the domain
                gdf = gdf.loc[gdf.intersects(self.gdf_dom.geometry[0])]

                gdf['geometry'] = gdf.geometry.simplify(self.cs_dom/2)
                #check the amount of overlap
                gdf['cgeometry'] = gdf.geometry.clip(shp_dom)
                gdf['per_area'] = gdf['cgeometry'].area/gdf['geometry'].area

                #filter the gdf according to the minimum overlap
                gdf = gdf.loc[gdf['per_area'] >= min_overlap]
                one_multipol = shapely.ops.unary_union(gdf.geometry)
                shp_dom = shapely.ops.unary_union([one_multipol,
                                                    shp_dom])
                
        # to smooth triangulation, simpler boundaries makes better mesh
        if mesh_area == 0:
            shp_dom = shapely.convex_hull(shp_dom).buffer(
                self.cs_dom*meshing_buff_mult, join_style='mitre')
        elif mesh_area == 1:
            shp_dom = shapely.oriented_envelope(shp_dom).buffer(
                self.cs_dom*meshing_buff_mult, join_style='mitre')  # could be avoided
        else:
            shp_dom = shapely.envelope(shp_dom).buffer(
                self.cs_dom*meshing_buff_mult, join_style='mitre')
        #TODO better to keep it in the gdf variable?
        self.shp_dom = shp_dom
        
    def create_domain_loop_from_poly(self):
        '''
        This function creates a loop from the polygon geometry of the domain.
        The function returns the index of the curve loop of the outer domain.
        
        Returns
        -------
        c_ind : int
            Index of the curve loop of the outer domain.
        '''
        
        assert self.shp_dom is not None, 'The domain shape is not defined'
        # lets add the outer points of the domain, later the lines, and finally area
        shp_dom_xy = list(zip(self.shp_dom.exterior.xy[0], self.shp_dom.exterior.xy[1]))
        p_ind = []
        for xy in shp_dom_xy[:-1]:  # to avoid repeated points
            ind = gmsh.model.geo.addPoint(xy[0], xy[1], 0, self.cs_dom)
            p_ind.append(ind)
        # p_ind = [x for x in range(1, i+1)]
        l_ind = []
        for i in range(len(p_ind)):
            ind = gmsh.model.geo.addLine(p_ind[i-1], p_ind[i])
            l_ind.append(ind)
        
        c_ind = []
        ind = gmsh.model.geo.addCurveLoop(l_ind, 1)
        c_ind.append(ind)
        self.c_ind = c_ind
        return c_ind

    def create_exponential_field(self, df_points, fac=1.1):
        '''
        This function creates an exponential field for the domain.
        
        Parameters
        ----------
        df_points : DataFrame
            Dataframe with the points to be added to the minimum cell size general field.
        fac : float, optional
            Factor that define rate of cell size growing for the exponential field.
            The default is 1.1.
            
        Returns
        -------
        ind_min_field : int
            Index of the minimum field.
        ind_exp_field : list
            List of the index of the exponential fields for each cell size.'''
        # Implementation of exponential field creation
        assert fac > 1, 'fac must be higher than 1, usually between 1.1 and 1.3'
        assert 'cs' in df_points.columns, 'df_points must have cell size column(cs)'
        gmsh.model.geo.synchronize()
        ind_exp_field = []
        #define a distance and exponential field for each cell size 
        for cs in df_points.cs.unique():
            df_points_cs = df_points.loc[df_points.cs==cs]
            ind_dist_field = gmsh.model.mesh.field.add("Distance")  # TODO check to get var
            gmsh.model.mesh.field.setNumbers(
                ind_dist_field, "PointsList", df_points_cs['id_gmsh'].to_list())  # [:1]
            ind_exp_field.append(gmsh.model.mesh.field.add("MathEval"))  # TODO check to get var
            d = ind_dist_field
            #define the function
            gmsh.model.mesh.field.setString(
                ind_exp_field[-1],
                "F", f'{cs}*{fac}^(Log(1+F{d}*2*Log({fac})/{cs})/Log({fac}))')
        # get the minimum cell size field from the intersection of all the fields
        ind_min_field = gmsh.model.mesh.field.add("Min")
        # set it a s gmsh field
        gmsh.model.mesh.field.setNumbers(ind_min_field, "FieldsList", ind_exp_field)
        self.ind_min_field = ind_min_field
        return ind_min_field, ind_exp_field

    def create_linear_threshold_field(self, df_points, fac=1.2):
        '''
        This function creates a linear field for the domain.
        
        Parameters
        ----------
        df_points : DataFrame
            Dataframe with the points to be added to the minimum cell size general field.
        fac : float, optional
            Factor that define rate of cell size growing for the linear field.
            The default is 1.2.
        
        Returns
        -------
        ind_min_field : int
            Index of the minimum field.
        ind_thres_field : list
            List of the index of the threshold fields for each cell size.'''

        # Implementation of linear field creation
        assert fac > 1, 'fac must be higher than 1, usually between 1.1 and 1.3'
        assert 'cs' in df_points.columns, 'df_points must have cell size column(cs)'
        gmsh.model.geo.synchronize()
        ind_thres_field = []
        #define a distance and exponential field for each cell size 
        for cs in df_points.cs.unique():
            df_points_cs = df_points.loc[df_points.cs==cs]
            ind_dist_field = gmsh.model.mesh.field.add("Distance")  # TODO check to get var
            gmsh.model.mesh.field.setNumbers(
                ind_dist_field, "PointsList", df_points_cs['id_gmsh'].to_list())
            ind_thres_field.append(gmsh.model.mesh.field.add("Threshold"))
            gmsh.model.mesh.field.setNumber(ind_thres_field[-1], "InField", ind_dist_field)
            gmsh.model.mesh.field.setNumber(ind_thres_field[-1], "SizeMin", cs)
            gmsh.model.mesh.field.setNumber(ind_thres_field[-1], "SizeMax", self.cs_dom)
            min_trans_cells = np.log(self.cs_dom/cs)/np.log(fac)
            min_trans_dist = cs*((fac**min_trans_cells)-1)/np.log(fac)
            gmsh.model.mesh.field.setNumber(ind_thres_field[-1], "DistMin", cs/2)
            gmsh.model.mesh.field.setNumber(ind_thres_field[-1], "DistMax", min_trans_dist/2)
        ind_min_field = gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.setNumbers(ind_min_field, "FieldsList", ind_thres_field)
        self.ind_min_field = ind_min_field
        # return ind_min_field, ind_thres_field
    
    def set_field(self):
        '''
        This function sets the minimum field as the background mesh field.
        '''
        # Implementation of field setting
        assert self.ind_min_field is not None, 'The minimum field is not defined'
        gmsh.model.mesh.field.setAsBackgroundMesh(self.ind_min_field)
        gmsh.model.geo.synchronize()

    def set_mesh_size_from_geometries(self, use_boundaries=False, use_points=False):
        '''
        This function defines if the mesh size from the geometries, or only field sizes
        because by default in gmsh the mesh size is set from the geometry and linearly
        interpolated.
        '''
        if use_boundaries:
            # Extend computation of mesh element sizes from the boundaries into the interior
            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
        else:
            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        if use_points:
            # Compute mesh element sizes from values given at geometry points
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
        else:
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        
        # gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    def add_internal_loops(self, loop_ids: list):
        '''
        This function adds internal loops to the domain surface.

        Parameters
        ----------
        loop_ids : list
            List of the internal loop ids to be added to the domain surface
        '''
        # Implementation of internal loop addition
        for loop_id in loop_ids:
            #check loop_id is not in the list
            assert loop_id not in self.loop_list_surface_dom, 'The loop was already added'
            self.loop_list_surface_dom.append(loop_id)
            print('Internal loop added successfully')

    def create_domain_surface(self):
        '''
        This function creates the domain surface.
        '''
        # Implementation of domain surface creation
        assert self.c_ind is not None, 'The domain loop is not defined'
        if self.loop_list_surface_dom is None:
            print('Warning: no internal loops defined')
        gmsh.model.geo.synchronize()
        self.ind_s_dom = gmsh.model.geo.addPlaneSurface(self.c_ind+self.loop_list_surface_dom)
        gmsh.model.geo.synchronize()
        return self.ind_s_dom

    def export_to_voronoi(self, ws, name: str, surface_ids: list, int_org_dom=True,
                          min_cell_overlap=0.5, triangle_vert_export=False):
        '''
        This function exports the domain mesh to a voronoi shapefile.
        
        Parameters
        ----------
        ws : str
            Path to the working directory.
        name : str
            Name of the voronoi shapefile.
        surface_ids : list
            List of the surface ids that you want to include in the final mesh exported.
        int_org_dom : bool, optional
            If True, the function filters the voronoi polygons that intersect the domain.
            The default is True.
        min_cell_overlap : float, optional
            Minimum overlap between the voronoi cells and the domain to be
            included in the final export when int_org_dom is True . The default is 0.5.
        triangle_vert_export : bool, optional
            If True, the function exports the triangle vertices to a shapefile. The default is False.

        '''                
        surf_tags = surface_ids
        # check if the domain surface id is included in the list of surfaces
        if self.ind_s_dom not in surf_tags:
            print('Warning: the domain surface is not included in the list of surfaces, it was added automatically')
            surf_tags.append(self.ind_s_dom)
        # if the domain surface has internal loops, warn that the loops surfaces must be included unless a hole is desired
        if self.loop_list_surface_dom != []:
            print('Warning: the domain surface has internal loops, they must be included in the list of surfaces unless a hole is desired')

        #get the coordinate of the triangle nodes, including the boundary of each surface
        nodeCoords = np.array([])
        for i in surf_tags:
            _, nodeCoords1, _ = gmsh.model.mesh.getNodes(2, i, includeBoundary=True)
            nodeCoords = np.concatenate((nodeCoords, nodeCoords1))
        if self.ind_embed_points != []:
            for i in self.ind_embed_points:
                _, nodeCoords1, _ = gmsh.model.mesh.getNodes(0, i, includeBoundary=True)
                nodeCoords = np.concatenate((nodeCoords, nodeCoords1))
        if self.ind_embed_lines != []:
            for i in self.ind_embed_lines:
                _, nodeCoords1, _ = gmsh.model.mesh.getNodes(1, i, includeBoundary=True)
                nodeCoords = np.concatenate((nodeCoords, nodeCoords1))
        #organize the coordinates in an array more suitable 
        triang_node_coords = nodeCoords.reshape(len(nodeCoords) // 3, 3)
        triang_node_coords = [[x[0], x[1]] for x in triang_node_coords]
        #delete duplicates caused by including the boundary of all surfaces
        triang_node_coords = np.unique(triang_node_coords, axis=0)
        #create the geopandas dataframe with the coordinates to use the geopandas voronoi function
        triang_points = [Point(i) for i in triang_node_coords]
        trian_p_gpd = gpd.GeoDataFrame(geometry=triang_points, crs=self.gdf_dom.crs)
        if triangle_vert_export:
            trian_p_gpd.to_file(os.path.join(ws, f'{name}_points.shp'))
        voro_gpd = trian_p_gpd.voronoi_polygons(tolerance=1e-3)
        #clip the voronoi polygons to keep the extent reasonable, and not distorting the cell centers too much in the boundary
        voro_gpd = voro_gpd.clip(self.shp_dom.buffer(self.cs_dom / 3, join_style='mitre'))
        if int_org_dom:
            voro_gpd = gpd.GeoDataFrame(geometry=voro_gpd.loc[voro_gpd.intersects(self.gdf_dom.geometry[0])], crs=self.gdf_dom.crs)
            voro_gpd_clip = voro_gpd.clip(self.gdf_dom.geometry[0]).geometry
            voro_gpd['per_area'] = voro_gpd_clip.area/voro_gpd.area
            voro_gpd = voro_gpd.loc[voro_gpd['per_area'] >= min_cell_overlap]
            #filter below threshol
        #export the voronoi polygons to a shapefile
        voro_gpd.to_file(os.path.join(ws, f'{name}.shp'))# TODO maybe using by default the mesh name
        # trian_p_gpd.to_file(os.path.join(ws, f'{name}_points.shp'))

        return voro_gpd

    def add_embedded_lines(self, id_line_list: list, surface_id: int = None):
        '''
        This function embeds lines in the domain surface.
        this is used when the cell centers are required to be aligned with the lines
        Parameters
        ----------
        id_line_list : list
            List of the line ids to be embedded in the domain surface.
        surface_id : int, optional
            Id of the surface where the lines will be embedded. The default is None.
        '''
        line_dim = 1
        surface_dim = 2
        if surface_id is None:
            surface_id = self.ind_s_dom
            print('Warning: the domain surface was used as the target surface')
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(line_dim, id_line_list, surface_dim, surface_id)
        self.ind_embed_lines += id_line_list

    def add_embedded_points(self, id_point_list: list, surface_id: int = None):
        '''
        This function embeds points in the domain surface.

        Parameters
        ----------
        id_point_list : list
            List of the point ids to be embedded in the domain surface.
        surface_id : int, optional
            Id of the surface where the points will be embedded. The default is None.

        '''
        point_dim = 0
        surface_dim = 2
        if surface_id is None:
            surface_id = self.ind_s_dom
            print('Warning: the domain surface was used as the target surface')
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(point_dim, id_point_list, surface_dim, surface_id)
        self.ind_embed_points += id_point_list