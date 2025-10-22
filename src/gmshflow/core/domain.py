"""Domain mesh creation and management for GMSHFlow."""

import os
from typing import List, Optional, Tuple

import geopandas as gpd
import gmsh
import numpy as np
import pandas as pd
import shapely
import shapely.ops
from shapely.geometry import Point


class GmshMeshDomain:
    """Domain mesh creation and management for GMSH modeling.

    This class handles the creation and configuration of mesh domains in GMSH,
    including domain geometry preparation, field setup, and mesh export functionality.

    Args:
        name: Name identifier for the domain.
        gdf_dom: GeoDataFrame containing the domain geometry (typically polygons).
        cs_dom: Default cell size for the domain. If None, cell size should be
            specified in a 'cs' column of the geodataframe.

    Attributes:
        name: Domain name identifier.
        cs_dom: Default cell size for meshing.
        gdf_dom: Domain geometry GeoDataFrame.
        gdf_list: List of additional geometries to include in domain shape.
        shp_dom: Final processed domain shape.
        c_ind: Index of the curve loop for the outer domain.
        ind_min_field: Index of the minimum field for mesh sizing.
        loop_list_surface_dom: List of internal loops for domain surface.
        ind_s_dom: Index of the domain surface.
        ind_embed_lines: List of embedded line indices.
        ind_embed_points: List of embedded point indices.

    Example:
        >>> domain_gdf = gpd.read_file("domain.shp")
        >>> domain = GmshMeshDomain("aquifer", domain_gdf, cs_dom=100.0)
        >>> domain.prepare_mesh_domain()
        >>> domain.create_domain_loop_from_poly()
    """
    def __init__(self, name: str, gdf_dom: gpd.GeoDataFrame, cs_dom: Optional[float] = None) -> None:
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

    def add_domain_polygon_geometry(self, gdf_dom_geom: gpd.GeoDataFrame) -> None:
        """Add additional polygon geometries to extend the domain.

        This function adds extra geometries to the domain shape definition.
        The geometries will be used when preparing the mesh domain to extend
        or modify the domain boundaries.

        Args:
            gdf_dom_geom: GeoDataFrame containing polygon geometries to add
                to the domain. Must contain only Polygon geometries.

        Raises:
            TypeError: If the input is not a GeoDataFrame.
            ValueError: If geometries are not all Polygon type.

        Example:
            >>> additional_polys = gpd.read_file("extensions.shp")
            >>> domain.add_domain_polygon_geometry(additional_polys)
        """

        if not isinstance(gdf_dom_geom, gpd.GeoDataFrame):
            raise TypeError(f'Expected GeoDataFrame, got {type(gdf_dom_geom).__name__}')

        if not all(gdf_dom_geom.geom_type == 'Polygon'):
            invalid_types = gdf_dom_geom.geom_type[gdf_dom_geom.geom_type != 'Polygon'].unique()
            raise ValueError(
                f'All geometries must be Polygon type. Found invalid types: {invalid_types}'
            )

        self.gdf_list.append(gdf_dom_geom)

    def prepare_mesh_domain(self, mesh_area: int = 1,
                                  gdf_list: List[gpd.GeoDataFrame] = None,
                                  min_overlap: float = 1.0,
                                  meshing_buff_mult: float = 1.0) -> None:
        """Prepare domain geometry for meshing with optional buffering and simplification.

        This function processes the domain geometry to create an optimal shape for meshing.
        It can simplify boundaries, add buffer zones, and incorporate additional geometries.

        Args:
            mesh_area: Domain boundary type:
                - 0: Convex hull (smoothest boundary)
                - 1: Oriented envelope (rectangular boundary aligned with data)
                - 2: Bounding box (axis-aligned rectangular boundary)
            gdf_list: List of additional GeoDataFrames to incorporate into domain shape.
                If None, uses self.gdf_list.
            min_overlap: Minimum overlap fraction (0-1) between domain and additional
                geometries required for inclusion.
            meshing_buff_mult: Buffer multiplier applied to cell size for domain boundary.

        Example:
            >>> domain.prepare_mesh_domain(mesh_area=1, meshing_buff_mult=1.5)
            >>> # Creates oriented envelope with 1.5x cell size buffer
        """

        if gdf_list is None:
            gdf_list = self.gdf_list

        shp_dom = self.gdf_dom.geometry[0].simplify(self.cs_dom/2)
        if gdf_list != []:
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

    def create_domain_loop_from_poly(self) -> List[int]:
        """Create GMSH curve loop from polygon geometry of the domain.

        Converts the domain polygon boundary into GMSH points, lines, and curve loops
        that can be used for surface creation and meshing.

        Returns:
            List containing the index of the curve loop for the outer domain.

        Raises:
            AssertionError: If domain shape has not been defined via prepare_mesh_domain().

        Example:
            >>> domain.prepare_mesh_domain()
            >>> loop_indices = domain.create_domain_loop_from_poly()
            >>> print(f"Created curve loop {loop_indices[0]}")
        """

        if self.shp_dom is None:
            raise RuntimeError(
                'Domain shape is not defined. Call prepare_mesh_domain() first.'
            )
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

    def create_exponential_field(self, df_points: pd.DataFrame, fac: float = 1.1) -> Tuple[int, List[int]]:
        """Create exponential mesh size field based on point locations.

        Sets up a GMSH field that grows mesh size exponentially with distance
        from specified points. This provides smooth size transitions while
        maintaining fine resolution near important features.

        Args:
            df_points: DataFrame with columns 'cs' (cell size) and 'id_gmsh'
                (GMSH point IDs) defining control points for mesh sizing.
            fac: Growth factor for exponential field (typically 1.1-1.3).
                Higher values create more aggressive size transitions.

        Returns:
            Tuple of (minimum_field_index, list_of_exponential_field_indices).

        Raises:
            AssertionError: If fac <= 1.0 or required columns are missing.

        Example:
            >>> control_points = pd.DataFrame({
            ...     'cs': [10, 20, 50],
            ...     'id_gmsh': [1, 2, 3]
            ... })
            >>> min_field, exp_fields = domain.create_exponential_field(control_points, fac=1.2)
        """
        # Implementation of exponential field creation
        if fac <= 1.0:
            raise ValueError(
                f'Growth factor must be > 1.0 (typically 1.1-1.3). Got: {fac}'
            )

        if 'cs' not in df_points.columns:
            available_cols = list(df_points.columns)
            raise ValueError(
                f"DataFrame must have 'cs' column for cell sizes. Available columns: {available_cols}"
            )
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
        if fac <= 1.0:
            raise ValueError(
                f'Growth factor must be > 1.0 (typically 1.1-1.3). Got: {fac}'
            )

        if 'cs' not in df_points.columns:
            available_cols = list(df_points.columns)
            raise ValueError(
                f"DataFrame must have 'cs' column for cell sizes. Available columns: {available_cols}"
            )
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
        if self.ind_min_field is None:
            raise RuntimeError(
                'Minimum field is not defined. Call create_exponential_field() or '
                'create_linear_threshold_field() first.'
            )
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
            # check loop_id is not in the list
            if loop_id in self.loop_list_surface_dom:
                raise ValueError(f'Loop ID {loop_id} has already been added to the domain surface')
            self.loop_list_surface_dom.append(loop_id)
            print('Internal loop added successfully')

    def create_domain_surface(self):
        '''
        This function creates the domain surface.
        '''
        # Implementation of domain surface creation
        if self.c_ind is None:
            raise RuntimeError(
                'Domain loop is not defined. Call create_domain_loop_from_poly() first.'
            )
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
        # Validate inputs before any GMSH operations
        if not surface_ids:
            raise ValueError("No valid surfaces provided for Voronoi export")

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
