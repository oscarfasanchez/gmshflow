# GMSHFLOW is a wrapper for GMSH that allows the user to create a mesh for a
# a MODFLOW6 groundwater flow model. 
# The user provides shapefiles in geopandas format for the domain, boundary conditions, and other features
# that are used to create the mesh.

# this script its a refactored version of the original script that was created by the author
# with the name of voro_barrier_gmsh.py but using OOP principles

# Import libraries
import geopandas as gpd
import numpy as np
import pandas as pd
import os
import gmsh
geo = gmsh.model.geo
import shapely
from shapely.geometry import Point

def merge_many_multilinestring_into_one_linestring(gdf):
    # This function merges many multilinestrings into one linestring
    # check if the geodataframe has multilinestrings
    assert any(gdf.geom_type == 'MultiLineString'), 'The geodataframe must have multilinestrings'
    # explode the multilinestrings
    gdf2 = gdf.explode().reset_index()
    # group the linestrings by the index
    gdf.geometry = gdf2.groupby('index')['geometry'].apply(
        lambda x: shapely.ops.linemerge(list(x), directed=True))
    #check that there are no Multilinestrings anymore and, that the number of lines is the same as the original
    assert all(gdf.geom_type == 'LineString'), 'The geodataframe must have only linestrings, then the algorithm didnt work'

    return gdf

class GmshModel:
    def __init__(self, name):
        gmsh.initialize()
        self.name = name
        gmsh.model.add(name)
        # self.geo = gmsh.model.geo

    def finalize(self):
        gmsh.finalize()

    def synchronize(self):
        gmsh.model.geo.synchronize()

    def generate_mesh(self, dimension=2):
        gmsh.model.mesh.generate(dimension)

    def write(self, filename):
        gmsh.write(filename)

    def run_gui(self):
        gmsh.fltk.run()
# TODO check if its suitable to include the following functions in the class
#  set_field, set_mesh_size_from_geometries, add_internal_loop, create_domain_surface, export_to_voronoi
#  add_embedded_lines, add_embedded_points, create_domain_loop_from_poly, create_exponential_field, create_linear_threshold_field
#  and adding in the other objects the code to pass that dara to the gmsh model here
#  for now I just want something that works, I will refactor later
# also write get functions for the internal variables of the classes to avoid direct access to them.    

#Also need to standarize variable names and define clearly the system to store important variables.

#define the class
class GmshMeshDomain:
    #define the constructor using the domain, the domain is a 
    #geopandas dataframe that contains a column with the cell size called 'cs' 
    def __init__(self, name, gdf_dom:gpd.GeoDataFrame, cs_dom:float = None):
        self.name = name
        self.cs_dom = cs_dom
        self.gdf_list=[]# list of geopandas dataframes with the geometries to add to the domain shape
        self.gdf_dom = gdf_dom
        self.c_ind = None # index of the curve loop of the outer domain
        self.shp_dom = None # final shape of the domain
        self.ind_min_field = None # index of the minimum field
        self.loop_list_surface_dom = [] # list of internal loops to define the surface of the domain
        self.ind_s_dom = None # index of the domain surface
        


    def add_domain_polygon_geometry(self, gdf_dom_geom):
        # This function adds extra geometries to the domain and adds them to the list of gdfs.
        #  The function returns a list of geopandas dataframes with the domain geometries.
        #check that gdf_dom is a geopandas dataframe polygon of one elementand that it has a column called cs
        if not isinstance(gdf_dom_geom, gpd.GeoDataFrame):
            raise TypeError('The geometry must be a geopandas dataframe')
        if not all(gdf_dom_geom.geom_type == 'Polygon'):
            raise ValueError('All geometries in the domain must be of type Polygon')
       
        self.gdf_list.append(gdf_dom_geom)


    
    def prepare_mesh_domain(self, mesh_area=1,
                                  gdf_list=[], min_overlap=1,
                                  meshing_buff_mult=1):
               
        # This function prepares the domain for meshing by either creating a buffer around the domain, simplifying
        # the domain, or both. The buffer is also affected for additional geopandas geometries that are provided and
        # maybe extend the domain. The function returns a geopandas dataframe with the domain geometry.
        # mesh_area: 0=convex hull, 1=oriented envelope, 2=bounding box
        
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
        
        assert self.shp_dom is not None, 'The domain shape is not defined'
        # lets add the outer points of the domain, later the lines, and finally area
        shp_dom_xy = list(zip(self.shp_dom.exterior.xy[0], self.shp_dom.exterior.xy[1]))
        p_ind = []
        for xy in shp_dom_xy[:-1]:  # to avoid repeated points
            ind = geo.addPoint(xy[0], xy[1], 0, self.cs_dom)
            p_ind.append(ind)
        # p_ind = [x for x in range(1, i+1)]
        l_ind = []
        for i in range(len(p_ind)):
            ind = geo.addLine(p_ind[i-1], p_ind[i])
            l_ind.append(ind)
        
        c_ind = []
        ind = geo.addCurveLoop(l_ind, 1)
        c_ind.append(ind)
        self.c_ind = c_ind
        return c_ind

    def create_exponential_field(self, df_points, fac=1.1):
        # Implementation of exponential field creation
        assert fac > 1, 'fac must be higher than 1, usually between 1.1 and 1.3'
        assert 'cs' in df_points.columns, 'df_points must have cell size column(cs)'
        gmsh.model.geo.synchronize()
        ind_thres_field = []
        #define a distance and exponential field for each cell size 
        for cs in df_points.cs.unique():
            df_points_cs = df_points.loc[df_points.cs==cs]
            ind_dist_field = gmsh.model.mesh.field.add("Distance")  # TODO check to get var
            gmsh.model.mesh.field.setNumbers(
                ind_dist_field, "PointsList", df_points_cs['id_gmsh'].to_list())  # [:1]
            ind_thres_field.append(gmsh.model.mesh.field.add("MathEval"))  # TODO check to get var
            d = ind_dist_field
            #define the function
            gmsh.model.mesh.field.setString(
                ind_thres_field[-1],
                "F", f'{cs}*{fac}^(Log(1+F{d}*2*Log({fac})/{cs})/Log({fac}))')
        # get the minimum cell size field from the intersection of all the fields
        ind_min_field = gmsh.model.mesh.field.add("Min")
        # set it a s gmsh field
        gmsh.model.mesh.field.setNumbers(ind_min_field, "FieldsList", ind_thres_field)
        self.ind_min_field = ind_min_field
        return ind_min_field, ind_thres_field

    def create_linear_threshold_field(self, df_points, fac=1.2):
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
        # Implementation of field setting
        assert self.ind_min_field is not None, 'The minimum field is not defined'
        gmsh.model.mesh.field.setAsBackgroundMesh(self.ind_min_field)
        gmsh.model.geo.synchronize()

    def set_mesh_size_from_geometries(self, use_boundaries=False, use_points=False):
        #by default in gmsh the mesh size is set from the geometry and linearly interpolated
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

    def add_internal_loops(self, loop_ids:list):
        # Implementation of internal loop addition
        for loop_id in loop_ids:
            #check loop_id is not in the list
            assert loop_id not in self.loop_list_surface_dom, 'The loop was already added'
            self.loop_list_surface_dom.append(loop_id)
            print('Internal loop added successfully')
  

    def create_domain_surface(self):
        # Implementation of domain surface creation
        assert self.c_ind is not None, 'The domain loop is not defined'
        if self.loop_list_surface_dom is None:
            print('Warning: no internal loops defined')
        gmsh.model.geo.synchronize()
        self.ind_s_dom = geo.addPlaneSurface(self.c_ind+self.loop_list_surface_dom)
        gmsh.model.geo.synchronize()
        return self.ind_s_dom

    def export_to_voronoi(self, ws, name:str, surface_ids:list, export_modflow=False, int_org_dom=True):
                        
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
        #organize the coordinates in an array more suitable 
        triang_node_coords = nodeCoords.reshape(len(nodeCoords) // 3, 3)
        triang_node_coords = [[x[0], x[1]] for x in triang_node_coords]
        #delete duplicates caused by including the boundary of all surfaces
        triang_node_coords = np.unique(triang_node_coords, axis=0)
        #create the geopandas dataframe with the coordinates to use the geopandas voronoi function
        triang_points = [Point(i) for i in triang_node_coords]
        trian_p_gpd = gpd.GeoDataFrame(geometry=triang_points, crs=self.gdf_dom.crs)
        voro_gpd = trian_p_gpd.voronoi_polygons(tolerance=1e-3)
        #clip the voronoi polygons to keep the extent reasonable, and not distorting the cell centers too much in the boundary
        voro_gpd = voro_gpd.clip(self.shp_dom.buffer(self.cs_dom / 3, join_style='mitre'))
        if int_org_dom:
            voro_gpd = voro_gpd.loc[voro_gpd.intersects(self.gdf_dom.geometry[0])]
        #export the voronoi polygons to a shapefile
        voro_gpd.to_file(os.path.join(ws, f'{name}.shp'))# TODO maybe using by default the mesh name
        # trian_p_gpd.to_file(os.path.join(ws, f'{name}_points.shp'))
        if export_modflow:
            import flopy
            verts, iverts = flopy.utils.cvfdutil.shapefile_to_cvfd(
                os.path.join(ws, f'{name}.shp'), verbose=True, duplicate_decimals=3, skip_hanging_node_check=True)
            gridprops_dic = flopy.utils.cvfdutil.get_disv_gridprops(verts, iverts, xcyc=None)
            return voro_gpd, gridprops_dic
        else:
            return voro_gpd
        
    

    def add_embedded_lines(self, id_line_list:list, surface_id:int==None):
        line_dim=1
        surface_dim=2
        if surface_id is None:
            surface_id = self.ind_s_dom
            print('Warning: the domain surface was used as the target surface')
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(line_dim, id_line_list, surface_dim, surface_id)
        pass

    def add_embedded_points(self, id_point_list:list, surface_id:int==None):
        point_dim=0
        surface_dim=2
        if surface_id is None:
            surface_id = self.ind_s_dom
            print('Warning: the domain surface was used as the target surface')
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(point_dim, id_point_list, surface_dim, surface_id)




class PolyGeometryHandler:
    def __init__(self, cs_poly=None):
        self.gdf_poly = None
        self.cs_poly = cs_poly
        self.gdf_coord = None
        #TODO include the following attributes inside self.gdf_poly
        self.s_ind = []

    def set_gpd_poly(self, gdf_poly: gpd.GeoDataFrame):
        # check it is a polygon dataframe
        assert all(gdf_poly.geom_type == 'Polygon'), 'All geometries must be of type Polygon'
        #check that it has a column called cs or cs_poly is not None
        assert 'cs' in gdf_poly.columns or self.cs_poly is not None, 'The geodataframe must have a cell size column or cs_poly must be defined'
        self.gdf_poly = gdf_poly
        #simplify the geometries
        if self.cs_poly is not None:
            self.gdf_poly.geometry = self.gdf_poly.geometry.simplify(self.cs_poly/2)
            self.gdf_poly['cs'] = self.cs_poly
        else:
            self.gdf_poly.geometry = self.gdf_poly.apply(lambda x: x.geometry.simplify(x.cs/2), axis=1)

        
        

    
    def create_loop_from_poly(self, def_surf=False):
        # This function creates a loop from a polygon geometry
        assert self.gdf_poly is not None, 'The polygon geometry is not defined'
        # lets get the points of the polygons, then the lines, and finally the area
        for i in self.gdf_poly.index:
            poly = self.gdf_poly.loc[i, 'geometry']
            poly_xy = list(zip(poly.exterior.xy[0], poly.exterior.xy[1]))
            p_ind = []
            for xy in poly_xy[:-1]:  # to avoid repeated points
                ind = geo.addPoint(xy[0], xy[1], 0, self.gdf_poly.loc[i, 'cs'])# TODO check
                p_ind.append(ind)
            
            l_ind = []
            for j in range(len(p_ind)):
                ind = geo.addLine(p_ind[j-1], p_ind[j])
                l_ind.append(ind)

            c_ind = []
            ind = geo.addCurveLoop(l_ind)
            c_ind.append(ind)
               

            if def_surf:
                ind_s = geo.addPlaneSurface(c_ind)
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
                ind = geo.addPoint(xy[0], xy[1], 0, self.cs_poly)
                p_ind_pos.append(ind)

            p_ind_neg = []
            for xy in off_neg_xy[:-1]:  # to avoid repeated points
                ind = geo.addPoint(xy[0], xy[1], 0, self.cs_poly)
                p_ind_neg.append(ind)

            # define lines
            l_ind_pos = []
            for j in range(len(p_ind_pos)):
                ind = geo.addLine(p_ind_pos[j-1], p_ind_pos[j])
                l_ind_pos.append(ind)

            l_ind_neg = []
            for i in range(len(p_ind_neg)):
                ind = geo.addLine(p_ind_neg[j-1], p_ind_neg[j])
                l_ind_neg.append(ind)

            # define loops
            loop_list_neg = l_ind_neg
            loop_list_pos = l_ind_pos

            c_ind_buf_neg.append(geo.addCurveLoop(loop_list_neg))
            c_ind_buf_pos.append(geo.addCurveLoop(loop_list_pos))
            self.c_ind.append([c_ind_buf_neg, c_ind_buf_pos])  
            print('Loop list created')
            if def_surf:
                # define surfaces
                ind_s_buff.append(geo.addPlaneSurface([c_ind_buf_pos[-1], c_ind_buf_neg[-1]]))
                ind_s_mid.append(geo.addPlaneSurface([c_ind_buf_neg[-1]]))
                self.s_ind.append([ind_s_mid, ind_s_buff])

                geo.synchronize()
                # set quads for buffer area
                geo.mesh.setRecombine(2, ind_s_buff[-1])
                geo.mesh.setAlgorithm(2, ind_s_buff[-1], 8)
                geo.mesh.setAlgorithm(2, ind_s_mid[-1], 6)
                self.gdf_poly.loc[self.gdf_poly.shape[0], ['layer', 'geometry']] = [f'off_pos{i}', off_pos]
                self.gdf_poly.loc[self.gdf_poly.shape[0], ['layer', 'geometry']] = [f'off_neg{i}', off_neg]
        print('finished: create_surfacegrid_from_buffer_poly')
        return self.gdf_poly, c_ind_buf_pos, c_ind_buf_neg, ind_s_buff, ind_s_mid

    def convert_to_points_for_threshold_fields(self) :
        # Implementation of conversion to points for threshold fields
        assert 'cs' in self.gdf_poly.columns, 'gdf must have cell size column(cs)'
        gdf2 = self.gdf_poly.copy()
        # apply to have different point densities in different segments of drn
        gdf2['geometry'] = self.gdf_poly.apply(lambda x: x.geometry.simplify(x.cs), axis=1) #  TODO cs_bar/2??
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
            lambda x: geo.addPoint(x.x, x.y, 0, x.cs), axis=1)
        self.gf_coord = gdf_coord
        return gdf_coord



    #TODO implement the following functions: delete...curves, delete...surfaces

class LineGeometryHandler:
    def __init__(self, cs_line=None):
        self.gdf_line = None
        self.cs_line = cs_line
        self.gdf_coord = None
        self.ind_s_buff = []
        self.c_ind_buf = []
        
        self.l_ind_list = []


    def set_gpd_line(self, gdf_line: gpd.GeoDataFrame):
        # check it is a line dataframe
        assert all((gdf_line.geom_type == 'LineString')|(gdf_line.geom_type == 'MultiLineString')), 'All geometries must be of type LineString'
        #check that it has a column called cs or cs_line is not None
        
            
        assert 'cs' in gdf_line.columns or self.cs_line is not None, 'The geodataframe must have a cell size column or cs_line must be defined'
        if any(gdf_line.geom_type == 'MultiLineString'):
            self.gdf_line = gdf_line.explode().reset_index()
        else:  
            self.gdf_line = gdf_line
        #simplify the geometries
        if self.cs_line is not None:
            self.gdf_line.geometry = self.gdf_line.geometry.simplify(self.cs_line/2)
            self.gdf_line['cs'] = self.cs_line
        else:
            self.gdf_line.geometry = self.gdf_line.apply(lambda x: x.geometry.simplify(x.cs/2), axis=1)

    def create_line_from_line(self):
        # This function creates a gmsh point and line from a line geometry
        assert self.gdf_line is not None, 'The line geometry is not defined'
        # lets get the points of the polygons, then the lines, and finally the area
        for i in self.gdf_line.index:
            line = self.gdf_line.loc[i, 'geometry']
            line_xy = list(zip(line.xy[0], line.xy[1]))
            p_ind = []
            for xy in line_xy[:]:  # to avoid repeated points
                ind = geo.addPoint(xy[0], xy[1], 0, self.gdf_line.loc[i, 'cs'])# TODO check
                p_ind.append(ind)
            
            l_ind = []
            for j in range(len(p_ind)-1):
                ind = geo.addLine(p_ind[j], p_ind[j+1])
                l_ind.append(ind)
            self.l_ind_list += l_ind
        # self.gdf _line['p_ind'] = p_ind

        return self.l_ind_list


    def create_surfacegrid_from_buffer_line(self, cs_thick:int = 1, simpl_fac=1.5, def_surf=True):
        # lets simplify the barrier to ease the quad meshing
        # Failed attempt to smooth angles, better just enforce 1 or 2 to avoid 
        # misuse of it with poor meshes
        print('starting: create_surfacegrid_from_buffer_line')
        print('Warning, this function expects that different feature lines are not intersecting')
        #to get decent quality mesh, the offset should be at least 1 cell and considerable simplification
        assert cs_thick >= 1 and cs_thick <= 2, 'this function was designed only for offset of 1 or 2 cells, if used for more it results in bad quality meshes'
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
                ind = geo.addPoint(xy[0], xy[1], 0, cs_line)
                p_ind_pos.append(ind)
            
            p_ind_neg = []
            for xy in off_neg_xy:
                ind = geo.addPoint(xy[0], xy[1], 0, cs_line)
                p_ind_neg.append(ind)
            # define lines
            l_ind_pos = []
            for j in range(len(p_ind_pos)-1):
                ind = geo.addLine(p_ind_pos[j], p_ind_pos[j+1])
                l_ind_pos.append(ind)

            l_ind_neg = []
            for j in range(len(p_ind_neg)-1):
                ind = geo.addLine(p_ind_neg[j], p_ind_neg[j+1])
                l_ind_neg.append(ind)
            
            l_strt = geo.addLine(p_ind_pos[0], p_ind_neg[0])
            l_end = geo.addLine(p_ind_neg[-1], p_ind_pos[-1])
            # define loops, carefully, in the right order
            loop_list = l_ind_neg+ [l_end]+ list(np.array(l_ind_pos[::-1])*-1)+ [l_strt]
            c_ind_buf = geo.addCurveLoop(loop_list)
            # define surfaces
            ind_s_buff = geo.addPlaneSurface([c_ind_buf])
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
                from statistics import mean
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
            #save the offset geometry on the original line gdf
            # self.gdf_line.loc[self.gdf_line.shape[0], ['layer', 'geometry', 'cs']] = [f'off_pos_{i}', off_pos, self.gdf_line.loc[i, 'cs']]
            # self.gdf_line.loc[self.gdf_line.shape[0], ['layer', 'geometry', 'cs']] = [f'off_neg_{i}', off_neg, self.gdf_line.loc[i, 'cs']]
            # self.ind_s_buff = ind_s_buff
            # self.c_ind_buf = c_ind_buf
            return self.c_ind_buf#ind_s_buff, c_ind_buf



    def convert_to_points_for_threshold_fields(self):
        # Implementation of conversion to points for threshold fields
        assert 'cs' in self.gdf_line.columns, 'gdf must have cell size column(cs)'
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
            lambda x: geo.addPoint(x.x, x.y, 0, x.cs), axis=1)
        self.gf_coord = gdf_coord
        return gdf_coord

class PointGeometryHandler:
    def __init__(self, cs_point=None):
        self.gdf_point = None
        self.cs_point = cs_point
        self.gdf_coord = None

    def set_gdf_point(self, gdf_point: gpd.GeoDataFrame):
        # check it is a point dataframe
        assert all(gdf_point.geom_type == 'Point'), 'All geometries must be of type Point'
        #check that it has a column called cs or cs_point is not None
        assert 'cs' in gdf_point.columns or self.cs_point is not None, 'The geodataframe must have a cell size column or cs_point must be defined'
        self.gdf_point = gdf_point
        #simplify the geometries
        if self.cs_point is not None:
            self.gdf_point['cs'] = self.cs_point
        
    def create_point_from_point(self, df_coord=False):
        # This function creates a gmsh point from a point geometry
        assert self.gdf_point is not None, 'The point geometry is not defined'
        # lets get the points of the polygons, then the lines, and finally the area
        p_ind = []
        for i in self.gdf_point.index:
            point = self.gdf_point.loc[i]
            ind = geo.addPoint(point.geometry.x, point.geometry.y, 0, self.gdf_point.loc[i, 'cs'])# TODO check
            self.gdf_point.loc[i, 'p_ind'] = ind
            p_ind.append(ind)
        if df_coord:
            self.gdf_coord = self.gdf_point.get_coordinates().reset_index(names='ind')
        return p_ind



if __name__ == "__main__":
    

    # Import all shapefiles
    wd_shp = os.path.join('..', '02_shp')
    wd_exe = os.path.join('..', '03_exe')
    wd_model = os.path.join('..', '01_model')

    # Domain
    gdf_dom = gpd.read_file(os.path.join(wd_shp, 'ad0.shp'))
    gdf_dom_0 = gdf_dom

    # River
    gdf_drn = gpd.read_file(os.path.join(wd_shp, 'drainage_network_custom.shp'))
    gdf_drn_0 = gdf_drn

    # Barrier
    gdf_bar = gpd.read_file(os.path.join(wd_shp, 'updated_split_barrier_obs.shp'))
    gdf_bar_merge = gpd.read_file(os.path.join(wd_shp, 'updated_merge_barrier.shp'))

    # Faults
    gdf_fault = gpd.read_file(os.path.join(wd_shp, 'faults.shp'))
    gdf_fault = gdf_fault.clip(gdf_dom.geometry[0])

    # Pits
    gdf_pit = gpd.read_file(os.path.join(wd_shp, 'pits_no_ovlp.shp'))  # 'pits.shp'
    

    # Define cell_size in each area
    cs_dom = 2000.0
    cs_riv = 125.0
    cs_bar = 31.25
    cs_pit = 125.0

    # Add cell size column
    gdf_drn['cs'] = cs_riv
    gdf_drn.loc[4, 'cs'] = cs_riv / 4
    gdf_pit['cs'] = cs_pit
    gdf_fault['cs'] = cs_pit

    # Initialize GmshModel and GmshMeshDomain
    gmsh_model = GmshModel("yandi_gmsh")
    mesh_domain = GmshMeshDomain("domain", gdf_dom, cs_dom)

    # Prepare mesh domain
    mesh_domain.prepare_mesh_domain(mesh_area=1, gdf_list=[gdf_drn, gdf_pit, gdf_fault], min_overlap=0.5)

    # Create domain loop
    c_ind = mesh_domain.create_domain_loop_from_poly()

    # Add domain polygon geometry

    pit = PolyGeometryHandler()
    pit.set_gpd_poly(gdf_pit)
    # c_ind_pit = pit.create_loop_from_poly(def_surf=True)

    # Add river line geometry
    riv = LineGeometryHandler()
    riv.set_gpd_line(gdf_drn)
    riv.create_line_from_line()

    # Add barrier line geometry
    bar = LineGeometryHandler(cs_bar)
    bar.set_gpd_line(gdf_bar_merge)
    # Create the line and surface for the buffer around the barrier(to get a cute barrier)
    c_ind_bar = bar.create_surfacegrid_from_buffer_line(cs_thick=1)

    # Add fault line geometry
    flt = LineGeometryHandler()
    flt.set_gpd_line(gdf_fault)
    flt.create_line_from_line()

    
    # Convert to points for threshold fields
    gdf_drn_coord = riv.convert_to_points_for_threshold_fields()
    gdf_pit_coord = pit.convert_to_points_for_threshold_fields()
    gdf_fault_coord = flt.convert_to_points_for_threshold_fields()
    gdf_bar_coord = bar.convert_to_points_for_threshold_fields()

    
    # Create exponential field
    gdf_thr_coord = pd.concat([gdf_pit_coord, gdf_drn_coord, gdf_fault_coord, gdf_bar_coord])#, 
    mesh_domain.create_exponential_field(gdf_thr_coord, fac=1.1)
    
    
    # set exponential field as background mesh
    mesh_domain.set_field()
   
    #deactivate the mesh size from the geometry points
    mesh_domain.set_mesh_size_from_geometries(use_boundaries=False, use_points=False)

    # Create domain surface
    #add the loops to the domain surface
    mesh_domain.add_internal_loops(c_ind_bar)
    # mesh_domain.add_internal_loops(c_ind_pit)
    #once added the internal loops, we can create the 
    ind_s_dom = mesh_domain.create_domain_surface()

    # Generate mesh
    print('Generating the mesh...')
    gmsh_model.generate_mesh()
    print('Mesh ready')

    # Run GUI 
    gmsh_model.run_gui()
    
    # Export to Voronoi
    surf_tags = [ind_s_dom] + bar.ind_s_buff#pit.s_ind#pit.s_ind bar.ind_s_buff
    voro = mesh_domain.export_to_voronoi(wd_shp, "gdf_voro", surf_tags, export_modflow=True)

    # Write mesh to file
    gmsh_model.write("yandi_gmsh.msh")

    

    # Finalize Gmsh
    gmsh_model.finalize()

    #TODO make a check for surface geometries intersecting, in that case, ask the user
    # to fix it
    #TODO get mesh quality metrics


