"""Core GMSH model wrapper for GMSHFlow."""

import gmsh
import pandas as pd


class GmshModel:
    '''
    Class to create a GMSH instance for meshing.

    Parameters
    ----------
    name : str
        Name of the GMSH model. 

    '''
    def __init__(self, name):
        '''
        Initializes the GMSH model.
        
        Parameters 
        ----------
        name : str
            Name of the GMSH model.
        '''
        gmsh.initialize()
        self.name = name
        gmsh.model.add(name)
        # self.geo = gmsh.model.geo

    def finalize(self):
        '''
        Finalizes the GMSH model.
        '''
        gmsh.finalize()

    def synchronize(self):
        '''
        Synchronizes the GMSH model.
        this is required before some operations
        '''
        gmsh.model.geo.synchronize()

    def generate_mesh(self, dimension=2):
        ''' 
        Generates the mesh for the GMSH model.'''
        gmsh.model.mesh.generate(dimension)

    def write(self, filename):
        '''
        Writes the GMSH model to a gmsh file.
        '''
        gmsh.write(filename)

    def run_gui(self):
        '''
        Runs the GMSH GUI
        this allow to see triangular mesh results and using the 
        GUI to modify the meshm or to export it to other formats.
        Also allows to see the mesh quality of the elements in a detailed way.
        '''
        gmsh.fltk.run()
        
    def get_triangular_quality(self):
        """
        Get the  follwowing quality measures of the element in the mesh: "minDetJac" and "maxDetJac"
        for the adaptively computed minimal and maximal Jacobian determinant,
        "minSJ" for the sampled minimal scaled jacobien, "minSICN" for the sampled
        minimal signed inverted condition number, "minSIGE" for the sampled signed
        inverted gradient error, "gamma" for the ratio of the inscribed to
        circumcribed sphere radius, "innerRadius" for the inner radius,
        "outerRadius" for the outerRadius, "minIsotropy" for the minimum isotropy
        measure, "angleShape" for the angle shape measure, "minEdge" for the
        minimum straight edge length, "maxEdge" for the maximum straight edge
        length.

        Returns
        -------
        df_qualities : Dataframe
            Datframe with the mesh quality of every triangle element in the mesh.

        """
        
        _, etags, _= gmsh.model.mesh.getElements(dim=2)
        # Get the following quality measures of the element in the mesh
        qualities = {
            "minSICN": gmsh.model.mesh.getElementQualities(etags[0], "minSICN"),  # minimal signed inverted condition number
            "minDetJac": gmsh.model.mesh.getElementQualities(etags[0], "minDetJac"),  # minimal Jacobian determinant
            "maxDetJac": gmsh.model.mesh.getElementQualities(etags[0], "maxDetJac"),  # maximal Jacobian determinant
            "minSJ": gmsh.model.mesh.getElementQualities(etags[0], "minSJ"),  # minimal scaled Jacobian
            "minSIGE": gmsh.model.mesh.getElementQualities(etags[0], "minSIGE"),  # minimal signed inverted gradient error
            "gamma": gmsh.model.mesh.getElementQualities(etags[0], "gamma"),  # ratio of the inscribed to circumcribed sphere radius
            "innerRadius": gmsh.model.mesh.getElementQualities(etags[0], "innerRadius"),  # inner radius
            "outerRadius": gmsh.model.mesh.getElementQualities(etags[0], "outerRadius"),  # outer radius
            "minIsotropy": gmsh.model.mesh.getElementQualities(etags[0], "minIsotropy"),  # minimum isotropy measure
            "angleShape": gmsh.model.mesh.getElementQualities(etags[0], "angleShape"),  # angle shape measure
            "minEdge": gmsh.model.mesh.getElementQualities(etags[0], "minEdge"),  # minimum straight edge length
            "maxEdge": gmsh.model.mesh.getElementQualities(etags[0], "maxEdge")  # maximum straight edge length
        }

        # Save everything on a dataframe
        df_qualities = pd.DataFrame(qualities)
        return df_qualities