"""Core GMSH model wrapper for GMSHFlow."""

from typing import Optional
import gmsh
import pandas as pd


class GmshModel:
    """GMSH model wrapper for mesh generation and management.

    This class provides a high-level interface to GMSH functionality,
    handling model initialization, mesh generation, and output operations.

    Args:
        name: Name identifier for the GMSH model.

    Attributes:
        name: The model name used in GMSH.

    Example:
        >>> model = GmshModel("groundwater_mesh")
        >>> model.generate_mesh(dimension=2)
        >>> model.write("output.msh")
        >>> model.finalize()
    """
    def __init__(self, name: str) -> None:
        """Initialize the GMSH model.
        
        Args:
            name: Name identifier for the GMSH model.
        """
        gmsh.initialize()
        self.name = name
        gmsh.model.add(name)
        # self.geo = gmsh.model.geo

    def finalize(self) -> None:
        """Finalize and cleanup the GMSH model.
        
        This should be called when done with the model to free resources.
        """
        gmsh.finalize()

    def synchronize(self) -> None:
        """Synchronize the GMSH model geometry.
        
        This is required before certain mesh operations to ensure
        all geometry changes are properly registered.
        """
        gmsh.model.geo.synchronize()

    def generate_mesh(self, dimension: int = 2) -> None:
        """Generate mesh for the GMSH model.
        
        Args:
            dimension: Mesh dimension (1=lines, 2=triangles, 3=tetrahedra).
        """
        gmsh.model.mesh.generate(dimension)

    def write(self, filename: str) -> None:
        """Write the GMSH model to a file.
        
        Args:
            filename: Output filename (typically with .msh extension).
        """
        gmsh.write(filename)

    def run_gui(self) -> None:
        """Launch the GMSH GUI for interactive mesh visualization.
        
        This allows visualization of triangular mesh results and provides
        GUI tools to modify the mesh or export to other formats. Also
        enables detailed mesh quality analysis.
        """
        gmsh.fltk.run()
        
    def get_triangular_quality(self) -> pd.DataFrame:
        """Get comprehensive quality measures for all triangular mesh elements.

        Computes various geometric quality metrics for each triangle element
        in the mesh to assess mesh quality and identify problematic elements.

        Returns:
            DataFrame with quality measures for every triangle element. Columns include:
            - minDetJac/maxDetJac: Minimal/maximal Jacobian determinant
            - minSJ: Sampled minimal scaled Jacobian  
            - minSICN: Minimal signed inverted condition number
            - minSIGE: Minimal signed inverted gradient error
            - gamma: Ratio of inscribed to circumscribed sphere radius
            - innerRadius/outerRadius: Inner and outer radius measures
            - minIsotropy: Minimum isotropy measure
            - angleShape: Angle shape measure
            - minEdge/maxEdge: Minimum and maximum straight edge lengths

        Example:
            >>> model = GmshModel("test")
            >>> model.generate_mesh()
            >>> quality_df = model.get_triangular_quality()
            >>> poor_elements = quality_df[quality_df['gamma'] < 0.3]
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