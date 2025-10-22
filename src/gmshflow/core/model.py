"""Core GMSH model wrapper for GMSHFlow."""

import gmsh
import pandas as pd


class GmshModel:
    """GMSH model wrapper with automatic resource management.

    This class provides a high-level interface to GMSH functionality with
    proper resource management via context manager protocol. GMSH resources
    are automatically initialized and cleaned up.

    Args:
        name: Name identifier for the GMSH model.

    Attributes:
        name: The model name used in GMSH.

    Example:
        >>> with GmshModel("groundwater_mesh") as model:
        ...     model.generate_mesh(dimension=2)
        ...     model.write("output.msh")
        ...     # Automatic cleanup when exiting context
    """
    def __init__(self, name: str) -> None:
        """Initialize the GMSH model wrapper.
        
        Args:
            name: Name identifier for the GMSH model.
            
        Note:
            GMSH is not initialized until entering the context manager.
            Use 'with GmshModel(name) as model:' for proper resource management.
        """
        self.name = name
        self._initialized = False

    def __enter__(self) -> "GmshModel":
        """Enter the context manager and initialize GMSH.
        
        Returns:
            The GmshModel instance for use in the context.
            
        Raises:
            RuntimeError: If GMSH is already initialized for this model.
        """
        if self._initialized:
            raise RuntimeError(f"GMSH model '{self.name}' is already initialized")

        gmsh.initialize()
        gmsh.model.add(self.name)
        self._initialized = True
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Exit the context manager and cleanup GMSH resources.
        
        This ensures GMSH is properly finalized even if an exception occurs
        during mesh operations.
        
        Args:
            exc_type: Exception type (if any).
            exc_val: Exception value (if any).
            exc_tb: Exception traceback (if any).
        """
        self.finalize()

    def finalize(self) -> None:
        """Finalize and cleanup the GMSH model safely.
        
        This method is idempotent and safe to call multiple times.
        It will only finalize GMSH if it was previously initialized.
        """
        if self._initialized:
            try:
                gmsh.finalize()
            except Exception as e:
                # Log warning but don't raise - cleanup should be robust
                print(f"Warning: GMSH finalization failed: {e}")
            finally:
                self._initialized = False

    def close(self) -> None:
        """Close the GMSH model (alias for finalize).
        
        Provides a common interface for resource cleanup.
        """
        self.finalize()

    def _ensure_initialized(self) -> None:
        """Ensure GMSH is initialized before operations.
        
        Raises:
            RuntimeError: If GMSH is not initialized (not in context manager).
        """
        if not self._initialized:
            raise RuntimeError(
                f"GMSH model '{self.name}' is not initialized. "
                "Use 'with GmshModel(name) as model:' pattern."
            )

    def synchronize(self) -> None:
        """Synchronize the GMSH model geometry.
        
        This is required before certain mesh operations to ensure
        all geometry changes are properly registered.
        
        Raises:
            RuntimeError: If GMSH is not initialized.
        """
        self._ensure_initialized()
        gmsh.model.geo.synchronize()

    def generate_mesh(self, dimension: int = 2) -> None:
        """Generate mesh for the GMSH model.
        
        Args:
            dimension: Mesh dimension (1=lines, 2=triangles, 3=tetrahedra).
            
        Raises:
            RuntimeError: If GMSH is not initialized.
        """
        self._ensure_initialized()
        gmsh.model.mesh.generate(dimension)

    def write(self, filename: str) -> None:
        """Write the GMSH model to a file.
        
        Args:
            filename: Output filename (typically with .msh extension).
            
        Raises:
            RuntimeError: If GMSH is not initialized.
        """
        self._ensure_initialized()
        gmsh.write(filename)

    def run_gui(self) -> None:
        """Launch the GMSH GUI for interactive mesh visualization.
        
        This allows visualization of triangular mesh results and provides
        GUI tools to modify the mesh or export to other formats. Also
        enables detailed mesh quality analysis.
        
        Raises:
            RuntimeError: If GMSH is not initialized.
        """
        self._ensure_initialized()
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

        Raises:
            RuntimeError: If GMSH is not initialized or no mesh exists.

        Example:
            >>> with GmshModel("test") as model:
            ...     model.generate_mesh()
            ...     quality_df = model.get_triangular_quality()
            ...     poor_elements = quality_df[quality_df['gamma'] < 0.3]
        """
        self._ensure_initialized()

        _, etags, _= gmsh.model.mesh.getElements(dim=2)
        
        # Check if mesh has elements
        if not etags or len(etags) == 0 or len(etags[0]) == 0:
            # Return empty DataFrame with expected columns
            columns = ['minSICN', 'minDetJac', 'maxDetJac', 'minSJ', 'minSIGE', 'gamma', 
                      'innerRadius', 'outerRadius', 'minIsotropy', 'angleShape', 'minEdge', 'maxEdge']
            return pd.DataFrame(columns=columns)
        
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
