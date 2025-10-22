"""Integration tests for GmshModel context manager with real GMSH."""

import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

@pytest.mark.gmsh
def test_model_context_manager_exception_handling():
    """Test context manager behavior with exceptions using real GMSH."""
    gmsh = pytest.importorskip("gmsh")
    from gmshflow.core.model import GmshModel
    
    # Test that GMSH is properly finalized even when exception occurs
    # We'll track the GMSH state before and after
    
    # Ensure GMSH is not initialized before test
    try:
        gmsh.finalize()  # Clean up any previous state
    except:
        pass  # Ignore if already finalized
    
    # Test exception handling with real GMSH
    gmsh_was_initialized = False
    try:
        with GmshModel("exception_test") as model:
            gmsh_was_initialized = True  # If we get here, GMSH was initialized
            assert model._initialized
            
            # Cause an exception
            raise ValueError("Test exception")
            
    except ValueError:
        pass  # Expected exception
    
    # Verify GMSH was properly cleaned up after exception
    assert gmsh_was_initialized, "GMSH should have been initialized"
    
    # Try to verify GMSH is finalized (this is tricky to test directly)
    # We can test by trying to use GMSH without context - should fail
    model2 = GmshModel("test_after_exception")
    with pytest.raises(RuntimeError, match="is not initialized"):
        model2.synchronize()  # Should fail because not in context

@pytest.mark.gmsh
def test_model_context_manager_real_gmsh():
    """Test context manager with real GMSH if available."""
    gmsh = pytest.importorskip("gmsh")
    from gmshflow.core.model import GmshModel
    
    # Test basic context manager lifecycle
    with GmshModel("integration_test") as model:
        assert model._initialized
        assert model.name == "integration_test"
        
        # Test that GMSH operations work
        model.synchronize()  # Should work in context
        
    # After context, should be finalized
    assert not model._initialized