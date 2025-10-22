"""Test configuration and shared fixtures."""

import sys
from pathlib import Path

import pytest

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

def gmsh_available():
    """Check if GMSH is available."""
    try:
        import gmsh
        gmsh.initialize()
        gmsh.finalize()
        return True
    except (ImportError, Exception):
        return False

# Mark for skipping GMSH tests when not available
gmsh_required = pytest.mark.skipif(
    not gmsh_available(),
    reason="GMSH not available"
)

@pytest.fixture
def gmsh_available_fixture():
    """Fixture to detect GMSH availability."""
    return gmsh_available()
