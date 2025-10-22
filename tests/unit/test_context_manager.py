"""Unit tests for GmshModel context manager (mocked)."""

import sys
from pathlib import Path
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))

from mock_gmsh import MockGmsh


class TestGmshModelContextManager:
    """Test GmshModel context manager with mocked GMSH."""

    @patch('gmshflow.core.model.gmsh', MockGmsh())
    def test_context_manager_lifecycle(self):
        """Test context manager initialization and cleanup."""
        from gmshflow.core.model import GmshModel

        model = GmshModel("test")
        assert not model._initialized

        with model as m:
            assert m is model
            assert model._initialized

        assert not model._initialized

    @patch('gmshflow.core.model.gmsh', MockGmsh())
    def test_operations_require_context(self):
        """Test that operations fail outside context manager."""
        from gmshflow.core.model import GmshModel

        model = GmshModel("test")

        with pytest.raises(RuntimeError, match="not initialized"):
            model.synchronize()

        with pytest.raises(RuntimeError, match="not initialized"):
            model.generate_mesh()

    @patch('gmshflow.core.model.gmsh', MockGmsh())
    def test_double_initialization_prevented(self):
        """Test that double initialization is prevented."""
        from gmshflow.core.model import GmshModel

        model = GmshModel("test")

        with model:
            with pytest.raises(RuntimeError, match="already initialized"):
                with model:
                    pass

    @patch('gmshflow.core.model.gmsh', MockGmsh())
    def test_safe_finalization(self):
        """Test that finalize is safe to call multiple times."""
        from gmshflow.core.model import GmshModel

        model = GmshModel("test")

        # Safe on uninitialized model
        model.finalize()
        model.close()

        with model:
            pass

        # Safe after context exit
        model.finalize()
        model.close()
