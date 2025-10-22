"""Unit tests for GmshModel validation without GMSH mocking."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))


class TestGmshModelContextManager:
    """Test GmshModel basic validation without complex GMSH mocking."""

    def test_model_initialization_validation(self):
        """Test basic model initialization."""
        from gmshflow.core.model import GmshModel

        model = GmshModel("test_model")
        assert model.name == "test_model"
        assert not model._initialized

    def test_operations_require_context_validation(self):
        """Test that operations fail outside context manager."""
        from gmshflow.core.model import GmshModel

        model = GmshModel("test")

        # These should fail because model is not initialized
        with pytest.raises(RuntimeError, match="not initialized"):
            model.synchronize()

        with pytest.raises(RuntimeError, match="is not initialized"):
            model.generate_mesh()

        with pytest.raises(RuntimeError, match="is not initialized"):
            model.write("test.msh")

    def test_double_initialization_validation(self):
        """Test double initialization prevention logic."""
        from gmshflow.core.model import GmshModel

        model = GmshModel("test")

        # Simulate initialized state
        model._initialized = True

        with pytest.raises(RuntimeError, match="already initialized"):
            model.__enter__()

    def test_safe_finalization_validation(self):
        """Test that finalize is safe to call multiple times."""
        from gmshflow.core.model import GmshModel

        model = GmshModel("test")

        # Should be safe on uninitialized model (no exceptions)
        model.finalize()
        model.close()

        # Should be safe to call multiple times
        model.finalize()
        model.close()

        # Model should remain uninitialized
        assert not model._initialized
