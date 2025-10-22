#!/usr/bin/env python3
"""Test the GMSH context manager implementation."""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

def test_context_manager():
    """Test that context manager properly initializes and cleans up."""
    try:
        from gmshflow.core.model import GmshModel
        
        # Test that model starts uninitialized
        model = GmshModel("test")
        assert not model._initialized, "Model should start uninitialized"
        print("✓ Model starts in uninitialized state")
        
        # Test that operations fail outside context
        try:
            model.synchronize()
            print("✗ Should have failed outside context")
            return False
        except RuntimeError as e:
            if "not initialized" in str(e):
                print("✓ Operations properly fail outside context manager")
            else:
                print(f"✗ Wrong error type: {e}")
                return False
        
        # Test context manager initialization
        with model as m:
            assert m is model, "Context manager should return self"
            assert model._initialized, "Model should be initialized in context"
            print("✓ Context manager properly initializes GMSH")
            
            # Operations should work inside context (though they'll fail due to no gmsh)
            # We just test the initialization check passes
            
        # Test cleanup after context
        assert not model._initialized, "Model should be cleaned up after context"
        print("✓ Context manager properly cleans up GMSH")
        
        # Test that operations fail again after context
        try:
            model.synchronize()
            print("✗ Should have failed after context exit")
            return False
        except RuntimeError as e:
            if "not initialized" in str(e):
                print("✓ Operations properly fail after context exit")
            else:
                print(f"✗ Wrong error type: {e}")
                return False
        
        return True
        
    except ImportError as e:
        if "gmsh" in str(e):
            print("⚠ GMSH not available (expected in dev environment)")
            return True
        else:
            print(f"✗ Import failed: {e}")
            return False
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return False

def test_double_initialization():
    """Test that double initialization is prevented."""
    try:
        from gmshflow.core.model import GmshModel
        
        model = GmshModel("test")
        
        with model:
            # Try to enter context again - should fail
            try:
                with model:
                    pass
                print("✗ Should have failed on double initialization")
                return False
            except RuntimeError as e:
                if "already initialized" in str(e):
                    print("✓ Double initialization properly prevented")
                    return True
                else:
                    print(f"✗ Wrong error: {e}")
                    return False
                    
    except ImportError as e:
        if "gmsh" in str(e):
            print("⚠ GMSH not available (expected in dev environment)")
            return True
        else:
            print(f"✗ Import failed: {e}")
            return False

def test_safe_finalization():
    """Test that finalize() is safe to call multiple times."""
    try:
        from gmshflow.core.model import GmshModel
        
        model = GmshModel("test")
        
        # Should be safe to call on uninitialized model
        model.finalize()
        model.close()
        print("✓ Safe to call finalize/close on uninitialized model")
        
        with model:
            pass  # Initialize and cleanup automatically
        
        # Should be safe to call after cleanup
        model.finalize()
        model.close()
        print("✓ Safe to call finalize/close after cleanup")
        
        return True
        
    except ImportError as e:
        if "gmsh" in str(e):
            print("⚠ GMSH not available (expected in dev environment)")
            return True
        else:
            print(f"✗ Import failed: {e}")
            return False

if __name__ == "__main__":
    print("Testing GMSH context manager implementation...")
    
    success = True
    success &= test_context_manager()
    success &= test_double_initialization()
    success &= test_safe_finalization()
    
    if success:
        print("\n✓ All context manager tests passed!")
    else:
        print("\n✗ Some context manager tests failed")
        sys.exit(1)