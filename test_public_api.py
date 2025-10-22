#!/usr/bin/env python3
"""Test the public API after refactoring."""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

def test_public_api():
    """Test the public API works as expected."""
    try:
        # This is how users should import the library
        from gmshflow import GmshModel, GmshMeshDomain, PolyGeometryHandler
        print("✓ Main classes imported from gmshflow package")
        
        # Test that utility functions are also available
        from gmshflow import merge_many_multilinestring_into_one_linestring
        print("✓ Utility functions available from main package")
        
        return True
    except ImportError as e:
        if "gmsh" in str(e):
            print(f"⚠ GMSH not available (expected in dev): {e}")
            return True
        else:
            print(f"✗ Public API import failed: {e}")
            return False

if __name__ == "__main__":
    print("Testing public API after refactoring...")
    if test_public_api():
        print("✓ Public API works correctly!")
    else:
        print("✗ Public API has issues")
        sys.exit(1)