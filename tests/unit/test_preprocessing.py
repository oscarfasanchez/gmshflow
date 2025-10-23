"""Unit tests for preprocessing utilities."""

import sys
from pathlib import Path

import geopandas as gpd
import pytest
from shapely.geometry import LineString, MultiLineString

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

# Import directly from preprocessing to avoid GMSH import chain
from gmshflow.utils.preprocessing import (
    merge_many_multilinestring_into_one_linestring,
    simplify_keeping_topology,
)


class TestMergeMultiLineString:
    """Test multilinestring merging function."""

    def test_merge_multilinestring_basic(self):
        """Test basic multilinestring merging."""
        # Create test multilinestring
        line1 = LineString([(0, 0), (1, 1)])
        line2 = LineString([(1, 1), (2, 2)])
        multiline = MultiLineString([line1, line2])

        gdf = gpd.GeoDataFrame({'id': [1]}, geometry=[multiline])

        result = merge_many_multilinestring_into_one_linestring(gdf)

        assert all(result.geom_type == 'LineString')
        assert len(result) == 1

    def test_no_multilinestring_raises_error(self):
        """Test that function raises error when no multilinestrings present."""
        line = LineString([(0, 0), (1, 1)])
        gdf = gpd.GeoDataFrame({'id': [1]}, geometry=[line])

        with pytest.raises(ValueError, match="GeoDataFrame must contain MultiLineString geometries"):
            merge_many_multilinestring_into_one_linestring(gdf)

class TestSimplifyKeepingTopology:
    """Test topology-preserving simplification."""

    def test_simplify_basic(self):
        """Test basic simplification."""
        # Create complex polygon that can be simplified
        coords = [(0, 0), (0.1, 0.05), (1, 0), (1, 1), (0, 1), (0, 0)]
        from shapely.geometry import Polygon
        poly = Polygon(coords)
        gdf = gpd.GeoDataFrame({'id': [1]}, geometry=[poly], crs='EPSG:4326')

        result = simplify_keeping_topology(gdf, cs=0.2)

        # Should return simplified geometry
        assert len(result) == 1
        assert result.iloc[0].geometry.is_valid

    def test_simplify_preserves_dataframe_structure(self):
        """Test that simplification preserves dataframe structure."""
        coords = [(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]
        from shapely.geometry import Polygon
        poly = Polygon(coords)
        gdf = gpd.GeoDataFrame({'name': ['test'], 'value': [42]}, geometry=[poly], crs='EPSG:4326')

        result = simplify_keeping_topology(gdf, cs=0.1)

        assert list(result.columns) == ['name', 'value', 'geometry']
        assert result['name'].iloc[0] == 'test'
        assert result['value'].iloc[0] == 42
