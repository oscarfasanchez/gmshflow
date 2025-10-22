"""Mock GMSH module for testing without GMSH dependency."""

class MockGmshModel:
    def add(self, name):
        pass

    def getEntities(self, dim):
        return [(1, 1), (1, 2)]  # Mock entity data

    class geo:
        @staticmethod
        def addPoint(x, y, z, size):
            return 1

        @staticmethod
        def addLine(p1, p2):
            return 1

        @staticmethod
        def addCurveLoop(lines):
            return 1

        @staticmethod
        def addPlaneSurface(loops):
            return 1

        @staticmethod
        def synchronize():
            pass

    class mesh:
        @staticmethod
        def generate(dim):
            pass

        @staticmethod
        def getElements(dim):
            return ([], [[1, 2, 3]], [])  # Mock mesh data

        @staticmethod
        def getElementQualities(elements, quality_type):
            return [0.8, 0.9, 0.7]  # Mock quality data

class MockGmsh:
    model = MockGmshModel()

    @staticmethod
    def initialize():
        pass

    @staticmethod
    def finalize():
        pass

    @staticmethod
    def write(filename):
        pass

    @staticmethod
    def is_initialized():
        return True

    class fltk:
        @staticmethod
        def run():
            pass
