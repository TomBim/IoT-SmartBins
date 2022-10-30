from math import sqrt
import numpy as np

class Intersection:
    _pos: float[2] #tupla position

    def __init__(self, x: float, y: float):
        self._pos = (x,y)

    def getPos(self) -> float[2]:
        return self._pos.copy()

class Street:
    _vec: Intersection[2]
    _length: float

    def __init__(self, A: Intersection, B: Intersection):
        self._vec = (A, B)
        point1 = A.getPos()
        point2 = B.getPos()
        self._length = sqrt((point1[0]-point2[0])**2 + (point1[1]-point2[1])**2)

    def getVector(self) -> Intersection[2]:
        return self._vec.copy()

    def getLength(self) -> float:
        return self._length

class Map:
    """
    """
    _intersections: list(Intersection)
    _streets: list(Street)
    _matrix: np.array #adjance matrix

    def __init__(self):
        self._intersections = None
        self._streets = None

    def addIntersection(self, x: float, y: float):
        self._intersections.append(Intersection(x,y))

    def addStreet(self, indexA: int, indexB: int):
        newStreet = Street(self._intersections[indexA], self._intersections[indexB])
        self._streets.append(newStreet)
        self._matrix

