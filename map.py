from __future__ import annotations
from math import sqrt
import numpy as np

class Intersection:
    def __init__(self, id: int, x: float, y: float):
        self._id = id
        self._pos = (x,y)
        self._neighbors: list[Intersection] = []
        self._dist_to_neighbors: list[float] = []

    def get_id(self) -> int:
        return self._id

    def get_pos(self) -> tuple[float, float]:
        return self._pos
    
    def _calculate_distance_to_neighbor(self, neighbor: Intersection) -> float:
        point1 = self.get_pos()
        point2 = neighbor.get_pos()

        return sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

    def add_neighbor(self, neighbor: Intersection, distance: float = None) -> None:
        self._neighbors.append(neighbor)
        self._dist_to_neighbors = distance if distance != None else self._calculate_distance_to_neighbor(neighbor)

    def get_neighbors(self) -> list[Intersection]:
        return self._neighbors

    def get_distances_to_neighbors(self) -> list[float]:
        return self._dist_to_neighbors

class Street:
    def __init__(self, A: Intersection, B: Intersection):
        self._vec = (A, B)
        point1 = A.get_pos()
        point2 = B.get_pos()
        self._length = sqrt((point1[0]-point2[0])**2 + (point1[1]-point2[1])**2)
        A.add_neighbor(B, self._length)
        B.add_neighbor(A, self._length)

    def get_vector(self) -> tuple(Intersection):
        return self._vec

    def get_length(self) -> float:
        return self._length

class Map:
    def __init__(self):
        self._intersections = []
        self._streets = []

    def add_intersection(self, x: float, y: float):
        self._intersections.append(Intersection(len(self._intersections),x,y))

    def add_street(self, indexA: int, indexB: int):
        new_street = Street(self._intersections[indexA], self._intersections[indexB])
        self._streets.append(new_street)

    def get_intersection_list(self) -> list[Intersection]:
        return self._intersections
    
    def get_intersection(self, index: int) -> Intersection:
        return self._intersections[index]

    def get_streets_list(self) -> list[Street]:
        return self._streets



if __name__ == "__main__":
    mapa = Map()
    positions = [(0.0, 0.0), (0.0, 1.0)]
    for position in positions:
        mapa.add_intersection(*position)

    mapa.add_street(0, 1)

    print(mapa.get_intersection(1).get_id())
    list(map(lambda x: print(x.get_vector()[0].get_pos(),"-", x.get_vector()[1].get_pos()), mapa.get_streets_list()))
