from __future__ import annotations
from math import sqrt
import numpy as np
import bin_lib.some_functions as fcs
from bin_lib.consts import *


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
        self._dist_to_neighbors.append(distance if distance != None else self._calculate_distance_to_neighbor(neighbor))

    def get_neighbors(self) -> list[Intersection]:
        return self._neighbors.copy()

    def get_distances_to_neighbors(self) -> list[float]:
        return self._dist_to_neighbors

class Street:
    def __init__(self, A: Intersection, B: Intersection, id: int):
        self._vec = (A, B)
        point1 = A.get_pos()
        point2 = B.get_pos()
        self._length = sqrt((point1[0]-point2[0])**2 + (point1[1]-point2[1])**2)
        A.add_neighbor(B, self._length)
        B.add_neighbor(A, self._length)
        self._intersctions_ids = (A.get_id(), B.get_id())
        self._id = id

    def get_vector(self) -> tuple[Intersection, Intersection]:
        return self._vec

    def get_length(self) -> float:
        return self._length

    def get_intersections_ids(self) -> tuple[int, int]:
        return self._intersctions_ids

    def get_id(self) -> int:
        return self._id

class Map:
    def __init__(self):
        self._intersections: list[Intersection] = []
        self._streets: list[Street] = []

    def add_intersection(self, x: float, y: float):
        self._intersections.append(Intersection(len(self._intersections),x,y))

    def add_street(self, indexA: int, indexB: int):
        new_street = Street(self._intersections[indexA], self._intersections[indexB], len(self._streets))
        self._streets.append(new_street)

    def get_intersections_list(self) -> list[Intersection]:
        """CAREFUL: the elements in the list are the original ones"""
        return self._intersections.copy()
    
    def get_intersection(self, index: int) -> Intersection:
        """CAREFUL: the return isn't a copy; it's the original Intersection"""
        return self._intersections[index]

    def get_streets_list(self) -> list[Street]:
        """CAREFUL: the elements in the list are the original ones"""
        return self._streets.copy()

    def get_street(self, index: int) -> Street:
        """CAREFUL: the return isn't a copy; it's the original Street"""
        return self._streets[index]

    def search_for_a_street(self, intersection_A: Intersection, intersection_B: Intersection) -> Street:
        A = intersection_A.get_id()
        B = intersection_B.get_id()
        if A == B:
            return None
        for s in self._streets:
            (a,b) = s.get_intersections_ids()
            if (a == A or a == B) and (b == A or b == B):
                return s
        return None

            


if __name__ == "__main__":
    mapa = Map()
    positions = [(0.0, 0.0), (0.0, 1.0)]
    for position in positions:
        mapa.add_intersection(*position)

    mapa.add_street(0, 1)

    print(mapa.get_intersection(1).get_id())
    list(map(lambda x: print(x.get_vector()[0].get_pos(),"-", x.get_vector()[1].get_pos()), mapa.get_streets_list()))

class Pos_Street:
    def __init__(self, street: Street, pos_in_street: float):
        self._street: Street = street
        self._pos_in_street: float = pos_in_street
        (A, B) = street.get_vector()
        A = A.get_pos()
        B = B.get_pos()
        x = A[0]*(1-pos_in_street) + B[0]*pos_in_street
        y = A[1]*(1-pos_in_street) + B[1]*pos_in_street
        self._pos_xy = (x,y)
        
    def get_street(self) -> Street:
        return self._street

    def get_pos_in_street(self) -> float:
        return self._pos_in_street

    def get_pos_xy(self) -> tuple[float, float]:
        return self._pos_xy

def read_map(file_intersections: str, file_streets: str, file_commertial_points: str) -> Map:
    """ function used to read a file with intersections and
    return the map (Map)

    Args:
        file_intersections (str):
            - name of the file with the intersections
            - must be one intersection per line, with
            x and y separated by a comma: "x,y"
        file_streets (str):
            - name of the file with the streets
            - must be one street per line, with A and B
            being the index of the intersections and separeted by
            comma: "A,B"
        file_commertial_points (str):
            - name of the file with the commertial points
            - must have the type (int) and the position
            - the position must be like this: "index,float"
                - index = index of the street
                - float = position in the street (0 -> 1)
            - so: "type,index,float"
            -idk if we will use this ??????????????????
        Return:
        map (Map)
    """

    mapa = Map()

    f = open(file_intersections, "r")
    for line in f:
        point = line.split(",")
        mapa.add_intersection(float(point[0]), float(point[1]))    
    f.close()

    f = open(file_streets, "r")
    for line in f:
        street = line.split(",")
        mapa.add_street(int(street[0]), int(street[1]))
    f.close()

    return mapa

def validate_a_map(mapa: Map) -> bool:
    intersections = mapa.get_intersections_list()
    streets = mapa.get_streets_list()

    # see if streets cross one another
    # and if neighbours are good
    Ls = len(streets)
    for i in range(Ls):
        (A0, B0) = streets[i].get_vector()
        good = False
        for possible_b0 in A0.get_neighbors():
            if possible_b0.get_id() == B0.get_id():
                good = True
        if not good:
            print("problem with neighbours")
            return False
        good = False
        for possible_a0 in B0.get_neighbors():
            if possible_a0.get_id() == A0.get_id():
                good = True
        if not good:
            print("problem with neighbours")
            return False

        A0 = A0.get_pos()
        B0 = B0.get_pos()
        (A0x,A0y) = A0
        (B0x,B0y) = B0
        (m0,n0,r0) = fcs.calculate_line_equation(A0,B0)
        for k in range(i+1, Ls):
            (A, B) = streets[k].get_vector()
            A = A.get_pos()
            B = B.get_pos()
            (Ax,Ay) = A
            (Bx,By) = B
            
            # how many common intersections in common:
            # 0 -> idk
            # 1 -> has no problem
            # 2 -> same street -> problem
            A0_A = fcs.calculate_distance(A0,A)
            B0_A = fcs.calculate_distance(B0,A)
            A0_B = fcs.calculate_distance(A0,B)
            B0_B = fcs.calculate_distance(B0,B)
            # 2 in common
            if abs(A0_A * A0_B) < EPSILON and abs(B0_A * B0_B) < EPSILON:
                return False
            # 0 in common
            if abs(A0_A * A0_B) > EPSILON and abs(B0_A * B0_B) > EPSILON:
                (m,n,r) = fcs.calculate_line_equation(A,B)
                A0_in_mnr = m*A0x + n*A0y - r
                B0_in_mnr = m*B0x + n*B0y - r
                A_in_m0n0r0 = m0*Ax + n0*Ay - r0
                B_in_m0n0r0 = m0*Bx + n0*By - r0
                # is there any point from the streets on the line of the other street?
                if A0_in_mnr * B0_in_mnr < EPSILON or A_in_m0n0r0 * B_in_m0n0r0 < EPSILON:
                    # crossing:
                    if A0_in_mnr * B0_in_mnr < -EPSILON and A_in_m0n0r0 * B_in_m0n0r0 < -EPSILON:
                        print(f'({i}, {k}) cross')
                        print(f'{i} = ({streets[i].get_vector()[0].get_pos()};{streets[i].get_vector()[1].get_pos()})')
                        print(f'{k} = ({streets[k].get_vector()[0].get_pos()};{streets[k].get_vector()[1].get_pos()})')
                        return False
                    # if it didn't cross, then the problem is if A0 inside AB, or A inside A0B0, or something like that
                    # otherwise, there is no problem with crossing or something like that
                    A0_B0 = fcs.calculate_distance(A0,B0)
                    A_B = fcs.calculate_distance(A,B)
                    if A0_A + A0_B <= A_B + EPSILON or B0_A + B0_B <= A_B + EPSILON or \
                        A0_A + B0_A <= A0_B0 + EPSILON or A0_B + B0_B <= A0_B0 + EPSILON:
                        if A0_A + A0_B <= A_B + EPSILON:
                            print(f'({A0} inside ({A};{B})')
                        if B0_A + B0_B <= A_B + EPSILON:
                            print(f'({B0} inside ({A};{B})')
                        if A0_A + B0_A <= A0_B0 + EPSILON:
                            print(f'({A} inside ({A0};{B0})')
                        if A0_B + B0_B <= A0_B0 + EPSILON:
                            print(f'({B} inside ({A0};{B0})')
                        return False

    return True