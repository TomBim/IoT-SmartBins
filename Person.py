from math import inf, sqrt
import Map
import heapq

from sklearn import neighbors
class Person:
    _id=None
    _pos=None
    _origin=None
    _origin_street=None
    _destination=None
    _destination_street=None
    _path=None #list of intersection points
    _fov=None
    _max_distance_carrying_trash=None #maximal distance that this person carries a trash
    _time_of_consumption=None #time needed for consuming
    _speed=None
    _has_trash=False #if consuming something, it's false
    _time_alive=None
    _distance_carrying_trash=None
    _map_=None

    def __init__(self, id, origin, destination, map_, *args):
        self._id = id
        self._origin = origin
        self._destination = destination
        self._map_ = map_
        self._path = self._path_planner()

    def _path_planner(self) -> list[Map.Intersection]:
        intersections_info = self._dijkstra()
        return self._path_constructor(intersections_info)

    def _calculate_distance(point1: tuple[float, float], point2: tuple[float, float]) -> float:
        return sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

    def _dijkstra(self):
        intersections_info = [{'distance_to_origin': inf, 'parent': None, 'closed': False} for _ in self._map_.get_intersection_list()]
        pq = []
        starting_intersections = self._origin_street.get_vec()

        for intersection in starting_intersections:
            distance = self._calculate_distance(self._origin, intersection.get_pos())
            intersections_info[intersection.get_id()]['distance_to_origin'] = distance
            heapq.heappush(pq, (distance, intersection))

        while len(pq) > 0:
            (distance, intersection) = heapq.heappop(pq)
            id = intersection.get_id()
            if not intersections_info[id]['closed']:
                intersections_info[id]['closed'] = True
                neighbors = intersection.get_neighbors()
                neighbors_dists = intersection.get_distances_to_neighbors()
                for neighbor, neighbor_dist in zip(neighbors, neighbors_dists):
                    neighbor_id = neighbor.get_id()
                    if intersections_info[neighbor_id]['distance_to_origin'] > distance + neighbor_dist:
                        neighbor_dist = distance + neighbor_dist
                        intersections_info[neighbor_id]['distance_to_origin'] = neighbor_dist
                        intersections_info[neighbor_id]['parent'] = id
                        heapq.heappush(pq, (neighbor_dist, neighbor))
        return intersections_info

    def _path_constructor(self, intersections_info):
        destination_intersections = self._destination_street.get_vec()
        intersection1_id = destination_intersections[0].get_id()
        intersection2_id = destination_intersections[1].get_id()
        dist_destination_intersection1 = self._calculate_distance(destination_intersections[0].get_pos(), self._destination)
        dist_destination_intersection2 = self._calculate_distance(destination_intersections[1].get_pos(), self._destination)
        dist_origin_intersection1 = intersections_info[intersection1_id]['distance_to_origin']
        dist_origin_intersection2 = intersections_info[intersection2_id]['distance_to_origin']
        last_intersection_id = intersection1_id \
            if dist_destination_intersection1 + dist_origin_intersection1 < dist_destination_intersection2 + dist_origin_intersection2 \
            else intersection2_id

        reversed_path = []
        intersection_id = last_intersection_id
        while intersections_info[intersection_id]['parent'] is not None:
            reversed_path.append(self._map_.get_intersection(intersection_id))
            intersection_id = intersections_info[intersection_id]['parent']
        
        return reversed_path[::-1]

    def get_id(self):
        return self._id

    def get_path(self):
        return self._path

    def get_origin(self):
        return self._origin

    def get_destination(self):
        return self._destination
    
    def set_origin(self, origin):
        self._origin = origin

    def set_destination(self, destination):
        self._destination = destination
