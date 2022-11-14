from math import inf, sqrt
from pickle import FALSE, TRUE
import map
import heapq
import random as rand
import math as m

class Person:
    """Private variables
        id: int
        pos: float[2]
        destination: Comertial_
        path: list(Intersection)
        fov: float
        max_distance_carrying_trash: float
        time_of_consumption: int (seconds)
        speed: float (m/s)
        has_trash: bool
            if consuming something, it's false
        time_alive: int (seconds)
        distance_carrying_trash: float

    """


    def __init__(self, id: int, origin_street: map.Pos_Street, destination, map_, *args):
        self._id = id
        self.pos = None
        self._origin = origin
        self._destination = destination
        self._map_ = map_
        self._path = self._path_planner()
        self._fov=None
        self._max_distance_carrying_trash=None
        self._time_of_consumption=None
        self._speed=None
        self._has_trash=False    
        self._time_alive=None
        self._distance_carrying_trash=None
        self._origin_street=None
        self._destination_street=None

    def _path_planner(self) -> list[map.Intersection]:
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

class Bin:
    """Private variables:
        pos: float[2]
        pos_street: [street,float]
        capacity: float (L)
        filling_rate: float (L/s)
        full: bool
        vol_trash: float (L)
        last_time_was_empted: int (seconds)
    """
    def __init__(self, pos: tuple[float, float], capacity: float):
        self._pos = pos
        self._capacity = capacity
        self._full = FALSE
        self._vol_trash = 0
        self._last_time_was_empted = -1
        self._filling_rate = 0

    def set_full(self):
        self._full = TRUE

    def set_capacity(self, new_capacity: float):
        self._capacity = new_capacity
    
    def set_pos(self, new_pos: tuple[float, float]):
        self._pos = new_pos

    def put_trash(self, vol_trash):
        self._vol_trash += vol_trash
        if(self._vol_trash > self._capacity):
            self._full = TRUE
        #update filling_rate

    def empty_bin(self):
        self._vol_trash = 0
        self._full = FALSE
        #self._last_time_was_empted = ?

class Trash_Potential:
    """Variables:
        prob_trash (float): probability of someone leaves the store with trash
        mu_vol (float): gaussian of volume of the trash
        sigma_vol (float): gaussian of volume of the trash
        mu_time_2 (float): gaussian 2 of time (see OBS)
        sigma_time_1 (float): gaussian 1 of time (see OBS)
        sigma_time_2 (float): gaussian 2 of time (see OBS)
        weight_1 (float): weight for the mean of the gaussians of time (see OBS)
        weight_2 (float): weight for the mean of the gaussians of time (see OBS)
        p (float): weight1/(weight1 + weight2). In the case if weight1 + weight2 = 0, it will be -1.

    OBS:
        We are considering, here, a function f for the distribution of the time of consumption.
        The time of consumption is the time that the person is consuming his food/whatever, i.e.,
        time after getting out of the store for needing a bin
        This function f will be the mean between 2 gaussians. That way, we can represent both peaks:
        people who eats in place/need a bin just after leaving the store; and people who leaves
        eating and will need a bin a little far from the store.
    """
    def __init__(self, prob_trash: float, mu_vol: float, sigma_vol: float, mu_time_2: float, sigmas_time: tuple[float], weights_time: tuple[float]):
        self._prob_trash = prob_trash
        self._mu_vol = mu_vol
        self._sigma_vol = sigma_vol
        self._sigma_time_1 = sigmas_time[0]
        self._mu_time_2 = mu_time_2
        self._sigma_time_2 = sigmas_time[1]
        self._weight_1 = weights_time[0]
        self._weight_2 = weights_time[1]
        if self._weight_1 + self._weight_2 != 0:
            self._p = self._weight_1 / (self._weight_1 + self._weight_2)
        else:
            self._p = -1

    def generate_trash(self) -> list[bool, tuple[float]]:
        """it should be called for each person that leaves the store.
        This function see if the person will leave with trash and 
        which are the time of consumption and the volume

        Return:
            bool: if the person leaves with trash
            tuple[float]: time of consumption and volume of trash:
                [0]: time
                [1]: volume
        """
        # trash?
        trash = (rand.random() < self._prob_trash)
        if not trash:
            return False

        # volume of trash
        vol = rand.gauss(mu=self._mu_vol, sigma=self._sigma_vol)
        while vol <= .01e-3: # 0.01 mL = (0.2 mm)³
            vol = rand.gauss(mu=self._mu_vol, sigma=self._sigma_vol)

        # time of consumption
        if self._p == -1:
            return [True, (vol, 0)]
        if rand.random() < self._p:
            time = m.fabs(rand.gauss(mu=0, sigma=self._sigma_time_1))
            return [True, (vol, time)]
        else:
            time = rand.gauss(mu=self._mu_time_2, sigma=self._sigma_time_2)
            while time <= 0:
                time = rand.gauss(mu=self._mu_time_2, sigma= self._sigma_time_2)
            return [True, (vol, time)]
   

class Commercial_Point:
    """Private variables:
        type (int):
            0: Food store
            1: Non-food store
            2: Job
        attractiveness (float):
            0 -> 1
            represents the potential of getting customers
        trash_generation_potential (Trash_Potential)
        pos (float[2])
        pos_street (Pos_Street)
    """

    def __init__(self, type: int, customer_potential: float, trash_generation_potential: Trash_Potential, street: map.Street, pos_in_street: float):
        self._type = type
        self._customer_potential = customer_potential
        self._trash_generation_potential = trash_generation_potential
        self._pos_street = map.Pos_Street(street, pos_in_street)
        A = street.get_vector()[0].get_pos()
        B = street.get_vector()[1].get_pos()
        x = A[0] + (B[0]-A[0])*pos_in_street
        y = A[1] + (B[1]-A[1])*pos_in_street
        self._pos = (x,y)

    def __init__(self, type: int, customer_potential: float, trash_generation_potential: Trash_Potential, pos_street: map.Pos_Street):
        self._type = type
        self._customer_potential = customer_potential
        self._trash_generation_potential = trash_generation_potential
        self._pos_street = pos_street
        A = pos_street.get_street.get_vector()[0].get_pos()
        B = pos_street.get_street.get_vector()[1].get_pos()
        x = A[0] + (B[0]-A[0])*pos_street.get_pos_in_street
        y = A[1] + (B[1]-A[1])*pos_street.get_pos_in_street
        self._pos = (x,y)

    def get_pos(self) -> tuple[float, float]:
        return self._pos_street.get_street

    def get_pos_in_street(self) -> float:
        return self._pos_street.get_pos_in_street

    def generate_trash(self) -> list[bool, tuple[float]]:
        """it should be called for each person that leaves the store.
        This function see if the person will leave with trash and 
        which are the time of consumption and the volume

        Return:
            bool: if the person leaves with trash
            tuple[float]: time of consumption and volume of trash:
                [0]: time
                [1]: volume
        """
        return self._trash_generation_potential.generate_trash()
