from math import inf, sqrt
from pickle import FALSE, TRUE
import heapq
import random as rand
import math as m
import bin_lib.map as map
import bin_lib.some_functions as fcts
import numpy as np

class Entity:
    def __init__(self, id: int, pos_street: map.Pos_Street) -> None:
        self._id = id
        self._pos_street = pos_street

    def get_id(self) -> int:
        return self._id

    def move_to(self, pos_street: map.Pos_Street) -> None:
        self._pos_street = pos_street

    def get_pos_street(self) -> map.Pos_Street:
        return self._pos_street

    def get_pos_xy(self) -> tuple[float, float]:
        return self._pos_street.get_pos_xy  
            

class Person(Entity):
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


    def __init__(self, id: int, origin_street: map.Pos_Street, destination_street: map.Pos_Street, map_: map.Map, *args):
        self._id = id
        self._pos_street = origin_street
        self._origin_street = origin_street
        self._destination_street = destination_street
        self._map_ = map_
        self._path = self._path_planner()
        self._fov=None
        self._max_distance_carrying_trash=None
        self._time_of_consumption=None
        self._speed: float=None
        self._has_trash=False    
        self._time_alive=None
        self._distance_carrying_trash=None
        self._distance_to_destination = self._calculate_distance_to_destination()

    def _path_planner(self) -> list[map.Intersection]:
        intersections_info = self._dijkstra()
        return self._path_constructor(intersections_info)

    def _calculate_distance(point1: tuple[float, float], point2: tuple[float, float]) -> float:
        return sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

    def _dijkstra(self):
        intersections_info = [{'distance_to_origin': inf, 'parent': None, 'closed': False} for _ in self._map_.get_intersections_list()]
        pq = []
        starting_intersections = self._origin_street.get_street().get_vector()

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
        destination_intersections = self._destination_street.get_street().get_vector()
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

    def _calculate_distance_to_destination(self):
        start = self.get_pos_xy()
        n_intersections = len(self._path)
        end = self._destination_street.get_pos_xy()
        d = 0
        for i in range(n_intersections):
            aux = self._path[i]
            d += self._calculate_distance(start, aux)
            start = aux
        return d + self._calculate_distance(start, end)

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

    def walk(self, time_walking: int) -> bool:
        """return if the person got to his destination"""
        d = self._speed * time_walking
        start = self._pos_street.get_pos_xy()

        # he arrives?
        if d < self._distance_to_destination:
            self._pos_street = self._destination_street
            return TRUE

        n_intersections = len(self._path)

        # he is already on the destination's street?
        if n_intersections == 0:
            s = self._pos_street.get_street()
            (A, B) = s.get_vector()
            (Ax, Ay) = A.get_pos()
            (Bx, By) = B.get_pos()
            if Ax != Bx: #se nao é reta vertical
                m = (self._destination_street.get_pos_xy()[0] - start[0])*(Bx-Ax)
                if m > 0: #apontam para o mesmo lugar
                    new_pos_in_street = self._pos_street.get_pos_in_street() + d / self._pos_street.get_street().get_length()
                else:#se apontam em sentidos opostos
                    new_pos_in_street = self._pos_street.get_pos_in_street() - d / self._pos_street.get_street().get_length()
            else: #se for vertical
                new_pos_in_street = self._pos_street.get_pos_in_street() + d / (By - Ay)
            self._pos_street = map.Pos_Street(self._pos_street.get_street, new_pos_in_street)
            return FALSE
        
        last_intersection = None
        while d > 0 and n_intersections > 0:
            next_intersection = self._path[0]
            distance_to_next_intersection = self._calculate_distance(start, next_intersection)
            # won't he stop at this street?
            if d > distance_to_next_intersection:
                d -= distance_to_next_intersection
                last_intersection = next_intersection
                self._path.pop(0)
                n_intersections -= 1
            #if he never changed of street - in this walking time
            elif last_intersection == None:
                s = self._pos_street.get_street()
                (A, B) = s.get_vector()
                (Ax, Ay) = A.get_pos()
                (Bx, By) = B.get_pos()
                if Ax != Bx: #se nao é reta vertical
                    m = (self._destination_street.get_pos_xy()[0] - start[0])*(Bx-Ax)
                    if m > 0: #apontam para o mesmo lugar
                        new_pos_in_street = self._pos_street.get_pos_in_street() + d / self._pos_street.get_street().get_length()
                    else:#se apontam em sentidos opostos
                        new_pos_in_street = self._pos_street.get_pos_in_street() - d / self._pos_street.get_street().get_length()
                else: #se for vertical
                    new_pos_in_street = self._pos_street.get_pos_in_street() + d / (By - Ay)
                self._pos_street = map.Pos_Street(self._pos_street.get_street, new_pos_in_street)
                return FALSE
            # if he changed of street
            else:
                s = self._map_.search_for_a_street(last_intersection, next_intersection)
                (A, B) = s.get_vector()
                if last_intersection.get_id() == A.get_id():
                    new_pos_in_street = d / s.get_length()
                else:
                    new_pos_in_street = 1 - d / s.get_length()
                self._pos_street = map.Pos_Street(s, new_pos_in_street)
                return FALSE
        
        # if it got it here, it's because n_intersections got to 0 before d
        # in this case, last_intersection shouldn't be None
        s = self._destination_street.get_street()
        (A, B) = s.get_vector()
        (Ax, Ay) = A.get_pos()
        (Bx, By) = B.get_pos()
        if Ax != Bx: #se nao é reta vertical
            m = (self._destination_street.get_pos_xy()[0] - start[0])*(Bx-Ax)
            if m > 0: #apontam para o mesmo lugar
                new_pos_in_street = self._pos_street.get_pos_in_street() + d / self._pos_street.get_street().get_length()
            else:#se apontam em sentidos opostos
                new_pos_in_street = self._pos_street.get_pos_in_street() - d / self._pos_street.get_street().get_length()
        else: #se for vertical
            new_pos_in_street = self._pos_street.get_pos_in_street() + d / (By - Ay)
        self._pos_street = map.Pos_Street(self._pos_street.get_street, new_pos_in_street)
        return FALSE

                



class Bin(Entity):
    """Private variables:
        id: int
        pos: float[2]
        pos_street: [street,float]
        capacity: float (L)
        filling_rate: float (L/s)
        full: bool
        vol_trash: float (L)
        last_time_was_empted: int (seconds)
    """
    def __init__(self, id: int, pos_street: map.Pos_Street, capacity: float):
        self._id = id
        self._pos_street = pos_street
        self._capacity = capacity
        self._full = FALSE
        self._vol_trash = 0
        self._last_time_was_empted = -1
        self._filling_rate = 0

    def set_full(self):
        self._full = TRUE

    def set_capacity(self, new_capacity: float):
        self._capacity = new_capacity

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
        mu_time_1 (float): gaussian 2 of time (see OBS)
        sigma_time_0 (float): gaussian 1 of time (see OBS)
        sigma_time_1 (float): gaussian 2 of time (see OBS)
        weight_0 (float): weight for the mean of the gaussians of time (see OBS)
        weight_1 (float): weight for the mean of the gaussians of time (see OBS)
        
        p (float): weight0/(weight0 + weight1). In the case if weight0 + weight1 = 0, it will be -1.
    OBS:
        We are considering, here, a function f for the distribution of the time of consumption.
        The time of consumption is the time that the person is consuming his food/whatever, i.e.,
        time after getting out of the store for needing a bin
        This function f will be the mean between 2 gaussians. That way, we can represent both peaks:
        people who eats in place/need a bin just after leaving the store; and people who leaves
        eating and will need a bin a little far from the store.
        Gaussian 0: mean = 0
        Gaussian 1: mean > 0
    """
    def __init__(self, prob_trash: float, payload: dict=None):
        """constructor of Trash_Potential
        Args:
            prob_trash (float): _description_
            payload (dict, optional): dictionary with the following arguments (necessary for prob_trash != 0):
                mu_vol (float): mu of gaussian of volume of the trash
                sigma_vol (float): sigma of gaussian of volume of the trash
                mu_time_1 (float): mu of gaussian 1 of time (see OBS)
                sigmas_time (tuple(float)): sigmas of gaussians of time (see OBS)
                weights (tuple(float)): weights for the mean of the gaussians of time (see OBS)        
        """
        self._prob_trash = prob_trash
        if not prob_trash == 0:
            self._mu_vol = payload['mu_vol']
            self._sigma_vol = payload['sigma_vol']

            sigmas_time = payload['sigmas_time']
            self._sigma_time_0 = sigmas_time[0]
            self._sigma_time_1 = sigmas_time[1]
            self._mu_time_1 = payload['mu_time_1']
            
            weights = payload['weights']
            self._weight_0 = weights[0]
            self._weight_1 = weights[1]
            
            if self._weight_0 + self._weight_1 != 0:
                self._p = self._weight_0 / (self._weight_0 + self._weight_1)
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
            time = m.fabs(rand.gauss(mu=0, sigma=self._sigma_time_0))
            return [True, (vol, time)]
        else:
            time = rand.gauss(mu=self._mu_time_1, sigma=self._sigma_time_1)
            while time <= 0:
                time = rand.gauss(mu=self._mu_time_1, sigma= self._sigma_time_1)
            return [True, (vol, time)]
   
class Commercial_Point(Entity):
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

    def __init__(self, type: int, customer_potential: float, trash_generation_potential: Trash_Potential, pos_street: map.Pos_Street):
        self._type = type
        self._customer_potential = customer_potential
        self._trash_generation_potential = trash_generation_potential
        self._pos_street = pos_street
        self._id = -1

    def get_type(self) -> int:
        return self._type

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


class Entities:
    def __init__(self) -> None:
        self._ppl: list[Person] = []
        self._com_points: list[Commercial_Point] = []
        self._bins: list[Bin] = []
        self._ppl_ids: list[int] = []
        self._bins_ids = []
        self._food_points = 0
        self._nfood_points = 0
        self._job_points = 0
        self._n_people = 0
        self._n_com_points = 0
        self._n_bins = 0
    
    def new_person(self, p: Person)
        self._ppl.append(p)
        self._ppl_ids.append(p.get_id)

    def _id2index(self, id_vec: list[int], id: int):
        return fcts.search_in_vec(id_vec, id)

    def move_bin(self, id_bin: int, new_pos_street: map.Pos_Street):
        self._bins[id_bin].move_to(new_pos_street)

    def make_a_person_walk(self, id_person: int, time_of_walk: int):
        p = self._ppl[fcts.search_in_vec(self._ppl_ids, id_person)]
        path = p.get_path()
        pos = p.get_pos_street()
