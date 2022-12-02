from __future__ import annotations
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
        return self._pos_street.get_pos_xy()
            

class Person(Entity):
    """Private variables
        id(int)
        pos_street (Pos_Street)
        destination_street (Pos_Street): 
        destination (Commertial_Point)
        path: list(Intersection)
        fov (float): vision distance
        max_time_carrying_trash: float
        time_of_consumption: int (seconds)
        speed: float (m/s)
        has_trash: bool
            if consuming something, it's false
        time_alive: int (seconds)
        time_carrying_trash: float

    """


    def __init__(self, id: int, origin: Commercial_Point, destination: Commercial_Point, map_: map.Map, payload: dict,*args):
        super().__init__(id, origin.get_pos_street())
        self._origin = origin
        self._origin_street = origin.get_pos_street()
        self._destination = destination
        self._destination_street = destination.get_pos_street()
        self._map_ = map_
        self._path = self._path_planner()
        self._fov: float= payload.get("fov")
        self._max_time_carrying_trash: int = payload.get("max_time_carrying_trash")
        self._trash_volume = payload["trash_volume"]
        self._time_of_consumption = payload["time_of_consumption"]
        self._has_trash: bool = (self._time_of_consumption == 0)
        self._speed: float = payload.get("speed")
        self._time_alive = 0
        self._time_carrying_trash = 0
        self._distance_to_destination = self._calculate_distance_to_destination()

    def _path_planner(self) -> list[map.Intersection]:
        intersections_info = self._dijkstra()
        return self._path_constructor(intersections_info)



    def _dijkstra(self):
        intersections_info = [{'distance_to_origin': inf, 'parent': None, 'closed': False} for _ in self._map_.get_intersections_list()]
        pq = []
        starting_intersections = self._origin_street.get_street().get_vector()

        for intersection in starting_intersections:
            distance = fcts.calculate_distance(self._origin.get_pos_xy(), intersection.get_pos())
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
        dist_destination_intersection1 = fcts.calculate_distance(destination_intersections[0].get_pos(), self._destination.get_pos_xy())
        dist_destination_intersection2 = fcts.calculate_distance(destination_intersections[1].get_pos(), self._destination.get_pos_xy())
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
            aux = self._path[i].get_pos()
            d += fcts.calculate_distance(start, aux)
            start = aux
        return d + fcts.calculate_distance(start, end)

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

    def get_time_to_finish_consumption(self):
        return self._time_of_consumption - self._time_alive

    def get_max_time_of_carrying_trash(self):
        return self._max_time_carrying_trash

    def walk(self, time_walking: int) -> bool:
        """return if the person got to his destination"""
        d = self._speed * time_walking
        start = self._pos_street.get_pos_xy()

        # he arrives?
        if d > self._distance_to_destination:
            self._pos_street = self._destination_street
            return TRUE

        n_intersections = len(self._path)
        self._time_alive += time_walking

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
            self._pos_street = map.Pos_Street(self._pos_street.get_street(), new_pos_in_street)
            return FALSE
        
        last_intersection = None
        while d > 0 and n_intersections > 0:
            next_intersection = self._path[0]
            distance_to_next_intersection = fcts.calculate_distance(start, next_intersection.get_pos())
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
                self._pos_street = map.Pos_Street(self._pos_street.get_street(), new_pos_in_street)
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
        self._pos_street = map.Pos_Street(self._pos_street.get_street(), new_pos_in_street)
        return FALSE

    def get_fov(self):
        return self._fov

    def get_trash_volume(self):
        return self._trash_volume

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
        super().__init__(id, pos_street)
        self._capacity = capacity
        self._full = FALSE
        self._vol_trash = 0
        self._last_time_was_empted = -1
        self._filling_rate = 0

    def is_full(self):
        return self._full

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
            return [False, (0,0)]

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
        super().__init__(-1, pos_street)
        self._type = type
        self._customer_potential = customer_potential
        self._trash_generation_potential = trash_generation_potential
        self._attractiveness = rand.random() + 1e-6 #cant be zero

    def get_attractiveness(self) -> float:
        return self._attractiveness

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

class Everything:
    def __init__(self, mapa: map.Map) -> None:
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
        self._last_persons_id = -1
        self._last_bins_id = -1
        self._com_points_attractiveness = np.array([])
        self._trash_in_the_strets = 0
        self._map = mapa
    
    def new_person(self):
        origin = rand.choices(range(self._food_points+self._nfood_points), self._com_points_attractiveness[0:(self._food_points + self._nfood_points)])
        origin = origin[0]
        destination = rand.choices(range(self._food_points+self._nfood_points+self._job_points), weights=self._com_points_attractiveness)
        while destination == origin:
            destination = rand.choices(range(self._food_points+self._nfood_points+self._job_points), weights=self._com_points_attractiveness)
        destination = destination[0]
        payload = {"origin": self._com_points[origin],
            "destination": self._com_points[destination]
        }

        [trash_bool, trash_info] = self._com_points[origin].generate_trash()
        if trash_bool:
            auxfov = abs(rand.gauss(mu=8, sigma=2))
            payload["fov"] = auxfov # between 3 and 13 meters
            aux = rand.gauss(mu=32*auxfov/8, sigma=18)
            if aux < 0:
                payload["max_time_carrying_trash"] = 0
            else:
                payload["max_time_carrying_trash"] = rand.gauss(mu=38, sigma=18)
            payload["speed"] = abs(rand.gauss(mu = 1.25, sigma = .25))
            payload["trash_volume"] = trash_info[0]
            payload["time_of_consumption"] = trash_info[1]

            id = self._next_unused_person_id()
            p = Person(id, payload["origin"], payload["destination"], self._map, payload)
            index = self._id2index(self._ppl_ids, id)
            self._ppl.insert(index, p)
            self._ppl_ids.insert(index, id)
            self._last_persons_id += 1

    def _next_unused_person_id(self):
        if self._last_persons_id > 1e9:
            self._last_persons_id = -1
            return 0
        return self._last_persons_id + 1

    def _next_unused_bin_id(self):
        if self._last_bins_id > 1e9:
            self._last_bins_id = -1
            return 0
        return self._last_bins_id + 1

    def _id2index(self, id_vec: list[int], id: int):
        return fcts.search_in_vec(id_vec, id)

    def remove_a_person(self, id: int):
        index = self._id2index(self._ppl_ids, id)
        self._ppl.pop(index)
        self._ppl_ids.pop(index)

    def move_bin(self, id_bin: int, new_pos_street: map.Pos_Street):
        index = self._id2index(self._bins_ids, id_bin)
        self._bins[index].move_to(new_pos_street)

    def remove_bin(self, id_bin: int):
        index = self._id2index(self._bins_ids, id_bin)
        self._bins[index].pop(index)
        self._bins_ids[index].pop(index)
    
    def new_bin(self, capacity, pos_street: map.Pos_Street):
        id = self._next_unused_bin_id()
        index = self._id2index(self._bins_ids, id)
        self._bins_ids.insert(index, id)
        self._bins.insert(index, Bin(id, pos_street, capacity))
        self._last_bins_id += 1
        self._n_bins += 1

    def make_a_person_walk(self, id_person: int, time_of_walk: int):
        p = self._ppl[fcts.search_in_vec(self._ppl_ids, id_person)]
        p.walk(time_of_walk)

    def new_com_point(self, type: int, customer_potential: float, trash_generation_potential: Trash_Potential, pos_street: map.Pos_Street):
        c = Commercial_Point(type, customer_potential, trash_generation_potential, pos_street)
        if type == 0:
            self._com_points.insert(self._food_points, c)
            self._com_points_attractiveness = np.insert(self._com_points_attractiveness, self._food_points, c.get_attractiveness())
            self._food_points += 1
        elif type == 1:
            self._com_points.insert(self._food_points+self._nfood_points, c)
            self._com_points_attractiveness = np.insert(self._com_points_attractiveness, self._food_points+self._nfood_points, c.get_attractiveness())
            self._nfood_points += 1
        else:
            self._job_points += 1
            self._com_points.append(c)
            self._com_points_attractiveness = np.append(self._com_points_attractiveness, c.get_attractiveness())

    def get_com_points(self):
        return self._com_points

    def get_com_points_attractiveness(self):
        return self._com_points_attractiveness

    def check_for_nearby_bins(self, id_person: int) -> tuple[bool, int]: 
        """
        This function checks if there is any bin closer than their fov
        If there is a bin nearby, it should return TRUE 
        Else, it should return FALSE
        """
        p = self._ppl[fcts.search_in_vec(self._ppl_ids, id_person)]
        fov = p.get_fov()

        for b in self._bins:
            dist_to_bin = fcts.calculate_distance(p.get_pos_xy(), b.get_pos_xy())
                    
            if dist_to_bin < fov and not b.is_full():
                return (True, b.get_id())

        return (False, -1)

    def update_people(self, TIME_STEP: int):
        got_to_destination: list[int] = []
        not_cool_ppl: list[int] = []
        for p in self._ppl:
            time_to_finish_consumption = p.get_time_to_finish_consumption()
            # if he won't need a bin in this step
            if time_to_finish_consumption > TIME_STEP:
                # if he gets to destination
                if p.walk(TIME_STEP):
                    got_to_destination.append(p.get_id())
            else:
                # if he gets to destination before needing a bin
                if p.walk(time_to_finish_consumption):
                    got_to_destination.append(p.get_id())
                # else, he needs to find walk searching for a bin
                else:
                    time_limit = p.get_max_time_of_carrying_trash()
                    t = time_to_finish_consumption
                    (found_a_bin, bin_id) = self.check_for_nearby_bins(p.get_id())
                    while time_limit > 0 and t < TIME_STEP and not found_a_bin:
                        if p.walk(1):
                            got_to_destination.append(p.get_id())
                        (found_a_bin, bin_id) = self.check_for_nearby_bins(p.get_id())
                        time_limit -= 1
                        t += 1
                    if found_a_bin:
                        self._bins[self._id2index(self._bins_ids, bin_id)].put_trash(p.get_trash_volume())
                    elif time_limit == 0:
                        (found_a_bin, bin_id) = self.check_for_nearby_bins(p.get_id())
                        if found_a_bin:
                            self._bins[self._id2index(self._bins_ids, bin_id)].put_trash(p.get_trash_volume())
                        else:
                            self._trash_in_the_strets += p.get_trash_volume()
                            not_cool_ppl.append(p.get_id())
        for id in got_to_destination:
            self.remove_a_person(id)
        for id in not_cool_ppl:
            self.remove_a_person(id)
        


