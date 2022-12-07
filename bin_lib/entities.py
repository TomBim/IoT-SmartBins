from __future__ import annotations
from math import inf, sqrt
import heapq
import random as rand
import math as m
from bin_lib.consts import *
import bin_lib.map as map
import bin_lib.some_functions as fcs
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
        pq: list[tuple[float, map.Intersection]] = []
        starting_intersections = self._origin_street.get_street().get_vector()

        for intersection in starting_intersections:
            distance = fcs.calculate_distance(self._origin.get_pos_xy(), intersection.get_pos())
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
        dist_destination_intersection1 = fcs.calculate_distance(destination_intersections[0].get_pos(), self._destination.get_pos_xy())
        dist_destination_intersection2 = fcs.calculate_distance(destination_intersections[1].get_pos(), self._destination.get_pos_xy())
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
        reversed_path.append(self._map_.get_intersection(intersection_id))
        
        return reversed_path[::-1]

    def _calculate_distance_to_destination(self):
        start = self.get_pos_xy()
        n_intersections = len(self._path)
        end = self._destination_street.get_pos_xy()
        d = 0
        for i in range(n_intersections):
            aux = self._path[i].get_pos()
            d += fcs.calculate_distance(start, aux)
            start = aux
        return d + fcs.calculate_distance(start, end)

    def get_path(self):
        return self._path

    def get_origin(self):
        return self._origin

    def get_destination(self):
        return self._destination

    def has_trash(self) -> bool:
        return self._has_trash
    
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
            return True

        n_intersections = len(self._path)
        self._time_alive += time_walking
        self._distance_to_destination -= d

        # he is already on the destination's street?
        if n_intersections == 0:
            s = self._pos_street.get_street()
            (A, B) = s.get_vector()
            (Ax, Ay) = A.get_pos()
            (Bx, By) = B.get_pos()
            if abs(Ax-Bx) > EPSILON: #se nao é reta vertical
                m = (self._destination_street.get_pos_xy()[0] - start[0])*(Bx-Ax)
                if m > 0: # apontam no mesmo sentido
                    new_pos_in_street = self._pos_street.get_pos_in_street() + d / s.get_length()
                else: # se apontam em sentidos opostos
                    new_pos_in_street = self._pos_street.get_pos_in_street() - d / s.get_length()
            else: #se for vertical
                new_pos_in_street = self._pos_street.get_pos_in_street() + d / (By - Ay)
            self._pos_street = map.Pos_Street(s, new_pos_in_street)
            return False
        
        # walk through the intersections
        last_intersection = None
        next_intersection = self._path[0] # n_intersections != 0 cause of return right above
        distance_to_next_intersection = fcs.calculate_distance(start, next_intersection.get_pos())
        while d > distance_to_next_intersection and n_intersections > 1:
            d -= distance_to_next_intersection
            last_intersection = next_intersection
            self._path.pop(0)
            n_intersections -= 1
            next_intersection = self._path[0]
            distance_to_next_intersection = fcs.calculate_distance(last_intersection.get_pos(), next_intersection.get_pos())
        if d > distance_to_next_intersection:
            d -= distance_to_next_intersection
            last_intersection = next_intersection
            self._path.pop(0)
            n_intersections = 0
        
        # walk on last street that he will use in this time step
        # more intersections to pass thourgh?
        if n_intersections > 0:
            next_intersection = self._path[0]
            # if he didn't change of street in this time step
            if last_intersection == None:
                s = self.get_pos_street().get_street()
                old_pos_in_street = self.get_pos_street().get_pos_in_street()
            # if he did change
            else:
                s = self._map_.search_for_a_street(last_intersection, next_intersection)
                start = last_intersection.get_pos()
                old_pos_in_street = fcs.calculate_distance(last_intersection.get_pos(), s.get_vector()[0].get_pos()) / s.get_length()
            destin = next_intersection.get_pos()
        # then he's at destination's street
        else:
            s = self._destination_street.get_street()
            destin = self._destination_street.get_pos_xy()
            # if he changed of street in this time step
            if not (last_intersection == None):
                start = last_intersection.get_pos()
                old_pos_in_street = fcs.calculate_distance(last_intersection.get_pos(), s.get_vector()[0].get_pos()) / s.get_length()
        (A, B) = s.get_vector()
        (Ax, Ay) = A.get_pos()
        (Bx, By) = B.get_pos()
        if abs(Ax-Bx) > EPSILON: # se nao é reta vertical
            m = (destin[0] - start[0])*(Bx-Ax)
            if m > 0: #apontam para o mesmo lugar
                new_pos_in_street = old_pos_in_street + d / s.get_length()
            else:#se apontam em sentidos opostos
                new_pos_in_street = old_pos_in_street - d / s.get_length()
        else: #se for vertical
            m = (destin[1] - start[1])*(By-Ay)
            if m > 0:
                new_pos_in_street = old_pos_in_street + d / (By - Ay)
            else:
                new_pos_in_street = old_pos_in_street - d / (By - Ay)
                
        # just to be cautious with float problems
        if new_pos_in_street >= 1:
            new_pos_in_street = 1 - EPSILON
        elif new_pos_in_street <= 0:
            new_pos_in_street = EPSILON
        self._pos_street = map.Pos_Street(s, new_pos_in_street)
        return False

    def get_fov(self):
        return self._fov

    def get_trash_volume(self):
        return self._trash_volume

    def set_true_trash(self):
        self._has_trash = True

    def increase_time_with_trash(self, increase: int):
        self._time_carrying_trash += increase
    
    def get_time_with_trash(self) -> int:
        return self._time_carrying_trash

    def get_speed(self) -> float:
        return self._speed

class Bin(Entity):
    """Private variables:
        id: int
        pos: float[2]
        pos_street: [street,float]
        capacity: float (L)
        filling_rate: float (L/s)
        full: bool
        vol_trash: float (L)
        last_time_was_emptied: int (seconds)
    """
    def __init__(self, id: int, pos_street: map.Pos_Street, capacity: float):
        super().__init__(id, pos_street)
        self._capacity = capacity
        self._full = False
        self._vol_trash = 0
        self._last_time_was_emptied = -1
        self._filling_rate = 0

    def is_full(self):
        return self._full

    def set_capacity(self, new_capacity: float):
        self._capacity = new_capacity

    def put_trash(self, vol_trash):
        self._vol_trash += vol_trash
        if(self._vol_trash > self._capacity):
            self._full = True
        #update filling_rate

    def empty_bin(self, time_now: int) -> float:
        """it empties a bin and returns the volume of trash the bin had b4 being emptied

        Returns:
            float: vol of trash b4 getting emptied
        """
        x = self._vol_trash
        self._vol_trash = 0
        self._full = False
        self._last_time_was_emptied = time_now
        return x

    def get_vol_trash(self):
        return self._vol_trash

    def get_percentage(self):
        return self._vol_trash / self._capacity

    def get_last_time_was_emptied(self):
        return self._last_time_was_emptied

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

class Street_Sweepers:
    def __init__(self, map: map.Map):
        self._trash_got_from_floor = 0
        self._total_cost = 0
        self._streets_cleaned = 0
        self._map_ = map
    
    def clean_streets(self, trash_on_floor_vol: float, trash_on_floor_pos: list[map.Pos_Street]):
        # streets_to_be_cleaned: list[tuple[int, int]] = []
        # for trash in trash_on_floor_pos:
        #     (intersection_A_id, intersection_B_id) = trash.get_street().get_intersctions_ids()
        #     good = True
        #     stop = False
        #     i = 0
        #     while i < len(streets_to_be_cleaned) and good and not stop:
        #         s = streets_to_be_cleaned[i]
        #         if intersection_A_id == s[0]:
        #             if intersection_B_id == s[1]:
        #                 good = False
        #                 stop = True
        #             elif intersection_B_id > s[1]:
        #                 stop = True
        #         elif intersection_A_id > s[0]:
        #             stop = True
        #         i += 1
        #     if stop and good:
        #         streets_to_be_cleaned.insert((intersection_A_id, intersection_B_id))

        streets_to_be_cleaned = {trash.get_street() for trash in trash_on_floor_pos}

        n = len(streets_to_be_cleaned)
        self._streets_cleaned += n
        self._total_cost += n * COST_OF_CLEANING_A_STREET
        self._trash_got_from_floor += trash_on_floor_vol
        
class Trash_Trucks:
    def __init__(self, map: map.Map, everything: Everything):
        self._trash_got_from_bins = 0
        self._total_cost = 0
        self._km_per_day = 0
        self._bins_emptied = 0
        self._map_ = map
        self._everything_ = everything
        self._days_used = 0

    def empty_bins(self, everything: Everything, time_now: int):
        bins_emptied: list[int] = []
        bins_list = everything.get_bins_list()
        for bin in bins_list:
            if bin.get_percentage >= MIN_PERCENTAGE_BIN:
                bins_emptied.append(bin.get_id())
                bin.empty_bin(time_now)
            elif time_now - bin.get_last_time_was_emptied() >= MAX_TIME_IGNORING_A_BIN:
                bins_emptied.append(bin.get_id())
                bin.empty_bin(time_now)
        if len(bins_emptied) > 0:
            total_dist = self._total_distance_traveled(bins_emptied)/1e3
            self._total_cost += total_dist * COST_TRUCK_PER_KM + COST_TRUCK_PER_DAY
            self._km_per_day = (self._km_per_day * self._days_used + total_dist) / (self._days_used + 1)
            self._days_used += 1

    def _dijkstra_intersections(self, start: map.Pos_Street):
        intersections_info = [{'distance_to_origin': inf, 'parent': None, 'closed': False} for _ in self._map_.get_intersections_list()]
        pq: list[tuple[float, map.Intersection]] = []
        starting_intersections = start.get_street().get_vector()

        for intersection in starting_intersections:
            distance = fcs.calculate_distance(start.get_pos_xy(), intersection.get_pos())
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

    def _next_bin(self, bins_to_empty: list[int], start_bin: int) -> tuple[int, float]:
        intersections_info = self._dijkstra_intersections(self._everything_.get_bin(start_bin))        
        next_bin_distance = inf
        next_bin = -1
        for bin in bins_to_empty:
            intersections_ids = self._everything_.get_bin(bin).get_pos_street().get_street().get_intersctions_ids()
            if intersections_info[intersections_ids[0]]['distance to origin'] > intersections_info[intersections_ids[1]]['distance to origin']:
                dist = intersections_info[intersections_ids[1]]['distance to origin']
            else:
                dist = intersections_info[intersections_ids[0]]['distance to origin']
            if dist < next_bin_distance:
                next_bin = bin
                next_bin_distance = dist
        return (next_bin, next_bin_distance)

    def _total_distance_traveled(self, bins_to_empty: list[int]):
        start = bins_to_empty[0]
        total_dist = 0
        while len(bins_to_empty) > 0:
            bins_to_empty.pop(0)
            (next_bin, next_bin_distance) = self._next_bin(self, bins_to_empty, start)
            total_dist += next_bin_distance
            for i in range(len(bins_to_empty)):
                if bins_to_empty[i] == next_bin:
                    bins_to_empty.pop(i)
                    break
        return total_dist   
        
class Everything:
    def __init__(self, mapa: map.Map) -> None:
        self._ppl: list[Person] = []
        self._com_points: list[Commercial_Point] = []
        self._bins: list[Bin] = []
        self._ppl_ids: list[int] = []
        self._bins_ids: list[int] = []
        self._food_points = 0
        self._nfood_points = 0
        self._job_points = 0
        self._n_people = 0
        self._n_com_points = 0
        self._n_bins = 0
        self._last_persons_id = -1
        self._last_bins_id = -1
        self._com_points_attractiveness = np.array([])
        self._trash_in_the_streets: float = 0
        self._map = mapa
        self._pos_trash_floor: list[map.Pos_Street] = []
        self._street_sweepers = Street_Sweepers(mapa)
        self._trash_trucks = Trash_Trucks(mapa, self)
    
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
                payload["max_time_carrying_trash"] = int(rand.gauss(mu=38, sigma=18))
            payload["speed"] = abs(rand.gauss(mu = 1.25, sigma = .25))
            payload["trash_volume"] = trash_info[0]
            payload["time_of_consumption"] = trash_info[1]

            id = self._next_unused_person_id()
            p = Person(id, payload["origin"], payload["destination"], self._map, payload)
            index = self._id2index(self._ppl_ids, id)
            self._ppl.insert(index, p)
            self._ppl_ids.insert(index, id)
            self._last_persons_id += 1
            self._n_people += 1

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
        return fcs.search_in_vec(id_vec, id)

    def remove_a_person(self, id: int):
        if len(self._ppl) != 0:
            index = self._id2index(self._ppl_ids, id)
            self._ppl.pop(index)
            self._ppl_ids.pop(index)
            self._n_people -= 1

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
        p = self._ppl[fcs.search_in_vec(self._ppl_ids, id_person)]
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

    def check_for_nearby_bins(self, p: Person) -> tuple[bool, int]: 
        """
        This function checks if there is any bin closer than their fov
        If there is a bin nearby, it should return True 
        Else, it should return False
        """
        fov = p.get_fov()

        for b in self._bins:
            dist_to_bin = fcs.calculate_distance(p.get_pos_xy(), b.get_pos_xy())
                    
            if dist_to_bin < fov and not b.is_full():
                return (True, b.get_id())

        return (False, -1)

    def _update_a_person(self, p: Person, TIME_STEP: int) -> bool:
        """

        Args:
            person (Person): _description_
            TIME_STEP (int): _description_

        Returns:
            bool: pop?
        """
        step = int(1.5 * p.get_fov() / p.get_speed())
        if step > TIME_STEP:
            step = TIME_STEP
        # if he didn't finish his comsumption
        if not p.has_trash():
            time_to_finish_consumption = p.get_time_to_finish_consumption()
            # if he won't need a bin in this step
            if time_to_finish_consumption > TIME_STEP:
                # if he gets to destination
                return p.walk(TIME_STEP)
        
            # make him walk while finish his consuption
            
            # if he gets to destination before needing a bin
            if p.walk(time_to_finish_consumption):
                return True
            # else, he needs to start walk while search for a bin
            p.set_true_trash()
            time_limit = p.get_max_time_of_carrying_trash()
            t = time_to_finish_consumption
            (found_a_bin, bin_id) = self.check_for_nearby_bins(p)
            while (time_limit - step >= 0) and (t + step <= TIME_STEP) and not found_a_bin:
                if p.walk(step):
                    return True
                (found_a_bin, bin_id) = self.check_for_nearby_bins(p)
                time_limit -= step
                t += step
                p.increase_time_with_trash(step)
            while (time_limit > 0) and (t < TIME_STEP) and not found_a_bin:
                if p.walk(1):
                    return True
                (found_a_bin, bin_id) = self.check_for_nearby_bins(p)
                time_limit -= 1
                t += 1
                p.increase_time_with_trash(1)
            if found_a_bin:
                self._bins[self._id2index(self._bins_ids, bin_id)].put_trash(p.get_trash_volume())
                return True
            if time_limit == 0:
                (found_a_bin, bin_id) = self.check_for_nearby_bins(p)
                if found_a_bin:
                    self._bins[self._id2index(self._bins_ids, bin_id)].put_trash(p.get_trash_volume())
                else:
                    self._trash_in_the_streets += p.get_trash_volume()
                    self._pos_trash_floor.append(p.get_pos_street())
                return True
            # else, t == TIME_STEP
            return False
        
        # finished his consumption already, before starting this
        time_limit = p.get_max_time_of_carrying_trash() - p.get_time_with_trash()
        t = 0
        (found_a_bin, bin_id) = self.check_for_nearby_bins(p)
        while (time_limit - step  >= 0) and (t + step < TIME_STEP) and not found_a_bin:
            if p.walk(step):
                return True
            (found_a_bin, bin_id) = self.check_for_nearby_bins(p)
            time_limit -= step
            t += step
            p.increase_time_with_trash(step)
        while time_limit > 0 and t < TIME_STEP and not found_a_bin:
            if p.walk(1):
                return True
            (found_a_bin, bin_id) = self.check_for_nearby_bins(p)
            time_limit -= 1
            t += 1
            p.increase_time_with_trash(1)
        if found_a_bin:
            self._bins[self._id2index(self._bins_ids, bin_id)].put_trash(p.get_trash_volume())
            return True
        if time_limit <= 0:
            (found_a_bin, bin_id) = self.check_for_nearby_bins(p)
            if found_a_bin:
                self._bins[self._id2index(self._bins_ids, bin_id)].put_trash(p.get_trash_volume())
            else:
                self._trash_in_the_streets += p.get_trash_volume()
                self._pos_trash_floor.append(p.get_pos_street())
            return True
        
        # he's with trash, but didn't lose all patience and time is gone
        return False
    
    def update_people(self, TIME_STEP: int):
        # n_pops = 0
        # for i in range(len(self._ppl)):
        #     p = self._ppl[i-n_pops]
        #     if self._update_a_person(p, TIME_STEP):
        #         self._ppl.pop(i-n_pops)
        #         self._ppl_ids.pop(i-n_pops)
        #         self._n_people -= 1
        #         n_pops += 1
        
        new_ppl = []
        new_ppl_ids = []
        for p in self._ppl:
            # antes = p.get_pos_xy()
            if not self._update_a_person(p, TIME_STEP):
            # if antes != p.get_pos_xy():
                new_ppl.append(p)
                new_ppl_ids.append(p.get_id())
        self._ppl = new_ppl
        self._ppl_ids = new_ppl_ids


        # new_ppl = []
        # new_ppl_ids = []
        # for p in self._ppl:
        #     antes = p.get_pos_xy()
        #     self._update_a_person(p, TIME_STEP)
        #     if antes != p.get_pos_xy():
        #         new_ppl.append(p)
        #         new_ppl_ids.append(p.get_id())
        # self._ppl = new_ppl
        # self._ppl_ids = new_ppl_ids
        
    def get_trash_in_the_street(self):
        return self._trash_in_the_streets
        
    def get_bin(self, id: int) -> Bin:
        return self._bins[self._id2index(self._bins_ids, id)]

    def get_bins_list(self) -> list[Bin]:
        return self._bins

    def sweep_streets(self):
        pass
    
    def empty_bins(self):
        pass