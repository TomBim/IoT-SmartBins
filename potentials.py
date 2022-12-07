import bin_lib.map as map
import bin_lib.some_functions as fcs
from bin_lib.consts import *
import numpy as np

class Potential:
    def __init__(self, mapa: map.Map, n_steps_in_stret: int, constant: float):
        """
        Args:
            mapa (map.Map): _description_
            n_steps_in_stret (int): exemples:
                0: updates only intersections
                3: updates intersections and 0.25, 0.5, 0.75 in the pos_in_street
        """
        self._pontual_charges: list[Pontual_Charge] = []
        self._modified_points = []
        self._n_steps_in_street = n_steps_in_stret
        self._dist_in_street_step = 1 / (n_steps_in_stret + 1)
        self._map_ = mapa
        self._potential_distributed_on_streets = np.zeros([])
        self._n_streets = len(mapa.get_streets_list())
        self._potential_distributed_on_streets = np.zeros([self._n_streets, (n_steps_in_stret + 2)])
        self._K = constant

    def new_charge(self, q, pos_street, mapa, id):
        self._pontual_charges.append(Pontual_Charge(q, pos_street, mapa, id))


    def update_distribution_of_potentials(self):
        for q in self._pontual_charges:
            q.update
        

class Pontual_Charge:
    def __init__(self, q: float, pos_street: map.Pos_Street, mapa: map.Map, id: int) -> None:
        self._q = q
        self._pos_street = pos_street
        self._need_update_distances = True
        self._need_update_potentials = True
        self._map_ = mapa
        self._id = id
        self.update()

    def set_q(self, new_q: float):
        self._q = new_q
        self._need_update_potentials = True

    def set_need_update(self, distances: bool, potentials: bool):
        """You can't set the neediness as false in this function. If you want to, use set_updated

        Args:
            distances (bool): set as needing an update on distances?
            potentials (bool): set as needing an update on potentials?
        """
        if distances:
            self._need_update_distances = True
            self._need_update_potentials = True
        elif potentials:
            self._need_update_potentials = True

    def set_updated(self, distances: bool, potentials: bool):
        """You can't set the neediness as true in this function. If you want to, use set_need_update

        Args:
            distances (bool): set as distances updated?
            potentials (bool): set as potentials updated?
        """
        if distances:
            self._need_update_distances = False
        if potentials:
            self._need_update_potentials = (False or self._need_update_distances)

    def update(self):
        if self._need_update_distances:
            self._intersections_distance = fcs.dijkstra(self._map_, self._pos_street)['distance to origin']
            self._need_update_distances = False
            for intersection_dist in self._intersections_distance:
                self._intersections_potentials = Ke * self._q / intersection_dist
            self._need_update_potentials = False
        elif self._need_update_potentials:
            for intersection_dist in self._intersections_distance:
                self._intersections_potentials = Ke * self._q / intersection_dist
            self._need_update_potentials = False

    def update_pot_distr(self, matrix: np.array):
        """update the matrix if it's necessary

        Args:
            matrix (np.array): matrix with the potential distribution
        """
        if self._need_update_distances or self._need_update_potentials

    def calculate_potential(self, pos_street)

class Gravity_Point

class Gravity_Potential

class Electrical_Point:

class Electrical_Potential:
