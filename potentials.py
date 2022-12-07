import bin_lib.map as map
import bin_lib.some_functions as fcs
from bin_lib.consts import *
import numpy as np

class Potential:
    def __init__(self, mapa: map.Map, n_steps_in_stret: int, constant: float, order: int):
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
        self._n_streets = len(mapa.get_streets_list())
        self._potential_distributed_on_streets = np.zeros([self._n_streets, (n_steps_in_stret + 2)])
        self._K = constant
        self._order = order

    def new_charge(self, q, pos_street, mapa, id):
        self._pontual_charges.append(Pontual_Charge(q, pos_street, mapa, id, self._potential_distributed_on_streets, self._K, self._order, self._n_steps_in_street))


    def update_distribution_of_potentials(self):
        for q in self._pontual_charges:
            q.update(self._n_steps_in_street)
        
class Pontual_Charge:
    def __init__(self, q: float, pos_street: map.Pos_Street, mapa: map.Map, id: int, total_pot_matrix: np.array, K: float, order: int, n_steps_in_street: int) -> None:
        self._q = q
        self._pos_street = pos_street
        self._map_ = mapa
        self._id = id
        self._total_pot_matrix = total_pot_matrix
        self._K = K
        self._order = order
        self._need_update_distances = True
        self._need_update_potentials = True
        self._intersections_distance = []
        self._pot_matrix = np.array([])
        self._n_steps_in_street = n_steps_in_street
        self.update(n_steps_in_street)

    def set_q(self, new_q: float):
        self._q = new_q
        self._need_update_potentials = True

    def set_position(self, new_pos_street: map.Pos_Street):
        self._pos_street = new_pos_street

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

    def update(self, n_steps_in_street: int):
        """_summary_

        Args:
            n_steps_in_stret (int): exemples:
                0: updates only intersections
                3: updates intersections and 0.25, 0.5, 0.75 in the pos_in_street
        """
        # needs update on vector of intersections_distance
        if self._need_update_distances:
            # if it's not initializing, we must take out the values from total matrix b4 continuing
            if len(self._intersections_distance) > 0:
                self._total_pot_matrix -= self._pot_matrix

            streets = self._map_.get_streets_list()
            n_streets = len(streets)

             # update intersections_distance
            intersections_info = fcs.dijkstra(self._map_, self._pos_street)
            # print(intersections_info)
            self._intersections_distance = np.array([])
            for i in range(len(intersections_info)):
                self._intersections_distance = np.append(self._intersections_distance,intersections_info[i]['distance_to_origin'])
            # print(np.array(intersections_info[0:-1]['distance to origin']))
            # self._intersections_distance = np.array(intersections_info[:]['distance to origin'])
            self._need_update_distances = False

            # update the matrix with potentials from this charge
            # first, we make a matrix with the distances
            this_streets_index = self._pos_street.get_street().get_id()
            m_pot_self = np.zeros([n_streets, (n_steps_in_street + 2)])
            for i in range(n_streets):
                # if it's this street, the distance gets lowers inside the street
                if i == this_streets_index:
                    length_street = streets[i].get_length()
                    intersections_ids = streets[i].get_intersctions_ids()
                    d0 = self._intersections_distance[intersections_ids[0]]
                    d1 = self._intersections_distance[intersections_ids[1]]
                    dist_step = length_street / (n_steps_in_street + 1)
                    k_max = n_steps_in_street + 1
                    for k in range(k_max + 1):
                        x = d0 - k*dist_step
                        y = d1 - (k_max-k)*dist_step
                        d = x if x>0 else y
                        if d < EPSILON:
                            d = EPSILON
                        m_pot_self[i,k] = d**self._order
                # general case: the distance inside a street is the smaller distance
                # (from the point in street to an intersection) + (from the intersection to the charge)
                else:
                    length_street = streets[i].get_length()
                    intersections_ids = streets[i].get_intersctions_ids()
                    d0 = self._intersections_distance[intersections_ids[0]]
                    d1 = self._intersections_distance[intersections_ids[1]]
                    dist_step = length_street / (n_steps_in_street + 1)
                    k_max = n_steps_in_street + 1
                    for k in range(k_max + 1):
                        x = d0 + k*dist_step
                        y = d1 + (k_max-k)*dist_step
                        d = x if x<y else y
                        if d < EPSILON:
                            d = EPSILON
                        m_pot_self[i,k] = d**self._order
            self._pot_matrix = self._K * self._q / m_pot_self
            self._need_update_potentials = False

            # add to the total matrix
            self._total_pot_matrix += self._pot_matrix
        # if we need to update only the local matrix of potentials
        elif self._need_update_potentials:
            # remove old local matrix from global matrix
            self._total_pot_matrix -= self._pot_matrix

            # calculate old K*q
            d = self._intersections_distance[self._map_.get_street(0).get_intersctions_ids()[0]]
            if d > EPSILON:
                old_Kq = self._pot_matrix[0,0] * d**self._order
            else:
                d = self._intersections_distance[self._map_.get_street(0).get_intersctions_ids()[1]]
                old_Kq = self._pot_matrix[0,(n_steps_in_street + 1)] * d**self._order

            # update local matrix and, afther, the global one
            self._pot_matrix = self._pot_matrix  * (self._K * self._q / old_Kq)
            self._need_update_potentials = False
            self._total_pot_matrix += self._pot_matrix