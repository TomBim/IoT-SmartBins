import bin_lib.map as map
import bin_lib.some_functions as fcs
from bin_lib.consts import *
import numpy as np
import queue

class Map_of_Potentials:
    def __init__(self, mapa: map.Map, n_steps_in_street: int, parent_map: Map_of_Potentials = None):
        self._mapa = mapa
        self._dim0 = len(mapa.get_streets_list())
        self._dim1 = n_steps_in_street + 2
        self._pot_matrix = np.matrix(np.zeros((self._dim0, self._dim1)))
        self._parent_map = parent_map

    def _update_parents(self, add: bool):
        """Add/subtract this map on/from the parent map (map with the totals).
        It should be called before (subtract) and after (add) this maps update

        Args:
            add (bool): add? if false, then it subtracts
        """
        if not self._parent_map == None:
            self._parent_map.son_asking_for_update(add, self._pot_matrix)
        
    def son_asking_for_update(self, add: bool, sons_matrix: np.matrix):
        """Its son is asking to add/subtrat him to/from this map.

        Args:
            add (bool): add? if false, then it subtracts
        """
        if add:
            self._pot_matrix += sons_matrix
        else:
            self._pot_matrix -= sons_matrix

    def update_charge(self, new_charge: float, old_charge: float):
        """updates the maps (this and its parents) changing the value of the charge.
        If you want to update the position of charge, use update_distances

        Args:
            new_charge (float): new value of the charge
            old_charge (float): old value of the charge
        """
        self._update_parents(False)
        self._pot_matrix = self._pot_matrix * new_charge/old_charge
        self._update_parents(True)

    def update_position(self, kq: float, order: int, new_pos_of_charge: map.Pos_Street):
        """updates this map of potential when the charge moved to another position

        Args:
            kq (float): constant k * charge's value
            order (int): order of potential
            new_pos_of_charge (map.Pos_Street): new position of the charge
        """
        self._update_parents(False)
        intersections_info = fcs.dijkstra(self._mapa, new_pos_of_charge)

        # with intersections_info, we can see the distances of the points on the street
        # remembering: each street will be divided in (n_steps_in_street+1) equal parts
        streets = self._mapa.get_streets_list()
        for i in range(self._dim0):
            s = streets[i]
            (A,B) = s.get_intersections_ids()

            # for each point in this division, update the potential using the
            # closest intersection
            division_size = s.get_length() / (self._dim1 - 1)
            dist_using_A = intersections_info[A]['distance_to_origin']
            dist_using_B = intersections_info[B]['distance_to_origin'] + s.get_length()
            for k in range(self._dim1):
                if dist_using_A < dist_using_B:
                    self._pot_matrix[i][k] = kq / dist_using_A**order
                else:
                    self._pot_matrix[i][k] = kq / dist_using_B**order
                dist_using_A += division_size
                dist_using_B -= division_size

        self._update_parents(True)                

class All_Potentials_and_Charges:
    def __init__(self, mapa: map.Map, n_divs_streets: int) -> None:
        self._mapa = mapa
        self._n_divs_streets = n_divs_streets
        self._list_of_potentials: list[Potential] = []
        self._list_of_potentials_orders: list[int] = []
        self._map_pot = Map_of_Potentials(mapa, n_divs_streets)

        self._list_of_charges: list[Pontual_Charge] = []
        self._list_of_charges_ids: list[int] = []
        self._n_charges: int = 0
        self._queue_of_ids_reopen: queue.Queue[int] = queue.Queue()

    def create_potential(self, order: int):
        #TODO: #1 update maps
        n_pots = len(self._list_of_potentials_orders)

        # if there isn't any potential yet, we can just add it
        if n_pots == 0:
            self._list_of_potentials.append(Potential(self._mapa,  self._map_pot, self._n_divs_streets, order))
            self._list_of_potentials_orders.append(order)
            return
        
        # if there is, we have to check if a potential with the same order doesn't already exists
        index = fcs.search_in_vec(self._list_of_potentials_orders, order)
        if index >= n_pots or self._list_of_potentials_orders[index] != order:
            self._list_of_potentials.insert(index, Potential(self._mapa, self._map_pot, self._n_divs_streets, order))
            self._list_of_potentials_orders.insert(index, order)

    def _get_index_of_potential(self, order: int) -> int:
        """gets the index of the order 'order' potential in the
        list 'list_of_potentials_orders'. If it doesn't find, it
        returns the opposite of the value returned by 
        some_functions.search_in_vec (index to the next greater value,
        so you can insert on this index) MINUS 1. The 'minus 1' is
        in case it's 0 but didn't find it.

        Args:
            order (int): order of the potential to be found

        Returns:
            int:
                if didn't find: the opposite of the value returned by 
                    some_functions.search_in_vec (index to the next greater
                    value, so you can insert on this index) MINUS 1 (-search - 1)
                else: the index on the list 'list_of_potentials_orders'
                
        """
        index = fcs.search_in_vec(self._list_of_potentials_orders, order)
        if index >= len(self._list_of_potentials_orders):
            return -index - 1
        
        aux = self._list_of_potentials_orders[index]
        if aux != order:
            return -index - 1
        return index
        
    def create_charge(self, q: float, order_of_potential: int, pos: map.Pos_Street) -> int:
        """create a pontual charge of value 'q' in the positin 'pos'
        that influentiates on the potential with order 'order_of_potential'

        Args:
            q (float): value of the charge
            order_of_potential (int): order of the potential influentiated by the charge
            pos (map.Pos_Street): position of the charge

        Returns:
            int: id of the charge
        """
        #TODO: #1 update maps
        id = self._next_unused_charges_id()
        index_to_insert = fcs.search_in_vec(self._list_of_charges_ids, id)
        index_of_pot = self._get_index_of_potential(order_of_potential)
        if index_of_pot < 0:
            self.create_potential(order_of_potential)
            index_of_pot = -(index_of_pot + 1)

        self._list_of_charges.insert(index_to_insert, Pontual_Charge(q, pos, self._mapa, id,
                                                                     self._list_of_potentials[index_of_pot].get_map_of_pot(),
                                                                     order_of_potential, self._n_divs_streets))
        self._list_of_charges_ids.insert(index_to_insert, id)

        self._n_charges += 1

        return id

    def _next_unused_charges_id(self) -> int:
        if self._n_charges == 0:
            return 0
        
        if not self._queue_of_ids_reopen.empty():
            return self._queue_of_ids_reopen.get()
        
        return self._list_of_charges_ids[-1] + 1

    def _get_index_of_charge(self, id):
        """gets the index of the charge with id 'id' in the
        list 'list_of_potentials_orders'. If it doesn't find, it
        returns the opposite of the value returned by 
        some_functions.search_in_vec (index to the next greater value,
        so you can insert on this index) MINUS 1. The 'minus 1' is
        in case it's 0 but didn't find it.

        Args:
            id (int): id of the charge to be found

        Returns:
            int:
                if didn't find: the opposite of the value returned by 
                    some_functions.search_in_vec (index to the next greater
                    value, so you can insert on this index) MINUS 1 (-search - 1)
                else: the index on the list 'list_of_potentials_orders'
                
        """
        index = fcs.search_in_vec(self._list_of_charges_ids, id)
        if index >= len(self._list_of_charges_ids):
            return -index - 1
        
        aux = self._list_of_charges_ids[index]
        if aux != id:
            return -index - 1
        return index

    def remove_charge(self, id):
        #TODO: #1 update maps
        index = self._get_index_of_charge(id)
        if index >= 0:
            self._list_of_charges.pop(index)
            self._list_of_charges_ids.pop(index)
            self._queue_of_ids_reopen.put(int)

    def move_charge(self, id: int, new_pos: map.Pos_Street):
        pass

class Potential:
    def __init__(self, mapa: map.Map, map_of_pot_total: Map_of_Potentials, n_divs_streets: int, order: int):
        """We don't need to create a constant since it can be inside the value of the carges (like k*q = Q)
        Args:
            mapa (map.Map): map of the city
            n_divs_streets (int): exemples:
                0: updates only intersections
                3: updates intersections and 0.25, 0.5, 0.75 in the pos_in_street
            order (int): order n of potential from V = k.q/d^n
        """
        self._pontual_charges: list[Pontual_Charge] = []
        self._modified_points = []
        self._n_divs_streets = n_divs_streets
        self._dist_in_street_step = 1 / (n_divs_streets + 1)
        self._map_ = mapa
        self._n_streets = len(mapa.get_streets_list())
        self._potential_distributed_on_streets = np.zeros([self._n_streets, (n_divs_streets + 2)])
        self._order = order
        self._map_of_pot = Map_of_Potentials(mapa, n_divs_streets, map_of_pot_total)

    def new_charge(self, q, pos_street, mapa, id):
        self._pontual_charges.append(Pontual_Charge(q, pos_street, mapa, id, self._potential_distributed_on_streets, self._K, self._order, self._n_divs_streets))


    def update_distribution_of_potentials(self):
        for q in self._pontual_charges:
            q.update(self._n_divs_streets)

    def get_map_of_pot(self) -> Map_of_Potentials:
        return self._map_of_pot
        
class Pontual_Charge:
    def __init__(self, q: float, pos_street: map.Pos_Street, mapa: map.Map, id: int, map_pot_parent: Map_of_Potentials, order: int, n_steps_in_street: int) -> None:
        """_summary_

        Args:
            q (float): charge's value
            pos_street (map.Pos_Street): charge's position
            mapa (map.Map): _description_
            id (int): charge's ID
            total_pot_matrix (np.array): map of the sum of every charge's potential (just this kind of potential)
            K (float): constant k from kq/d^n
            order (int): constant n from kq/d^n
            n_divisions_each_street (int): exemples:
                0: updates only intersections
                3: updates intersections and 0.25, 0.5, 0.75 in the pos_in_street

        Other private variables:
            intersections_distance (np.array): vector of the distances from the charge to each intersection
            need_update_distances (bool): need to update the vector of distances to intersections?
            pot_matrix (np.array): map of this charge's potential
            need_update_potentials (bool): need to update the charge's potential map?
        """
        self._q = q
        self._pos_street = pos_street
        self._map_ = mapa
        self._id = id
        self._K = K
        self._order = order
        self._n_divs_streets = n_steps_in_street

        self._total_pot_matrix = total_pot_matrix
        self._pot_matrix = np.array([])
        self._need_update_potentials = True

        self._intersections_distance = []
        self._need_update_distances = True
        
        self.update(n_steps_in_street)

    def set_q(self, new_q: float):
        self._q = new_q
        self._need_update_potentials = True

    def set_position(self, new_pos_street: map.Pos_Street):
        self._pos_street = new_pos_street

    def set_need_update(self, distances: bool, potentials: bool):
        """You can't set the neediness as false in this function. If you want to do it, use set_updated

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
        """You can't set the neediness as true in this function. If you want to do it, use set_need_update

        Args:
            distances (bool): set as distances updated?
            potentials (bool): set as potentials updated?
        """
        if distances:
            self._need_update_distances = False
        if potentials:
            self._need_update_potentials = (False or self._need_update_distances)

    def update(self, n_steps_in_street: int):
        """update the distances vector and the potentials map if needed

        Args:
            n_divisions_each_street (int): exemples:
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
                    intersections_ids = streets[i].get_intersections_ids()
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
                    intersections_ids = streets[i].get_intersections_ids()
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
            d = self._intersections_distance[self._map_.get_street(0).get_intersections_ids()[0]]
            if d > EPSILON:
                old_Kq = self._pot_matrix[0,0] * d**self._order
            else:
                d = self._intersections_distance[self._map_.get_street(0).get_intersections_ids()[1]]
                old_Kq = self._pot_matrix[0,(n_steps_in_street + 1)] * d**self._order

            # update local matrix and, afther, the global one
            self._pot_matrix = self._pot_matrix  * (self._K * self._q / old_Kq)
            self._need_update_potentials = False
            self._total_pot_matrix += self._pot_matrix