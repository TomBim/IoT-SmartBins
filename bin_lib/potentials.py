from __future__ import annotations
import numpy as np
import queue
import math as m
import copy

import bin_lib.map as map
import bin_lib.some_functions as fcs
import bin_lib.entities as entities
from bin_lib.consts import *


class Map_of_Potentials:
    def __init__(self, mapa: map.Map, order: int, map_empty: bool = True,
                 pot_map: list[np.array] = [], empty_street: list[bool] = []):
        """Makes deepcopy of pot_map and empty_street if those are given.

        Args:
            mapa (map.Map): map of the city
            order (int): order of potential
            map_empty (bool, optional): is this map of pot empty?. Defaults to True.
            pot_map (list[np.array], optional): distribution of potential on the map of the city. Defaults to [].
            empty_street (list[bool], optional): is the street 'i' empty? Defaults to [].
        """
        self._mapa = mapa
        self._order = order
        self._div_size = SIZE_OF_DIVISION[order-1]
        self._n_streets = len(mapa.get_streets_list())

        if map_empty:
            self._map_empty = True
            for i in range(self._n_streets):
                self._pot_map.append(np.array([]))
                self._empty_street.append(True)
        else:
            self._pot_map: list[np.array] = copy.deepcopy(pot_map)
            self._empty_street: list[bool] = copy.deepcopy(empty_street)
            self._map_empty = False

    def __add__(self, y) -> Map_of_Potentials:
        """Used for helping with Map1+Map2 and returns Map_res whose _pot_map is the sum of
        both's _pot_map's.
        In the case of different orders, it uses the functions "interpolate(new_order)" to let both maps
        with the same shape.
        CAREFUL: if you do s = x+y, you will lose the reference to old s.
        CAREFUL: this function considers EPSILON_POT. Therefore, it goes to zero if a sum < EPSILON_POT
        OBS: (x+y) is a new map_of_potential, don't need any copy/deepcopy

        Args:
            y (_type_): should be an instance of Map_of_Potentials (can be subclass)

        Raises:
            TypeError: if order isn't the same
            TypeError: if both arguments of the sum aren't Map_of_Potentials

        Returns:
            Map_of_Potentials: map with the sum of pot_map's
        """
        if isinstance(y, Map_of_Potentials):
            # if any of them is empty, it's trivial
            if self.is_map_empty() and y.is_map_empty():
                return Map_of_Potentials(self._mapa, self._order)
            
            if self.is_map_empty():
                return Map_of_Potentials(y._mapa, y._order, False, y._pot_map, y._empty_street)
            
            if y.is_map_empty():
                return Map_of_Potentials(self._mapa, self._order, False, self._pot_map, self._empty_street)
            
            # if none of them is empty, let's calculate it!!!
            # if they are from the same order of potential
            if y.get_order() == self._order:
                res = Map_of_Potentials(self._mapa, self._order)
                for s in range(self._n_streets):
                    x_empty = self._empty_street[s]
                    y_empty = y._empty_street[s]
                    if x_empty:
                        if not y_empty:
                            res._pot_map[s] = y._pot_map[s]
                            res._empty_street[s] = False
                    else:
                        if y_empty:
                            res._pot_map[s] = self._pot_map[s]
                            res._empty_street[s] = False
                        else:
                            # After making the sum, we need to see,
                            # if the result doesn't empty any street.
                            # Since empting street is good for performance,
                            # we don't need to change the not empty streets
                            aux = self._pot_map[s] + y._pot_map[s]
                            if np.allclose(aux, 0, atol=EPSILON_POT):
                                res._pot_map[s] = 0*aux
                                res._empty_street[s] = True
                            else:
                                res._pot_map[s] = aux
                                res._empty_street[s] = False
                    if np.all(res._empty_street):
                        res._map_empty = True
                return res
            # if its not the same order
            else:
                # we need to interpolate the lower order, so we can make the sum
                if y.get_order() < self._order:
                    a = self
                    b = y.interpolate(self._order)
                else:
                    a = self.interpolate(y.get_order())
                    b = y
                return a+b
        raise TypeError('Addition available only for same type')
    
    def __radd__(self, y) -> Map_of_Potentials:
        return self.__add__(y)
    
    def __mul__(self, y) -> Map_of_Potentials:
        """Helps with a*Map, returning a Map whose _pot_map = a * Map._pot_map.
        CAREFUL: if you do p = x*y, you will lose the reference to the old p.
        OBS: don't need any copy/deepcopy, the function creates a new map_of_potential
        OBS: it doesn't consider EPSILON_POT. If you want to, do it after. But in case of
            updating charge's value, it isn't so cool, because you will lose information
            and won't be able to just multiply by new_q/old_q

        Args:
            y (int/float): must be a scalar

        Raises:
            TypeError: if y isn't a scalar (int/float)

        Returns:
            Map_of_Potentials: p._pot_map = y * Map._pot_map (updates booleans too)
        """
        if isinstance(y, (int, float)):
            if y == 0:
                return Map_of_Potentials(self._mapa, self._order)
            
            res = Map_of_Potentials(self._mapa, self._order)
            for s in range(self._n_streets):
                if not self._empty_street[s]:
                    res._pot_map[s] = y * self._pot_map[s]
                    res._empty_street[s] = False
            return res
        raise TypeError('Multiplication available only by scalars, not lists/vectors/others')
    
    def __rmul__(self, y) -> Map_of_Potentials:
        return self.__mul__(y)

    def get_order(self) -> int:
        return self._order

    def is_map_empty(self) -> bool:
        return self._map_empty
    
    def _clear_map(self):
        empty_map = Map_of_Potentials(self._mapa, self._order)
        self._pot_map = empty_map._pot_map
        self._empty_street = empty_map._empty_street
        self._map_empty = True
    
    def get_n_steps(self, street: map.Street) -> int:
        """returns the number of steps in a specific street. It does it using the SIZE_OF_DIVISION
        from the order of the map. Counts the point 0.

        Args:
            street (map.Street): the street whose n_steps you want

        Returns:
            int: number of steps in the street (counting the zero)
        """
        return m.ceil(street.get_length() / self._div_size) + 1

    def interpolate(self, new_order) -> Map_of_Potentials:
        interpol = Map_of_Potentials(self._mapa, new_order)
        if not self._map_empty:
            interpol._map_empty = False
            streets = self._mapa.get_streets_list()
            for i in range(self._n_streets):
                if not self._empty_street[i]:
                    interpol._empty_street[i] = False
                    s = streets[i]
                    n_steps_higher = interpol.get_n_steps(s)
                    n_steps_lower = self.get_n_steps(s)
                    ratio = n_steps_higher / n_steps_lower
                    step_lower = -1
                    for step_higher in range(n_steps_higher):
                        # let's make "mean = a*x + b"
                        if not step_lower == step_higher // ratio:
                            step_lower = step_higher // ratio
                            a = self._pot_map[i][step_lower+1] - self._pot_map[i][step_lower]
                            b = self._pot_map[i][step_lower]
                        mean = a * (step_higher % ratio) + b
                        if abs(mean) < EPSILON_POT:
                            mean = 0
                        interpol._pot_map[i][step_higher] = mean
        return interpol
                    

class Map_of_Pot_One_Order(Map_of_Potentials):
    def __init__(self, mapa: map.Map, order: int):
        """Each order of potential has its own map, so it's easier to
        update and make sums, subtractions etc.

        Args:
            mapa (map.Map): map of the city
            order (int): order of the potential represented by this map
        """
        super().__init__(mapa, order)

    def charge_asking_update(self, diff: Map_of_Potentials):
        """Updates this map of potential considering only the difference
        to add (diff map = new charges' map - old charges' map).

        Args:
            diff (Map_of_Potentials): _description_
        """
        res = self + diff
        self._pot_map = res._pot_map
        self._empty_street = res._empty_street
        self._map_empty = res._map_empty

class Map_of_Pot_for_Charges(Map_of_Potentials):
    def __init__(self, mapa: map.Map, order: int, parent_map: Map_of_Pot_One_Order, intersections_dists: np.array,
                 pos_of_charge: map.Pos_Street):
        """Each charge has its own maps: one for each order. That way, it's easier to update when 
        we change the charge's value or we move it, etc.

        Args:
            mapa (map.Map): map of the city
            order (int): order which this map represents
            parent_map (Map_of_Pot_One_Order): map of pot of the same order
            pos_of_charge (map.Pos_Street): position of the charge
        """
        super().__init__(mapa, order)
        self._parent_map = parent_map
        self._intersections_dists = intersections_dists
        self._pos = pos_of_charge

    def update_charges_value(self, new_q: float, old_q: float):
        """FOR FIRST UPDATE, USE first_update function.
        updates the maps (this and its parents) changing the value of the charge.
        If you want to update the position of charge, use update_distances

        Args:
            new_charge (float): new value of the charge
            old_charge (float): old value of the charge
        """
        # if map was empty, then need to calculates the distances
        if old_q == 0 or self._map_empty:
            self.first_update(new_q, self._pos)
            return
        
        # if new_q is less than EPSILON_CHARGE, then we put 0, so the map stay empty
        if abs(new_q) < EPSILON_CHARGE:
            new_q = 0
        
        # updating
        ratio = new_q / old_q
        new_map = self * ratio
        diff = self * (ratio - 1)
        self._pot_map = new_map._pot_map
        self._parent_map.charge_asking_update(diff)

        # if new_q == 0, we need to make the map empty
        if ratio == 0:
            self._empty_street = new_map._empty_street
            self._map_empty = True

    def update_charges_position(self, q: float, new_pos_of_charge: map.Pos_Street,
                                first_update: bool = False):
        """FOR FIRST UPDATE, USE first_update function.
        updates this map of potential when the charge moved to another position.
        ATTENTION: the actual division size used here is a little different:
            n_steps = ceil(street_length / division_size_from_constants  +  1)
            actual_div_size = street_length / n_steps

        Args:
            q (float): charge's value
            new_pos_of_charge (map.Pos_Street): new position of the charge
            intersections_dists (list[float]): intersections' distances to the new_pos
            first_update (bool): this is the first update? (if yes, subtracting from parent isn't necessary).
                For first update, it's better using first_update function. The first_update function uses this
                function, that's why we have this bool. But it's not recomended to call this fction for first
                updates, because the fction can be changed.
        """
        self._pos = new_pos_of_charge
        if not first_update:
            old_map = Map_of_Potentials(self._mapa, self._order, False, self._pot_map, self._empty_street)

        # with intersections_dists, we can see the distances of the points on the street
        max_range = MAXIMAL_RANGE[self._order]
        streets = self._mapa.get_streets_list()
        charges_street = new_pos_of_charge.get_street().get_id()
        for i in range(self._n_streets):
            s = streets[i]
            (A,B) = s.get_intersections_ids()

            # is this the street where the charge is?
            if i == charges_street:
                # it will start from A to B. Therefore, it subtracts division on A
                # and sums it on B, always choosing the positive one as the real distance.
                dist_from_A = self._intersections_dists[A]
                dist_from_B = self._intersections_dists[B] - s.get_length()
                n_steps = self.get_n_steps(s)
                division_size = s.get_length() / (n_steps-1)
                dists_this_street = np.zeros(n_steps)
                for k in range(n_steps):
                    d = max(dist_from_A, dist_from_B)   # choose the positive one
                    if d < max_range:
                        if d < EPSILON_DIST:
                            d = EPSILON_DIST
                        dists_this_street[k] = q / d**self._order                    
                    dist_from_A -= division_size
                    dist_from_B += division_size
                self._pot_map[i] = dists_this_street
                self._empty_street[i] = False

            # it isn't the charge's street
            else:
                dist_from_A = self._intersections_dists[A]
                dist_from_B = self._intersections_dists[B]
                
                # if the entire street is out of the maximal range
                if dist_from_A > max_range or dist_from_B > max_range:
                    self._empty_street[i] = True
                
                # if the street isn't out of range    
                else:
                    # it will start from A, going to B. Therefore, it
                    # sums division on A and subtracts it from B, always
                    # choosing the smaller one as the real distance
                    dist_using_A = dist_from_A
                    dist_using_B = dist_from_B + s.get_length()
                    n_steps = self.get_n_steps(s)
                    division_size = s.get_length() / (n_steps-1)
                    dists_this_street = np.zeros(n_steps)
                    for k in range(n_steps):
                        d = min(dist_using_A, dist_from_B)  # chooses the smaller one
                        if d < max_range:
                            if d < EPSILON_DIST:
                                d = EPSILON_DIST
                            dists_this_street[k] = q / d**self._order
                        dist_using_A += division_size
                        dist_using_B -= division_size
                    self._empty_street[i] = False
                    self._pot_map = dists_this_street
        
        if first_update:
            self._parent_map.charge_asking_update(self)
        else:
            self._parent_map.charge_asking_update(self - old_map)
      
    def first_update(self, q: float, pos_of_charge: map.Pos_Street):
        self.update_charges_position(q, pos_of_charge, True)

class All_Pots_and_Charges:
    def __init__(self, mapa: map.Map) -> None:
        """Concentrates everything concerning potentials: every map of potential,
        every charge, potentials' orders, etc.
        You can use 2 objects of this class if you want separate trash generators from trash getters

        Args:
            mapa (map.Map): map of the city
        """
        self._mapa = mapa
        self._list_of_potentials: list[Potential] = []
        self._list_of_potentials_orders: list[int] = []
        self._maps_of_pot_one_order: list[Map_of_Pot_One_Order] = []
        self._create_potentials()

        self._list_of_charges: list[Pontual_Charge] = []
        self._list_of_charges_ids: list[int] = []
        self._n_charges: int = 0
        self._queue_of_ids_reopen: queue.Queue[int] = queue.Queue()

    def _create_potentials(self):
        """Creates the types of potential: V = kq/d^order, order = {1,2,3,...,MAXIMAL_ORDER}. Since k is used
        only for transforming mesure units, it considers k=1.

        Args:
            order (int): order from kq/d^order
        """
        # Don't need any updates on maps since we don't change anything on the charges
        for o in range(MAXIMAL_ORDER):
            aux = Potential(self._mapa, o+1)
            self._list_of_potentials.append(aux)
            self._list_of_potentials_orders.append(o+1)
            self._maps_of_pot_one_order.append(aux.get_map_of_pot())

        # n_pots = len(self._list_of_potentials_orders)

        # # if there isn't any potential yet, we can just add it
        # if n_pots == 0:
        #     aux = Potential(self._mapa, order)
        #     self._list_of_potentials.append(aux)
        #     self._list_of_potentials_orders.append(order)
        #     self._maps_of_pot_one_order.append(aux.get_map_of_pot())
        #     return
        
        # # if there is, we have to check if a potential with the same order doesn't already exists
        # index = fcs.search_in_vec(self._list_of_potentials_orders, order)
        # if index >= n_pots or self._list_of_potentials_orders[index] != order:
        #     aux = Potential(self._mapa, order)
        #     self._list_of_potentials.insert(index, aux)
        #     self._list_of_potentials_orders.insert(index, order)
        #     self._maps_of_pot_one_order.insert()

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
        
    def create_charge(self, pos: map.Pos_Street) -> int:
        """creates a pontual charge in the position 'pos'

        Args:
            pos (map.Pos_Street): position of the charge

        Returns:
            int: id of the charge created
        """
        # Since we don't put any q to this charge, we don't need to update maps
        id = self._next_unused_charges_id()
        index_to_insert = fcs.search_in_vec(self._list_of_charges_ids, id)
        
        self._list_of_charges.insert(index_to_insert, Pontual_Charge(pos, id, self._mapa, self._maps_of_pot_one_order))
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

    def remove_charge(self, id: int):
        """Deletes the charge 'id'

        Args:
            id (int): id of the charge to be deleted
        """
        index = self._get_index_of_charge(id)
        if index >= 0:
            charge = self._list_of_charges[index]
            charge.set_qs_to_zero()
            self._list_of_charges.pop(index)
            self._list_of_charges_ids.pop(index)
            self._queue_of_ids_reopen.put(id)

    def move_charge(self, id: int, new_pos: map.Pos_Street):
        """Moves the charge 'id' to a new position (new_pos) in the map.

        Args:
            id (int): id of the carge to be moved
            new_pos (map.Pos_Street): position to where charge goes
        """
        index = self._get_index_of_charge(id)
        if index < 0:
            print('ID not found')
            return
        self._list_of_charges[index].move_to(new_pos)

class Potential:
    def __init__(self, mapa: map.Map, order: int):
        """Represents the Potential of order 'order'. 
        We don't need to create a constant since it can be inside the value of the carges (like k*q = Q).

        Args:
            mapa (map.Map): map of the city
            order (int): order n of potential from V = k.q/d^n
        """
        self._mapa = mapa
        self._map_of_pot = Map_of_Pot_One_Order(mapa, order)
        self._order = order

    def get_map_of_pot(self) -> Map_of_Pot_One_Order:
        return self._map_of_pot
    
    def get_order(self) -> int:
        return self._order
        
class Pontual_Charge:
    def __init__(self, pos_street: map.Pos_Street, id: int, mapa: map.Map) -> None:
        """Model for each pontual charge, considering every type of charge.
        Therefore, this pontual charge has a specific value of charge for each order.

        Args:
            qs (float): charges' values for each order (put in same order, pls)
            orders (int): constants n from kq/d^n
            pos_street (map.Pos_Street): charge's position
            id (int): charge's ID (used to identify this 'object')
            n_divs_street (int): exemples:
                0: updates only intersections
                3: updates intersections and 0.25, 0.5, 0.75 in the pos_in_street
            maps_pot_parent (Map_of_Potentials): maps of the sums of every charge's potentials of order 'orders'
            mapa (map.Map): city's map

        Other private variables:
            pot_matrix (np.array): map of this charge's potential
            intersections_dists (np.array): distance to each intersection
        """
        self._pos_street = pos_street
        self._id = id
        self._mapa = mapa
        
        self._qs: list[float] = []
        self._orders: list[int] = []
        self._n_orders = 0
        self._maps_of_pot: list[Map_of_Pot_for_Charges] = []
        self._maps_of_pot_parent: list[Map_of_Pot_One_Order] = []

        # store the distance from charge to each intersection
        self._intersetions_dists = np.array([])
        aux = fcs.dijkstra(mapa, pos_street)
        n_intersecs = len(mapa.get_intersections_list())
        for i in range(n_intersecs):
            np.append(self._intersetions_dists, aux[i]['distance_to_origin'])

    def new_order_of_pot(self, q_for_new_order: float, new_order: int, map_of_pot_parent: Map_of_Pot_One_Order):
        for i in range(self._n_orders):
            # if already exists this order, it just updates q
            if self._orders[i] == new_order:
                self._maps_of_pot[i].update_charges_value(q_for_new_order / self._qs[i])
                self._qs[i] = q_for_new_order
                return

        # if doesn't exists, les create it!!!
        self._qs.append(q_for_new_order)
        self._orders.append(new_order)
        self._maps_of_pot_parent.append(map_of_pot_parent)

        # creating map of pot for this order
        aux = Map_of_Pot_for_Charges(self._mapa, new_order, map_of_pot_parent, self._intersetions_dists)
        aux.first_update(q_for_new_order, self._pos_street)
        self._maps_of_pot.append(aux)

    def set_q_for_only_one_pot(self, new_q: float, order: int) -> None:
        # try to find the index on the vectors
        for i in range(self._n_orders):
            if self._orders[i] == order:
                self._maps_of_pot[i].update_charges_value(new_q, self._qs[i])
                self._qs[i] = new_q
                return

        # if found: good. Ohterwise, it has to append this order,
        # but we don't have the map_of_pot_parent. So it will raise an error, for now #TODO #5
        raise ValueError('It doesnt exist this order')
    
    def set_qs(self, new_qs: list[float], orders: list[int]) -> None:
        for i in range(len(new_qs)):
            self.set_q_for_only_one_pot(new_qs[i], orders[i])

    def move_to(self, new_pos_street: map.Pos_Street):
        """It changes the position of the charge to new_pos_street, 
        updating the maps.

        Args:
            new_pos_street (map.Pos_Street): the position to which the charge is going
        """
        # change pos
        self._pos_street = new_pos_street
        
        # update the vector with distances to intersections
        aux = fcs.dijkstra(self._mapa, new_pos_street)
        for i in range(len(self._mapa.get_intersections_list())):
            self._intersetions_dists[i] = aux[i]['distance_to_origin']

        # update maps of potentials
        for i in range(len(self._orders)):
            self._maps_of_pot[i].update_charges_position(self._qs[i], new_pos_street)

    def get_id(self) -> int:
        return self._id

    def get_orders(self) -> list[int]:
        """returns all orders. The list is copied, so you can
        change it at will.

        Returns:
            list[int]: list of orders (copied, not original)
        """
        return self._orders.copy()

    def set_qs_to_zero(self):
        """Set all the charges to zero. Already updating maps.
        """
        for order in range(MAXIMAL_ORDER):
            self.set_q_for_only_one_pot(0, order)

    def _update(self):
        """DONT USE THIS FUNCION. ITS OLD. update the distances vector and the potentials map if needed

        Args:
            n_divisions_each_street (int): exemples:
                0: updates only intersections
                3: updates intersections and 0.25, 0.5, 0.75 in the pos_in_street
        """
        # TODO 3: this is a non-functioning function
        raise TypeError('Dont use this function, its old')
    
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
                        if d < EPSILON_DIST:
                            d = EPSILON_DIST
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
                        if d < EPSILON_DIST:
                            d = EPSILON_DIST
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
            if d > EPSILON_DIST:
                old_Kq = self._pot_matrix[0,0] * d**self._order
            else:
                d = self._intersections_distance[self._map_.get_street(0).get_intersections_ids()[1]]
                old_Kq = self._pot_matrix[0,(n_steps_in_street + 1)] * d**self._order

            # update local matrix and, afther, the global one
            self._pot_matrix = self._pot_matrix  * (self._K * self._q / old_Kq)
            self._need_update_potentials = False
            self._total_pot_matrix += self._pot_matrix

def Pot_Analyer(mapa: map.Map, all_pots_and_charges: All_Pots_and_Charges, everything: entities.Everything):
    bins_list = everything.get_bins_list()
    for bin in bins_list:
        bin.get_filling_rate() * bin.get_capacity()
