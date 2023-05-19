from __future__ import annotations
import math as m
import numpy as np
import copy

from bin_lib import *
from bin_lib import map
import bin_lib.map as map
from bin_lib.consts import *

raise NotImplementedError("Don't use this for now")

class Heat_Flux_Map:
    def __init__(self, mapa: map.Map, flux_map: list[np.array] = [], same_value_street: list[bool] = []) -> None:
        self._mapa = mapa
        self._n_streets = len(mapa.get_streets_list())

        if len(flux_map) > 0:
            if len(flux_map) != self._n_streets or len(same_value_street) != self._n_streets:
                raise ValueError("flux_map and same_value_street should have length = number of strets")
            self._flux_map = copy.deepcopy(flux_map)
            self._same_value_street = copy.deepcopy(same_value_street)
        else:
            self._flux_map = []
            self._same_value_street = []
            for i in range(self._n_streets):
                self._flux_map.append(np.array[0])
                self._same_value_street.append(True)

    def __add__(self, y: Heat_Flux_Map) -> Heat_Flux_Map:
        if not isinstance(y, Heat_Flux_Map):
            raise TypeError("Sum available only between 2 Heat_Flux_Map")
        res = Heat_Flux_Map(self._mapa)
        for i in range(self._n_streets):
            both_one_value = (self._same_value_street[i] and y._same_value_street[i])
            sum_this_str = self._flux_map[i] + y._flux_map[i]
            # after sum, les see if the result didn't become one same value
            # in the whole street
            if not both_one_value:
                aux = sum_this_str[0]
                aux2 = sum_this_str - aux
                if np.allclose(aux2, 0, rtol=EPSILON_TEMPTR):
                    res._same_value_street[i] = True
                    res._flux_map[i] = np.array([aux])
                else:
                    res._same_value_street[i] = False
                    res._flux_map[i] = aux2
        return res
    
    def __radd__(self, y: Heat_Flux_Map) -> Heat_Flux_Map:
        return self.__add__(y)
    
    def __mul__(self, y: float) -> Heat_Flux_Map:
        if not isinstance(y, (float,int)):
            raise TypeError("Multiplication available only by constants (float/int)")
        
        res = Heat_Flux_Map(self._mapa)
        if y == 0 or abs(y) < EPSILON_TEMPTR**2:
            return res
        
        for i in range(self._n_streets):
            res._same_value_street[i] = self._same_value_street[i]
            res._flux_map[i] = self._flux_map[i] * y

        return res

    def __rmul__(self, y: float) -> Heat_Flux_Map:
        return self.__mul__(y)
        
class HFMap_One_Source(Heat_Flux_Map):
    def __init__(self, mapa: map.Map, source: Temptr_Source, gnd_source: Temptr_Source, total_map: HFMap_Total) -> None:
        super().__init__(mapa)
        self._source = source
        self._gnd = gnd_source
        self._total_map = total_map

    def refresh(self):
        if self._gnd == None:
            return

        # TODO: gnd in extremity
        if self._gnd.is_in_extremity():
            raise NotImplementedError("Ground shouldn't be in extremity. For now, it's not available.")
        
        if self._source.get_sources_div()          == self._gnd.get_sources_div() and \
           self._source.get_sources_str().get_id() == self._gnd.get_sources_str().get_id():
            old_HFMap = Heat_Flux_Map(self._mapa, self._flux_map, self._same_value_street)
            empty_map = Heat_Flux_Map(self._mapa)
            self._flux_map = empty_map._flux_map
            self._same_value_street = empty_map._same_value_street
            self._total_map.update(old_HFMap,self)
            return

        # initializing some useful variables
        n_intersecs = len(self._mapa.get_intersections_list())
        adj_matrix_as_array = self._mapa.get_adj_matrix().getA()
        aux = self._source.get_sources_str().get_vector()
        intersecA_source_id = aux[0].get_id()
        intersecB_source_id = aux[1].get_id()
        aux = self._gnd.get_sources_str().get_vector()
        intersecA_gnd_id = aux[0].get_id()
        intersecB_gnd_id = aux[1].get_id()
        old_HFMap = Heat_Flux_Map(self._mapa, self._flux_map, self._same_value_street)

        # First we need to remove the conexion between the
        # intersections A_B from source and ground.
        # But we should do it only if they are not in extremities.
        if not self._source.is_in_extremity():
            adj_matrix_as_array[intersecA_source_id][intersecB_source_id] = 0
            adj_matrix_as_array[intersecB_source_id][intersecA_source_id] = 0
        if not self._gnd.is_in_extremity():
            adj_matrix_as_array[intersecA_gnd_id][intersecB_gnd_id] = 0
            adj_matrix_as_array[intersecB_gnd_id][intersecA_gnd_id] = 0

        # If the source isn't in extremity, we created a new point,
        # so we should add the point to the matrix
        if not self._source.is_in_extremity():
            aux = np.zeros(n_intersecs)
            aux2 = (self._source.get_sources_div() *
                    get_size_of_div(self._source.get_sources_str()))
            aux[intersecA_source_id] = aux2
            aux[intersecB_source_id] = self._source.get_sources_str().get_length() - aux2
            np.column_stack((adj_matrix_as_array,aux))
            aux = np.append(aux,0)
            np.row_stack((adj_matrix_as_array,aux))
        
        # now, we need to make all the rows as negatives and
        # the main diagonal as the sum of his row/column, but positive
        for i in range(n_intersecs+self._source.is_in_extremity()):
            adj_matrix_as_array[i] = -abs(adj_matrix_as_array[i])
            adj_matrix_as_array[i][i] = -np.sum(adj_matrix_as_array[i])

        # But, the gnd's neighbours need is a little different.
        # The main diagonal needs to add the resistance between the
        # neighbour and the gnd
        aux = (self._gnd.get_sources_div() *
               get_size_of_div(self._gnd.get_sources_str()))
        adj_matrix_as_array[intersecA_gnd_id][intersecA_gnd_id] += aux
        adj_matrix_as_array[intersecB_gnd_id][intersecB_gnd_id] += \
                self._gnd.get_sources_str().get_length() - aux
        
        # to finish it, we need to add a row and a column for the
        # flux between the source and the gnd
        if not self._source.is_in_extremity():
            aux = np.zeros(n_intersecs+1)
            aux[-1] = 1
        else:
            aux = np.zeros(n_intersecs)
            if self._source.get_sources_div() == 0:
                aux[intersecA_source_id] = 1
            else:
                aux[intersecB_source_id] = 1
        np.column_stack((adj_matrix_as_array,aux))
        aux = np.append(aux,0)
        np.row_stack((adj_matrix_as_array,aux))

        # MATRIX IS GOOD!!!!!
        # now we need to intert it and get each intersections's temperature
        # plus the flux between source and gnd
        # equation will be: matrix * temps = res
        # therefore, we need: temps = matrix^(-1) * res
        matrix = np.matrix(adj_matrix_as_array)
        res = np.zeros(n_intersecs+1+self._source.is_in_extremity())
        res[-1] = self._source.get_temptr() - self._gnd.get_temptr()
        temps = matrix.getI() * res

        # with 'temps', it's trivial to calculate the fluxes
        sources_str_id = self._source.get_sources_str().get_id()
        gnds_str_id = self._gnd.get_sources_str().get_id()
        for i in range(self._n_streets):
            s = self._mapa.get_street(i)
            [A,B] = s.get_vector()
            
            if i == gnds_str_id:
                self._same_value_street[i] = False
                dist_from_A = self._gnd.get_sources_div() * get_size_of_div(s)
                flux_from_A = temps[A.get_id()] / dist_from_A
                flux_from_B = -temps[B.get_id()] / (s.get_length() - dist_from_A)
                n_divs = get_n_divs(s)
                aux = np.zeros(n_divs)
                gnds_div = self._gnd.get_sources_div()
                for div in range(n_divs):
                    aux[div] = flux_from_A if div < gnds_div else flux_from_B
                self._flux_map[i] = aux
            elif i == sources_str_id and not self._source.is_in_extremity():
                dT_source_gnd = self._source.get_temptr() - self._gnd.get_temptr()
                self._same_value_street[i] = False
                dist_from_A = self._gnd.get_sources_div() * get_size_of_div(s)
                flux_from_A = (temps[A.get_id()] - dT_source_gnd) / dist_from_A
                flux_from_B = (dT_source_gnd - temps[B.get_id()]) / (s.get_length() - dist_from_A)
                n_divs = get_n_divs(s)
                aux = np.zeros(n_divs)
                gnds_div = self._gnd.get_sources_div()
                for div in range(n_divs):
                    aux[div] = flux_from_A if div < gnds_div else flux_from_B
                self._flux_map[i] = aux
            else:
                self._same_value_street[i] = True
                self._flux_map[i] = np.array((temps[A.get_id()] - temps[B.get_id()]) / s.get_length())

        # update the total map
        self._total_map.update(old_HFMap, self)

    def set_new_gnd(self, new_gnd: Temptr_Source):
        self._gnd = new_gnd
        self.refresh()
    
class HFMap_Total(Heat_Flux_Map):
    def __init__(self, mapa: map.Map) -> None:
        super().__init__(mapa)

    def update(self, old_HFMap: Heat_Flux_Map, new_HFMap: Heat_Flux_Map):
        res: Heat_Flux_Map = self - new_HFMap + old_HFMap
        self._flux_map = res._flux_map
        self._same_value_street = res._same_value_street

    def get_flux_a_street(self, s: map.Street) -> dict['same_value': bool, 'flux': np.array]:
        i = s.get_id()
        dictio = {'same_value': self._same_value_street[i],
                  'flux': self._flux_map[i]}
        return dictio

class Temptr_Source:
    def __init__(self, mapa: map.Map, pos_str: map.Pos_Street, temptr: float, id: int,
                 heat_flux_map_total: HFMap_Total, gnd_source: Temptr_Source = None) -> None:
        self._temperature = temptr
        self._heat_flux_map = HFMap_One_Source(mapa, self, None, heat_flux_map_total)
        self.set_position(pos_str)
        if not gnd_source == None:
            self.set_new_gnd(gnd_source)
    
    def get_id(self) -> int:
        return self._id

    def set_temptr(self, new_temptr: float):
        self._temperature = new_temptr
    
    def get_temptr(self) -> float:
        return self._temperature
    
    def get_sources_div(self) -> int:
        return self._sources_div
    
    def get_sources_str(self) -> map.Street:
        return self._sources_str
    
    def is_in_extremity(self) -> bool:
        return self._extremity
    
    def set_position(self, new_pos_str: map.Pos_Street):
        self._sources_str = new_pos_str.get_street()
        self._pos_str = new_pos_str
        self._sources_div = get_div_from_pos_str(new_pos_str, point=True)
        self._extremity = (self._sources_div == 0 or
                           self._sources_div == get_n_divs(self._sources_str))
        
    def set_new_gnd(self, new_gnd: Temptr_Source):
        if new_gnd.get_id() == self._id:
            self._am_i_gnd = True
        self._heat_flux_map.set_new_gnd(new_gnd)

class Temptr_Map:
    def __init__(self, mapa: map.Map) -> None:
        self._mapa = mapa
        self._t_map: list[np.array] = []
        self._n_streets = len(mapa.get_streets_list())
        self._updated_streets = np.empty(self._n_streets, dtype=bool)

        for i in range(self._n_streets):
            s = mapa.get_street(i)
            self._t_map.append(np.zeros(get_n_divs(s)+1))

    def update(self, hf_map: HFMap_Total, gnd_source: Temptr_Source):
        self._updated_streets *= False
        s = gnd_source.get_sources_str()
        [A,B] = s.get_vector()
        aux = hf_map.get_flux_a_street(s)
        gnds_div = gnd_source.get_sources_div()
        resist = get_size_of_div(s)
        t_gnd = gnd_source.get_temptr()

            


    def _update_rec(self, hf_map: HFMap_Total, intersec: map.Intersection, temptr: float):
        pass

    def _equalize_intersec_temptr(self, intersec: map.Intersection, temptr: float):
        neighbors = intersec.get_neighbors()
        for neighbor in neighbors:
            s = self._mapa.search_for_a_street(intersec, neighbor)
            A = s.get_vector()[0]
            if A.get_id() == intersec.get_id():
                self._t_map[s.get_id()][0] = temptr
            else:
                self._t_map[s.get_id()][-1] = temptr


class Temperature_Handler:
    def __init__(self, mapa: map.Map) -> None:
        self._mapa = mapa
        self._temptr_map: list[np.array[float]] = []
        self._source_flags: list[np.array[bool]] = []
        self._sources_list: list[(int,int)] = []
        self._sources: list[Temptr_Source] = []
        self._n_sources: int = 0
        self._gnd_source_id: int = -1
        self._total_heat_flux = HFMap_Total(self._mapa)


        streets = mapa.get_streets_list()
        n_str = len(streets)
        for i in range(n_str):
            s = streets[i]
            n_divs = get_n_divs(s)
            aux = np.empty(n_divs, dtype=bool)
            aux *= False
            self._source_flags.append(aux)
            self._temptr_map.append(np.zeros(n_divs))

    def new_temptr_source(self, pos_str: map.Pos_Street, temptr: float):
        aux = Temptr_Source(self._mapa, pos_str, temptr, self._total_heat_flux)
        self._sources.append(aux)
        self._n_sources += 1
        if self._gnd_source_id == -1:
            if not aux.is_in_extremity():
                self._gnd_source_id = self._n_sources - 1
                self._set_gnd()
        else:
            aux.set_new_gnd(self._sources[self._gnd_source_id])

    def _find_new_gnd(self) -> int:
        aux = self._gnd_source_id
        if self._n_sources > 0:
            for i in range(self._n_sources):
                if not self._sources[i].is_in_extremity():
                    self._gnd_source_id = i
                    return i
        self._gnd_source_id = -1
        return -1
    
    def _set_gnd(self):
        if self._gnd_source_id < 0:
            return
        
        new_gnd = self._sources[self._gnd_source_id]
        for i in range(self._n_sources):
            self._sources[i].set_new_gnd(new_gnd)





    def set_temperature_pos_str(self, pos_str: map.Pos_Street, new_temp: float):
        """Set the temperature of a division in a street to 'new_temp'.
        The division is where the 'pos_str' is.

        Args:
            pos_str (map.Pos_Street): where to set the temperature
            new_temp (float): new value of temperature
        """
        s = pos_str.get_street()
        s_id = s.get_id()
        div = get_div_from_pos_str(pos_str, point=True)
        self._temptr_map[s_id][div] = new_temp
        if not self._source_flags[s_id][div]:
            self._source_flags[s_id][div] = True
            self._sources_list.append((s_id,div))

    def set_temperature_street(self, street: map.Street, new_temp: float):
        """Set the temperature of a whole street to a new value.

        Args:
            street (map.Street): street used in this functions
            new_temp (float): the new temperature for the whole street
        """
        n_divs = get_n_divs(street)
        s_id = street.get_id()
        self._temptr_map[s_id] = np.zeros(n_divs) + new_temp
        for k in range(n_divs):
            if not self._source_flags[s_id][k]:
                self._source_flags[s_id][k] = True
                self._sources_list.append((s_id,k))

    def get_temperature_pos_str(self, pos_str: map.Pos_Street) -> float:
        s = pos_str.get_street()
        div = get_div_from_pos_str(pos_str,point=True)
        return self._temptr_map[s.get_id()][div]

    def get_temperature_street(self, street: map.Street) -> np.array:
        """Get the array of the temperature in the street.
        CAREFUL: it's the original, not a copy.

        Args:
            street (map.Street): street whose temperature you want

        Returns:
            np.array: array of temperatures in the street
        """
        return self._temptr_map[street.get_id()]

    def heat_a_pos_str(self, pos_street: map.Pos_Street, dT: float):
        s = pos_street.get_street()
        s_id = s.get_id()
        div = get_div_from_pos_str(pos_street, point=True)
        self._temptr_map[s_id][div] += dT
        if not self._source_flags[s_id][div]:
            self._source_flags[s_id][div] = True
            self._sources_list.append((s_id,div))

    def heat_a_street(self, street: map.Street, dT: float):
        s_id = street.get_id()
        self._temptr_map[s_id] += dT
        n_divs = get_n_divs(street)
        for k in range(n_divs):
            if not self._source_flags[s_id][k]:
                self._source_flags[s_id][k] = True
                self._sources_list.append((s_id,k))
    
    def thermal_balance(self):
        streets = self._mapa.get_streets_list()
        for i in range(len(streets)):
            s = streets[i]
            


def get_n_divs(street: map.Street) -> int:
    return int(street.get_length() / SIZE_OF_DIV_TEMPTR)

def get_size_of_div(street: map.Street) -> float:
    return street.get_length() / get_n_divs(street)

def get_div_from_pos_str(pos_str: map.Pos_Street, point=False) -> int:
    """Given a position in a street, the function returns the division
    where the pos_str belongs. In case of point=True, the function will
    return the closest point (considering a margin of MAX_DIST_EXTREMITY
    for the extremities), not the division.

    Args:
        pos_str (map.Pos_Street): _description_
        point (bool, optional): _description_. Defaults to False.

    Returns:
        int: _description_
    """
    n_divs = get_n_divs(pos_str.get_street())

    if point:
        dist_A =   pos_str.get_pos_in_street()   * pos_str.get_street().get_length()
        dist_B = (1-pos_str.get_pos_in_street()) * pos_str.get_street().get_length()

        if dist_A < MAX_DIST_EXTREMITY:
            return 0
        if dist_B < MAX_DIST_EXTREMITY:
            return n_divs
        div = np.round(pos_str.get_pos_in_street() * n_divs)
        if div == 0:
            return 1
        if div == n_divs:
            return n_divs-1
        return div
    else:
        div = int(pos_str.get_pos_in_street() * n_divs)
        if div == n_divs:
            return n_divs-1
        return div

