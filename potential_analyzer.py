import numpy as np

from bin_lib import *
from bin_lib.consts import *
import bin_lib.potentials as pot

class Trash_Reports:
    def __init__(self, mapa: map.Map) -> None:
        self._mapa = mapa
        n_streets = len(mapa.get_streets_list())
        self._reps_dirty_street = np.zeros([3,n_streets])
        self._reps_clean_street = np.zeros(n_streets)

    def new_rep(self, street: map.Street, str_sweeper: bool = True,
                time_since_cleaned_strs: int = 0, dirty_lvl: int = 2, clean = False):
        """for now, we will only use for street sweepers

        Args:
            street (map.Street): _description_
            str_sweeper (bool, optional): _description_. Defaults to True.
            time_since_cleaned_strs (int, optional): _description_. Defaults to 0.
            dirty_lvl (int, optional): how much dirty is it? 1 to 3 (low, medium or high). Defaults to 2
            clean (bool, optional): it's clean?. Defaults to False.
        """
        str_id = street.get_id()
        if not str_sweeper:
            raise ValueError("Reports only for street sweepers")
        if dirty_lvl < 1 or dirty_lvl > 3:
            raise ValueError("Dirty level should be between 1 and 3")
        if clean:
            self._reps_clean_street[str_id] += 1
        else:
            self._reps_dirty_street[dirty_lvl-1][str_id] += 1
        # if str_sweeper:
        #     weight = STR_SWEEPER_WEIGHT
        #     time_since_cleaned_strs = PERIOD_CLEAN_STREETS
        # else:
        #     weight = 1
        #     time_since_cleaned_strs = int(time_since_cleaned_strs / DAY + 0.5)
        # if clean:
        #     weight *= time_since_cleaned_strs
        # else
        #     weight *= (PERIOD_CLEAN_STREETS // DAY) - time_since_cleaned_strs

    def get_trash_lvl(self) -> np.array:
        """Based on reports, it returns the level of dirtiness in each street.
        The result is given between 0 and 1, where 0 is clean and 1 is very dirty.

        Returns:
            np.array: lvl of dirtiness by street between [0,1]
        """
        n_streets = len(self._mapa.get_streets_list())
        mean_lvl = np.zeros(n_streets)
        for i in range(n_streets):
            lvl1 = self._reps_dirty_street[0,i]
            lvl2 = self._reps_dirty_street[1,i]
            lvl3 = self._reps_dirty_street[2,i]
            clean = self._reps_clean_street[i]
            if lvl1 > 0 or lvl2 > 0 or lvl3 > 0:
                mean_lvl[i] = (lvl1 + 2*lvl2 + 3*lvl3) / (clean + lvl1 + lvl2 + lvl3) / 3
        return mean_lvl
    

def pot_analyzer(bins_pot: pot.All_Pots_and_Charges, coms_pot: pot.All_Pots_and_Charges,
                 reps: Trash_Reports)