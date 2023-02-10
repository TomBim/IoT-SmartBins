from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import bin_lib.map as map
import bin_lib.entities as entities
import bin_lib.some_functions as fcs
from bin_lib.consts import *
from solution import *
import potentials

NUMBER_OF_SIMULATIONS = 1
TIME_OF_SIMULATION = 3*DAY
TIME_STEP = 60
FILE_INTERSECTIONS = "bin_lib/intersections.txt"
FILE_STRETS = "bin_lib/streets.txt"
FILE_COMMERTIAL_POINTS = ""

# start: read map
mapa = map.read_map(FILE_INTERSECTIONS, FILE_STRETS, FILE_COMMERTIAL_POINTS)
# print(map.validate_a_map(mapa))
# fcs.plot_map(mapa)
# com_points = fcs.create_rand_com_points(mapa, (5,2,4), mapa)
# fcs.plot_entities(com_points)
# plt.show()
everything = entities.Everything(mapa)
fcs.create_rand_com_points(mapa, (5,2,4), everything)
fcs.create_rand_bins(mapa, everything)

# create potentials
com_points_1st = potentials.Potential(mapa, 3, 100, 1)
com_points_2nd = potentials.Potential(mapa, 5, 1e5, 2)
bins_2nd = potentials.Potential(mapa, 5, 1e6, 2)

# create charges
for cp in everything.get_com_points():
    com_points_1st.new_charge(10, cp.get_pos_street(), mapa, cp.get_id())
    com_points_2nd.new_charge(100, cp.get_pos_street(), mapa, cp.get_id())
for b in everything.get_bins_list():
    bins_2nd.new_charge(10, b.get_pos_street(), mapa, b.get_id())

# sei la
filling_rate_average = np.array([])


NUMBER_OF_STEPS = TIME_OF_SIMULATION // TIME_STEP

# x = np.array([])
# v = np.array([])
# bin0 = np.array([])
# bin1 = np.array([])
# bin2 = np.array([])
for sims in range(NUMBER_OF_SIMULATIONS):
    filling_rate_average = np.zeros(len(everything.get_bins_list()))
    for t in range(NUMBER_OF_STEPS):
        time_now = t*TIME_STEP
        everything.update_people(TIME_STEP, time_now)
        fcs.create_rand_ppl(mapa, everything, TIME_STEP)

        if ((int)(time_now) % TIME_TO_CLEAN_STREETS) == 0:
            everything.sweep_streets()
        
        if ((int)(time_now) % TIME_TO_EMPTY_BINS) == 0:
            everything.empty_bins(time_now)
            
            for i in range(1,len(everything.get_bins_list())):
                filling_rate_average[i] = (filling_rate_average[i] + everything.get_bins_list()[i].get_filling_rate()) / 2

    final_potential_map = assert_bins_pot(everything, filling_rate_average, com_points_1st, com_points_2nd, bins_2nd)
    realocate_bins(everything, final_potential_map)
        
        # print(len(everything._ppl))
        # x = np.append(x, time_now/DAY)
        # v = np.append(v, everything.get_trash_in_the_street())
        # bin0 = np.append(bin0, everything.get_bin(0).get_vol_trash())
        # bin1 = np.append(bin1, everything.get_bin(1).get_vol_trash())
        # bin2 = np.append(bin2, everything.get_bin(2).get_vol_trash())




# print(len(everything.get_bins_list()))
# print(everything._last_bins_id)
# print(len(everything.get_com_points()))

# plt.figure(figsize=(5,2.7), layout='constrained')
# plt.scatter(x,v,s=12, label='chao')
# plt.scatter(x,bin0,s=12, label='bin0')
# plt.scatter(x,bin1,s=12, label='bin1')
# plt.scatter(x,bin2,s=12, label='bin2')
# plt.xlabel("Days")
# plt.legend()
# plt.grid(visible=True)
# plt.show()