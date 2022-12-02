from __future__ import annotations
import matplotlib.pyplot as plt
import bin_lib.map as map
import bin_lib.entities as entities
import bin_lib.some_functions as fcs

NUMBER_OF_SIMULATIONS = 3
TIME_OF_SIMULATION = 3600*365*5
TIME_STEP = 60
FILE_INTERSECTIONS = "bin_lib/intersections.txt"
FILE_STRETS = "bin_lib/streets.txt"
FILE_COMMERTIAL_POINTS = ""

# start: read map
mapa = map.read_map(FILE_INTERSECTIONS, FILE_STRETS, FILE_COMMERTIAL_POINTS)
# fcs.plot_map(mapa)
# com_points = fcs.create_rand_com_points(mapa, (5,2,4))
# fcs.plot_entities(com_points)
# plt.show()
everything = entities.Everything(mapa)
fcs.create_rand_com_points(mapa, (5,2,4), everything)
fcs.create_rand_bins(mapa, everything)

NUMBER_OF_STEPS = TIME_OF_SIMULATION // TIME_STEP

for sims in range(NUMBER_OF_SIMULATIONS):
    for t in range(NUMBER_OF_STEPS):
        everything.update_people(TIME_STEP)
        fcs.create_rand_ppl(mapa, everything, TIME_STEP)

