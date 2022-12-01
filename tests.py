from bin_lib import *

file_intersections = "bin_lib/intersection_points.txt"
file_streets = "bin_lib/streets.txt"

def map_test(files: list[str]):
    mapa = map.read_map(files[0], files[1], "")
    intersections = mapa.get_intersection_list()
    streets = mapa.get_streets_list()
    for x in streets:
        x.get_length()
    some_functions.plot_map(intersections, streets)