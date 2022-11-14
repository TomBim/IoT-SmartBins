import map
import entities

file_intersections = "intersection_points.txt"
file_streets = "streets.txt"

def map_test(files: list(str)):
    mapa: map.Map = map.read_map(files[0], files[1])
    streets_list = mapa.get_intersection_list
    for x in streets_list:
        x.get_length