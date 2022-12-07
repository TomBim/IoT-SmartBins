from __future__ import annotations
import random as rand
import numpy as np
import math as m
from math import sqrt, inf
import heapq
import matplotlib as mpl
import matplotlib.pyplot as plt
import bin_lib.map as map
import bin_lib.entities as entities
from bin_lib.consts import *

def dijkstra(mapa: map.Map, origin: map.Pos_Street):
    intersections_info = [{'distance_to_origin': inf, 'parent': None, 'closed': False} for _ in mapa.get_intersections_list()]
    pq: list[tuple[float, map.Intersection]] = []
    starting_intersections = origin.get_street().get_vector()

    for intersection in starting_intersections:
        distance = calculate_distance(origin.get_pos_xy(), intersection.get_pos())
        intersections_info[intersection.get_id()]['distance_to_origin'] = distance
        heapq.heappush(pq, (distance, intersection))

    while len(pq) > 0:
        (distance, intersection) = heapq.heappop(pq)
        id = intersection.get_id()
        if not intersections_info[id]['closed']:
            intersections_info[id]['closed'] = True
            neighbors = intersection.get_neighbors()
            neighbors_dists = intersection.get_distances_to_neighbors()
            for neighbor, neighbor_dist in zip(neighbors, neighbors_dists):
                neighbor_id = neighbor.get_id()
                if intersections_info[neighbor_id]['distance_to_origin'] > distance + neighbor_dist:
                    neighbor_dist = distance + neighbor_dist
                    intersections_info[neighbor_id]['distance_to_origin'] = neighbor_dist
                    intersections_info[neighbor_id]['parent'] = id
                    heapq.heappush(pq, (neighbor_dist, neighbor))
    return intersections_info

def calculate_line_equation(point1: tuple[float, float], point2: tuple[float, float]) -> tuple[float, float, float]:
    """calculates m, n, r:
    mx + ny = r

    Args:
        point1 (tuple[float]): _description_
        point2 (tuple[float]): _description_

    Returns:
        list[float, float, float]: [m,n,r]
    """
    (x1,y1) = point1
    (x2,y2) = point2
    deltax = x2-x1
    deltay = y2-y1
    # if vertical line
    if deltax == 0:
        return (1,0,x1)
    # if horizontal line
    if deltay == 0:
        return (0,1,y1)
    # normal line
    return (1/deltax, -1/deltay, (x1/deltax - y1/deltay))

def calculate_distance(point1: tuple[float, float], point2: tuple[float, float]) -> float:
    return sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

def search_in_vec_aux(vec: list, value: int, i_min: int, i_max: int) -> int:
    i = int((i_max + i_min) / 2)
    aux = vec[i]
    if aux == value:
        return i
    if i_min == i_max:
        if aux < value:
            return i_min+1
        return i_min
    if value > aux:
        return search_in_vec_aux(vec, value, i+1, i_max)
    return search_in_vec_aux(vec, value, i_min, i)

def search_in_vec(vec: list, value: int) -> int:
    """search a value in a vector and returns
    the index. The vector must be sorted already.
    If value isn't found, it will return the index 
    of the first greater id, so you can insert xD.

    Args:
        vec (list): vector in which we will make the search
        valeu (int): value searched

    Returns:
        int: index of the value searched
    """
    if len(vec) == 0:
        return 0
    return search_in_vec_aux(vec, value, 0, len(vec)-1)

def plot_map(mapa: map.Map) -> None:
    plt.figure(figsize=(5,2.7), layout='constrained')
  
    # plot streets
    streets = mapa.get_streets_list()
    for s in streets:
        A, B = s.get_vector()
        A = A.get_pos()
        B = B.get_pos()
        x = np.array([A[0], B[0]])
        y = np.array([A[1], B[1]])
        plt.plot(x,y, label="")

    # plot intersections
    intersections = mapa.get_intersections_list()
    x = np.array([])
    y = np.array([])
    for i in intersections:
        pos = i.get_pos()
        x = np.append(x, pos[0])
        y = np.append(y, pos[1])
    plt.plot(x,y, '.', markersize=10, color='k')
        
def plot_entities(com_points: list[entities.Commercial_Point]=None, bins: list[entities.Bin]=None, people: list[entities.Person]=None):
    if not com_points == None:
        xf = np.array([])
        yf = np.array([])
        xnf = np.array([])
        ynf = np.array([])
        xj = np.array([])
        yj = np.array([])
        for c in com_points:
            t = c.get_type()
            pos = c.get_pos()
            if t == 0:
                xf = np.append(xf, pos[0])
                yf = np.append(yf, pos[1])
            elif t == 1:
                xnf = np.append(xnf, pos[0])
                ynf = np.append(ynf, pos[1])
            else:
                xj = np.append(xj, pos[0])
                yj = np.append(yj, pos[1])                
        plt.plot(xf, yf,'x', color='b', label="Restaurantes", markersize=15)
        plt.plot(xnf, ynf,'*', color='b', label="Lojas", markersize=10)
        plt.plot(xj, yj,'^', color='b', label="Industrias", markersize=10)


    if not bins == None:
        x = np.array([])
        y = np.array([])
        for b in bins:
            pos = b.get_pos()
            x = np.append(x, pos[0])
            y = np.append(y, pos[1])
        plt.plot(x,y,'*', 'r', label="Lixeiras")

    if not people == None:
        x = np.array([])
        y = np.array([])
        for p in people:
            pos = p.get_pos()
            x = np.append(x, pos[0])
            y = np.append(y, pos[1])
        plt.plot(x,y,'*', 'j', label="Pessoas")
    
    plt.legend()

def create_rand_points(n_points: int, max_range: float) -> list[tuple[float]]:
    p = []
    for n in range(n_points):
        p.append((rand.uniform(-max_range,max_range),rand.uniform(-max_range,max_range)))
    return p

def create_rand_streets_aux(points: list [tuple[float]], s: list[tuple[int]], s0: tuple[int]) -> bool:
    """see if the street crosses with some other stret
    obs:
        we will use equations for the lines like below:
            px + qy = r
    Args:
        s (list[tuple[int]]): _description_
        s0 (tuple[int]]): _description_
    """

    L = len(s)

    # build equation of the line s0
    a0 = points[s0[0]]
    b0 = points[s0[1]]
    a0x = a0[0]
    a0y = a0[1]
    b0x = b0[0]
    b0y = b0[1]
    if a0x == b0x:
        p0 = 1
        q0 = 0
        r0 = a0x
    elif a0y == b0y:
        p0 = 0
        q0 = 1
        r0 = a0y
    else:
        p0 = 1 / (b0x - a0x)
        q0 = - 1 / (b0y - a0y)
        r0 = a0x/(b0x - a0x) - a0y/(b0y - a0y)

    # for each street, see if s0 crosses it
    for i in range(L):
        st = s[i]

        # if they have extremity in common, they don't cross each other
        if not (st[0] == s0[0] or st[0] == s0[1] or st[1] == s0[0] or st[1] == s0[1]):
            a = points[st[0]]
            b = points[st[1]]
            ax = a[0]
            ay = a[1]
            bx = b[0]
            by = b[1]

            # if both points from st are above/below s0, they don't cross each other
            ka = p0*ax + q0*ay - r0
            kb = p0*bx + q0*by - r0
            if not ka*kb > 0:                
                # build equation of the line st
                if ax == bx:
                    p = 1
                    q = 0
                    r = ax
                elif ay == by:
                    p = 0
                    q = 1
                    r = ay
                else:
                    p = 1 / (bx - ax)
                    q = - 1 / (by - ay)
                    r = ax/(bx - ax) - ay/(by - ay)

                # if both points from s0 are above/below st, they don't cross each other
                # otherwise, they will cross
                ka = p*a0x + q*a0y - r
                kb = p*b0x + q*b0y - r
                if not ka*kb > 0:
                    return True
    return False
        
def create_rand_streets(points: list[tuple[float]]) -> list[tuple[int]]:
    s: list[tuple[int]] = []
    leng = len(points)
    # n_streets won't necessarily be the number of streets. It will be just a guide number.
    # We take care of the case in which it's not possible to create all these streets after.
    n_streets = rand.randrange(int((leng*(leng-1))/2))
    print(f"previsto: {n_streets}")
    if n_streets > 0:        
        tried_matrix = np.zeros((leng,leng))
        for i in range(leng):
            tried_matrix[i][i] = 1
        number_available = leng*(leng-1)
        for n in range(n_streets):
            a = rand.randrange(leng)
            b = rand.randrange(leng)
            tried = tried_matrix[a][b]
            while (number_available > 0) and (a == b or tried or create_rand_streets_aux(points, s, (a,b))):
                if not tried:
                        tried_matrix[a][b] = 1
                        tried_matrix[b][a] = 1
                        number_available -= 2
                a = rand.randrange(leng)
                b = rand.randrange(leng)
                tried = tried_matrix[a][b]
            if number_available > 0:
                s.append((a,b))
                tried_matrix[a][b] = 1
                tried_matrix[b][a] = 1
                number_available -= 2
            else:
                print(f"only: {n}")
        print(f"n available: {number_available}")
        s.sort()
    return s

# TODO: see if these gaussians are good
def create_rand_com_points(mapa: map.Map, weight_for_com_points: tuple[int], everything: entities.Everything):
    """Generate random commercial points

    Args:
        mapa (map.Map): map
        weight_for_com_points (tuple[int]): weight for types
            [0]: food
            [1]: non-food
            [2]: job

    Returns:
        list[entities.Commercial_Point]: list of the randomly generated commercial points
    """
    streets = mapa.get_streets_list()
    n_streets = len(streets)
    n_com_points = 2 + rand.randrange(n_streets*3)
    random_types = rand.choices((0,1,2), weights=weight_for_com_points, k=n_com_points)
    for i in range(n_com_points):
        # pick a type and generate random trash potential for the commertial point
        random_types_i = random_types[i]
        payload = []
        if random_types_i == 0: # food
            prob_trash = m.sqrt(rand.random()) # higher probability, in general
            payload = {
                "mu_vol": 1e-3 + 2*rand.random(), # between 1mL and 2L
                "sigma_vol": 0.3,
                "mu_time_1": 30 + 150*rand.random(), # between 30s and 3min
                "sigmas_time": (0.3, rand.randrange(60,99)/100), # 0.9 and 0.60-0.99
                "weights": (rand.random(), rand.random()) # rly no idea if its good
            }
            trash_pot = entities.Trash_Potential(prob_trash,payload)
        elif random_types_i == 1: # non-food
            prob_trash = (rand.random())**2 # smaller probability, in general
            payload = {
                "mu_vol": 1e-3 + 2*rand.random(), # between 1mL and 2L
                "sigma_vol": 0.8,
                "mu_time_1": 0, # for non-food, the gaussian with non-zero mean isn't important
                "sigmas_time": (0.9, 1), # 0.9
                "weights": (1,0) # gaussian with mean != 0 isn't necessary here
            }
            trash_pot = entities.Trash_Potential(prob_trash,payload)
        else: # job
            trash_pot = entities.Trash_Potential(prob_trash=0)
        
        # generate random position
        street = streets[rand.randrange(n_streets)]
        rand_pos_street = map.Pos_Street(street, EPSILON + (1-2*EPSILON)*rand.random())

        # put in the vector
        customer_potential = rand.random()

        everything.new_com_point(random_types_i,customer_potential,trash_pot,rand_pos_street)
 
def generate_random_map(n_points: int, max_range: float, file_intersections: str, file_streets: str, weight_for_com_points: tuple[int], file_com_points: str):
    """generates the files of random intersections, streets, and commertial points

    Args:
        n_points (int): number of intersection points
        max_range (float): maximal size of the map
            map can go from (-max_range,-max_range) to (max_range,max_range)
        file_intersections (str): name of file where the function will put the random intersection points
        file_streets (str): name of file where the function will put the random streets
        weight_for_com_points (tuple[int]): weight considered while generating the type of commertial points
        file_com_points (str): name of file where the function will put the random commertial points
    """
    p = create_rand_points(n_points, max_range)
    s = create_rand_streets(p)

    f = open(file_intersections, "w")
    mat = open("intersections.m", "w")
    mat.write("intersections = [")
    for point in p:
        f.write(str(point[0])+",")
        f.write(str(point[1])+"\n")
        mat.write(str(point[0])+",")
        mat.write(str(point[1])+";\n")
    f.close()
    mat.write("]")
    mat.close()

    f = open(file_streets,"w")
    mat = open("streets.m","w")
    mat.write("streets = [")
    for street in s:
        f.write(str(street[0])+",")
        f.write(str(street[1])+"\n")
        mat.write(str(street[0])+",")
        mat.write(str(street[1])+";\n")
    f.close
    mat.write("]")
    mat.close()

def create_rand_bins(mapa: map.Map, everything: entities.Everything):
    """Generate random bin locations

    Args:
        mapa (map.Map): map

    Returns:
        list[entities.Bin]: list of the randomly generated bins
    """
    streets = mapa.get_streets_list()
    n_streets = len(streets)
    n_bins = n_streets // 3 + rand.randrange(n_streets*2)

    bins_pos_xy: list[tuple[float, float]] = []
    
    bin_capacity = 20
    street = streets[rand.randrange(n_streets)]
    rand_pos_street = map.Pos_Street(street, rand.random())
    everything.new_bin(bin_capacity, rand_pos_street)
    bins_pos_xy.append(rand_pos_street.get_pos_xy())
    n = 1
    tries = 0
    while n < n_bins and tries < 3*n_bins:
        # generate random position
        street = streets[rand.randrange(n_streets)]
        rand_pos_street = map.Pos_Street(street, rand.random())

        # generate fixed / random bin capacity 
        # bin_capacity = 20
        # bin_capacity = 50 * rand.random()
        pos_xy = rand_pos_street.get_pos_xy()
        # put bin in the vector if there is no one else too close
        good = True
        L = len(bins_pos_xy)
        k = 0
        while good and k < L:
            if calculate_distance(bins_pos_xy[k], pos_xy) < 50:
                good = False
            k += 1
        if good:
            everything.new_bin(bin_capacity, rand_pos_street)
            tries = -1
            n += 1
            bins_pos_xy.append(pos_xy)
        tries += 1

def create_rand_ppl(mapa: map.Map, everything: entities.Everything, TIME_STEP: int):
    """create random ppl using poisson for getting the number of ppl that will be generated
    Obs:
        Should be called every step

    Args:
        mapa (map.Map): _description_
        everything (entities.Everything): _description_
        TIME_STEP (int): constant
    """
    attractiveness_of_com_points = everything.get_com_points_attractiveness()
    n_com_points = len(attractiveness_of_com_points)

    # lambda for poisson: depends on number of com. points and their attractiveness
    lambida = (n_com_points * attractiveness_of_com_points.mean())/150
    n_ppl = np.random.poisson(lambida)

    # create these ppl
    for i in range(n_ppl):
        everything.new_person()