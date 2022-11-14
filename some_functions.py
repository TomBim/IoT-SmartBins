import random as rand
import numpy as np

def create_rand_points(n_points: int, max_distance: float) -> list[tuple[float]]:
    p = []
    for n in range(n_points):
        p.append((rand.uniform(-max_distance,max_distance),rand.uniform(-max_distance,max_distance)))
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

def generate_random_map(n_points: int, max_distance: float, file_intersections: str, file_streets: str):
    p = create_rand_points(n_points, max_distance)
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

generate_random_map(10,1000,"intersections.txt", "streets.txt")
print("ok")