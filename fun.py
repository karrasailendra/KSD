import numpy as np
import math
KP = 15
def perpendicular_distance(x1, y1, a, b, c):
    d = abs((a * x1 + b * y1 + c)) / (math.sqrt(a * a + b * b))
    return d


def lineFromPoints(ox, oy):  # you might face difficulty for 90 slope
    if ox[1] == ox[0]:
        a = 0
        b = 1
        c = - ox[1]
    else:
        m = (oy[1] - oy[0]) / (ox[1] - ox[0])
        a = 1
        b = -m
        c = -(oy[0] - ox[0] * m)

    return a, b, c

def calc_potential_field(gx, gy, ox, oy, reso, rr, area):

    xw = area
    yw = area
    minx = 0
    miny = 0
    print(ox)

    pmap = np.zeros([xw, yw])

    for ix in range(xw):
        print(ix)
        x = ix * reso + minx
        print(x)
        for iy in range(yw):
            y = iy * reso + miny
            ug = np.zeros(len(gx))
            uo = np.zeros(len(ox))
            for k in range(len(gx)):
                ug[k] = calc_attractive_potential(x, y, gx[k], gy[k])
            ug1 = sum(ug)
            #            breakpoint()
            for l in range(len(ox)):
                uo[l] = calc_repulsive_potential(x, y, ox[l], oy[l], rr)
            #            breakpoint()
            uo1 = sum(uo)
            uf = ug1 + uo1
            pmap[ix, iy] = uf

    return pmap, minx, miny


def calc_attractive_potential(x, y, gx, gy):
    if np.hypot(x - gx, y - gy) > 20:
        return 0
    elif np.hypot(x - gx, y - gy) == 0:
        return 200
    else:
        if 50 * (0.5 * KP * (1 / np.hypot(x - gx, y - gy))) > 200:
            return 200
        else:
            return 50 * (0.5 * KP * (1 / np.hypot(x - gx, y - gy)))


def calc_repulsive_potential(x, y, ox, oy, rr):

    a, b, c = lineFromPoints(ox, oy)
    minoy = min(ox)
    minox = min(oy)
    maxoy = max(ox)
    maxox = max(oy)
    dq = perpendicular_distance(x, y, a, b, c)
    if minox - 20 > x or maxox + 20 < x or minoy - 20 > y or maxoy + 20 < y:
        return 0.0
    else:
        # if dq <= rr:
        if dq <= 0.1:
            dq = 0.1
            return -200
            # if -50*(0.5 * KP * (1/dq))< -200:
            #     return -200
            # else:
        elif -50 * (0.5 * KP * (1 / dq)) < -200:
            return -200
        else:
            return -50 * (0.5 * KP * (1 / dq))
