import numpy as np
import matplotlib.pyplot as plt
import random
import math

# Parameters
KP = 15  # attractive potential gain
ETA = 100.0  # repulsive potential gain
AREA_WIDTH = 30.0  # potential area width [m]
from sklearn.cluster import DBSCAN
# import sys
#
# sys.path.append("")

# from file import function
show_animation = True
import fun

def levy(x, y):
    z = np.zeros([1, 2])
    beta = 1.6
    num = math.gamma(1 + beta) * np.sin(np.pi * beta / 2)
    den = math.gamma((1 + beta) / 2) * beta * 2 ** ((beta - 1) / 2)
    sigma_u = (num / den) ** (1 / beta)
    k = 0  # random.randint(-50, 50)
    u = np.random.normal(k, sigma_u ** 2, 1)
    v = np.random.normal(k, 1, 1)
    z[0][0] = u[0] / (abs(v[0]) ** (1 / beta))

    theta = np.pi * random.random()
    x = x + 30 * z[0][0] * math.cos(theta)
    y = y + 30 * z[0][0] * math.sin(theta)
    out = [[x, y]]
    return out

def levyheading():

    beta = 1.6
    num = math.gamma(1 + beta) * np.sin(np.pi * beta / 2)
    den = math.gamma((1 + beta) / 2) * beta * 2 ** ((beta - 1) / 2)
    sigma_u = (num / den) ** (1 / beta)
    k = 0  # random.randint(-50, 50)
    u = np.random.normal(k, sigma_u ** 2, 1)
    v = np.random.normal(k, 1, 1)
    step_size = abs(120*u[0] / (abs(v[0]) ** (1 / beta)))
    randnum = [-1, 1]
    theta = randnum[random.randint(0, 1)]*np.pi * random.random()
    print(theta)
    if theta < 0:
        theta = 2 * np.pi + theta

    return theta, step_size

def get_motion_model():
    # dx, dy
    motion = [[1, 0],
              [0, 1],
              [-1, 0],
              [0, -1],
              [-1, -1],
              [-1, 1],
              [1, -1],
              [1, 1]]

    return motion


def opp(gx, gy, rx, ry,j):
    theta = math.atan2(gy - ry, gx - rx)
    dist = np.hypot(gy - ry, gx - rx)
    out = [[rx + dist * np.sin(theta), ry - dist * np.cos(theta)]]
    return out

def repulsive(distance, heading, j, ix, iy, nbot):
    restheta = 0
    near = 0
    for k in range(nbot):
        if distance[j][k] < 50 and [j] != [k]:
            theta1 = heading[k] - heading[j]
            if theta1 < 0:
                theta1 = 2 * np.pi + theta1
            restheta = restheta + theta1
            near = near + 1
            # ix = ix + factor*distance[j][k] * math.cos(theta1)

    out = 0        # iy = iy + factor*distance[j][k] * math.sin(theta1)
    if near > 0:
        restheta = heading[j] + (restheta / near) + np.pi/2
        if restheta < 0:
            restheta = 2 * np.pi + restheta
        out = restheta

    return out
def resheading(heading,repulsion,levy,ix,iy):
    heading_factor = 0.5
    repulsion_factor = 0.25
    levy_factor = 0.25
    if repulsion == 0:
        repulsion_factor = 0

    rescostheta = (heading_factor*math.cos(heading)+ repulsion_factor*math.cos(repulsion)+ levy_factor*math.cos(levy))
    ressintheta = (heading_factor*math.sin(heading)+ repulsion_factor*math.sin(repulsion)+ levy_factor*math.sin(levy))
    resheading = math.atan2(ressintheta, rescostheta)
    ix = ix + 1 * rescostheta
    iy = iy + 1 * ressintheta
    return ix, iy, resheading

def repulsive_new(distance, heading, j, ix, iy, nbot,rx,ry):
    costheta = 0
    sintheta = 0
    near = 0
    total_dist = 0
    for k in range(nbot):
        if distance[j][k] < 100 and [j] != [k]:
            theta = math.atan2(ry[j]-ry[k], rx[j]-rx[k])
            costheta = (1 / distance[j][k]) * (math.cos(theta)) + costheta
            sintheta = (1 / distance[j][k]) * (math.sin(theta)) + sintheta
            total_dist = (1 / distance[j][k]) + total_dist
            near = near + 1
            # ix = ix + factor*distance[j][k] * math.cos(theta1)

    out = 0        # iy = iy + factor*distance[j][k] * math.sin(theta1)
    if near > 0:
        costheta = ((costheta / total_dist))
        sintheta = ((sintheta / total_dist))
        restheta = math.atan2(sintheta,costheta)
    else:
        # costheta = np.pi/2
        # sintheta = 0
        restheta = heading[j]

    return restheta

def global_variable(distance, j,rpos,grp):
    costheta = 0
    sintheta = 0
    total_dist = 0
    for k in [h for h in grp if h != j]:
        # print(k)
        theta = math.atan2(rpos[j][1]-rpos[k][1], rpos[j][0]-rpos[k][0])
        costheta = (1 / distance[j][k]) * (math.cos(theta)) + costheta
        sintheta = (1 / distance[j][k]) * (math.sin(theta)) + sintheta
        total_dist = (1 / distance[j][k]) + total_dist

    costheta = ((costheta / total_dist))
    sintheta = ((sintheta / total_dist))
    restheta = math.atan2(sintheta, costheta)

    return restheta

class motionplanning:

    def __init__(self, ox, oy, reso, rr, motion, rx, ry, pmap, j, nbot, heading, iter, rpos, grp):
        # global robgoal
        global levy_heading
        global step_size
        global current_distance
        global bot

        ix = rx[j]
        iy = ry[j]
        factor = 0.0
        # breakpoint()
        # if np.hypot(goalx - ix, goaly - iy) > 2:
        globalheading = np.zeros(nbot)
        if iter != 0:
            for ki in range(len(grp)):
                for kj in grp[ki][0]:
                    globalheading[kj] = global_variable(distance, kj, rpos, grp[ki][0])

        if pmap[int(ix), int(iy)] != 0:
            # breakpoint()
            if pmap[int(ix), int(iy)] < 0:
                self.heading = heading[j] + np.pi/2
                self.ix = ix + 2*math.cos(self.heading)
                self.iy = iy + 2*math.sin(self.heading)
                if pmap[int(ix), int(iy)] >= pmap[int(self.ix), int(self.iy)]:
                    self.heading = heading[j] - np.pi / 2
                    self.ix = ix + 2*math.cos(self.heading)
                    self.iy = iy + 2*math.sin(self.heading)
                # if pmap[int(self.ix), int(self.iy)] < 0:
                #     self.heading = heading[j] - np.pi / 2
                #     self.ix = ix + 10 * math.cos(self.heading)
                #     self.iy = iy + 10 * math.sin(self.heading)
            else:
                maxp = -1 * float("inf")
                minp = float("inf")
                maxix, maxiy = 1, -1
                for i, _ in enumerate(motion):
                    inx = int(ix + motion[i][0])
                    iny = int(iy + motion[i][1])
                    if inx >= len(pmap) or iny >= len(pmap[0]):
                        p = -1 * float("inf")  # outside area
                    else:
                        p = pmap[inx, iny]
                    if minp > p:
                        minp = p
                    if maxp < p:
                        maxp = p
                        maxix = inx
                        maxiy = iny

                self.ix = maxix
                self.iy = maxiy
                self.heading = heading[j]+np.pi/2
                while self.heading > 2*np.pi:
                    self.heading = self.heading - 2*np.pi
                current_distance[j] = current_distance[j] + 2

        elif abs(current_distance[j] - step_size[j]) < 3:

            levy_heading[j], step_size[j] = levyheading()
            # robgoal[j] = levy(ix, iy)[0]
            # self.ix = ix
            # self.iy = iy
            self.heading = levy_heading[j]
            self.ix = ix + 1 * math.cos(self.heading)
            self.iy = iy + 1 * math.sin(self.heading)
            # print(levy_heading[j])
            # print("next")
            current_distance[j] = 0

        else:
            # theta = levy_heading
            # self.heading = 0* heading[j] + 1 * repulsive_new(distance, heading, j, ix, iy, nbot, rx, ry) + levy_heading[j]
            # if self.heading < 0:
            #     self.heading = 2*np.pi+self.heading
            # while self.heading > 2*np.pi:
            #     self.heading = self.heading - 2*np.pi
            if iter == 0:
                self.ix, self.iy, self.heading = resheading(heading[j], repulsive_new(distance, heading, j, ix, iy, nbot,rx,ry), levy_heading[j], ix, iy)
            else:
                self.ix, self.iy, self.heading = resheading(heading[j], globalheading[j], levy_heading[j], ix, iy)
                print(levy_heading)

            current_distance[j] = current_distance[j] + 0.5

# def clusgoal(grp, rpos):
#     global bot
#     for j in grp:
#         bot[j].globalheading = global_variable(distance,j,rpos,grp)


def potential_field_planning(sx, sy, gx, gy, ox, oy, reso, rr, nbot, area):
    # calc potential field
    global rmap
    global bmap
    global stepsize_mat
    global distance
    global levy_heading
    global step_size
    global current_distance
    global bot

    pmap, minx, miny = fun.calc_potential_field(gx, gy, ox, oy, reso, rr, area)
    # breakpoint()

    # search path

    if show_animation:
        draw_heatmap(pmap)
        plt.show

    robgoal = np.zeros([nbot, 2])
    step_size = np.zeros(nbot)
    levy_heading = np.zeros(nbot)

    # intialization
    bot = np.zeros(nbot)
    bot = list(bot)
    ix, iy, rx, ry = np.zeros(nbot), np.zeros(nbot), np.zeros(nbot), np.zeros(nbot)
    heading = np.zeros(nbot)
    for rb in range(nbot):
        ix[rb] = round((sx[rb] - minx) / reso)
        iy[rb] = round((sy[rb] - miny) / reso)
        rx[rb] = sx[rb]
        ry[rb] = sy[rb]

    motion = get_motion_model()
    for i in range(nbot):
        robgoal[i] = levy(rx[i], ry[i])[0]
        levy_heading[i], step_size[i] = levyheading()
    stepsize_mat = []
    iter = 1000
    current_distance = np.zeros([nbot])
    rpos = np.zeros([nbot,2])
    distance = np.zeros([nbot,nbot])
    for ri in range(nbot):
        for rj in range(nbot):
            distance[ri][rj] = np.hypot(ry[ri] - ry[rj], rx[ri] - rx[rj])
    grp = 0
    for i in range(iter):
        # print(i)


        for j in range(nbot):
            bot[j] = motionplanning(ox, oy, reso, rr, motion, rx, ry, pmap, j,nbot, heading,i,  rpos, grp)
            rx[j] = bot[j].ix
            ry[j] = bot[j].iy
            heading[j] = bot[j].heading
            bmap[int(rx[j]), int(ry[j])]=1
            rmap[j][int(rx[j])][int(ry[j])] = 1

        for w in range(nbot):
            rpos[w][0]=rx[w]
            rpos[w][1]=ry[w]
        clustering = DBSCAN(eps=100, min_samples=3).fit(rpos)
        grp = []
        numclus = max(clustering.labels_)
        for w1 in range(numclus + 1):
            # grp[j] = np.where(clustering.labels_ == j)
            grp.append(np.where(clustering.labels_ == w1))
            # clusgoal(grp[w1][0], rpos)


        for ri in range(nbot):
            for rj in range(nbot):
                distance[ri][rj] = np.hypot(bot[ri].ix-bot[rj].ix, bot[ri].iy-bot[rj].iy)
        stepsize_mat.append(step_size[0])
        if show_animation:
            for j in range(nbot):
                plt.plot(bot[j].ix, bot[j].iy, ',')
                # if j == 0:
                #     plt.pause(0.001)
                # else:
                #     plt.plot(bot[j].ix, bot[j].iy, ",g")
                    # plt.pause(0.001)
                # print(j)
    #                plt.ylim((0, 500))
    #                plt.xlim((0, 500))
    #                plt.pause(0.01)

    print("Goal!!")

    return rx, ry


def draw_heatmap(data):
    data = np.array(data).T
    z_min, z_max = data.min(), data.max()
    plt.pcolormesh(data, vmin=z_min, vmax=z_max, cmap='RdBu')
    # plt.show


def main():
    global rmap
    global bmap

    print("potential_field_planning start")

    #    sx = 10.0  # start x position [m]
    #    sy = 20.0  # start y positon [m]
    gx = []  # goal x position [m]
    gy = []  # goal y position [m]
    num_vic = 4
    grid_size = 1  # potential grid size [m] reso
    robot_radius = 3.0  # robot radius [m]
    n = 4  # number of robots
    #    obs = 1; # number of obstacles
    area = 800
    ox = []
    oy = []
    bmap = np.zeros([area, area])
    rmap = [[[0 for k in range(area)] for j in range(area)] for i in range(n)]
    rmap = np.asarray(rmap)
    #
    #    for i in range(obs):
    #        ox.append([random.random()*area, random.random()*area]);
    #        oy.append([random.random()*area, random.random()*area]);

    oy = [[0, 0], [0, 800], [800, 800], [800, 0]]
    ox = [[0, 800], [800, 800], [800, 0], [0, 0]]

    # oy = [[0, int((1 / 3) * area)], [int((1 / 3) * area), int((1 / 3) * area)], [0.5 * area, area],
    #       [int((1 / 3) * area), int((1 / 3) * area)], [0, 0], [0, 800], [800, 800], [800, 0]]
    # ox = [[int(0.8 * area), int(0.8 * area)], [0, int((1 / 3) * area)], [int((1 / 3) * area), int((1 / 3) * area)],
    #       [int(0.9 * area), area], [0, 800], [800, 800], [800, 0], [0, 0]]
    #    ox = [[3, 3],[3, 3]]
    #    oy = [[3, 0],[5, 8]]
    # print(ox)

    for i in range(num_vic):
        gx.append(random.random() * area);
        gy.append(random.random() * area);
    #    ox = [30.0, 10.0, 40.0, 75.0]  # obstacle x position list [m]
    #    oy = [20.0, 30.0, 52.0, 80.0]  # obstacle y position list [m]

    if show_animation:
        plt.grid(True)
        plt.axis("equal")
    xr = []
    yr = []
    for i in range(int(n)):
        xr.append(50+10*random.random())
        # xr.append(600+10*random.random())
        yr.append(200+10*random.random())
        # yr.append(600 + 10 * random.random())
        # xr.append(random.random() * area)
        # yr.append(random.random() * area)
    #        yr.append(random.random()*random.randrange(-1, 2, 2)*area);
    #    breakpoint()
    # path generation
    rx, ry = potential_field_planning(
        xr, yr, gx, gy, ox, oy, grid_size, robot_radius, n, area)
    if show_animation:
        plt.show()

    plt.figure(2)
    draw_heatmap(bmap)


if __name__ == '__main__':
    print(__file__ + " start!!")
    main()
    rmap
    bmap
    distance
    print(__file__ + " Done!!")
    print(np.count_nonzero(bmap == 1))
