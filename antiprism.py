from util import *
import numpy as np


def computedist(x):
    return 2 * np.log((x + np.sqrt((9 + 4 * np.sqrt(2)) * x ** 4 / 32
                            - np.sqrt(2) * x ** 2 / 4 + 1)) /
               np.sqrt((1 - (4 + np.sqrt(2)) * x ** 2 / 8) ** 2))

def solvescale(sidelen):
    xmin = 0
    antiprismcr = 0.822664388008
    xmax = 1/antiprismcr
    while True:
        xmid = (xmax+xmin)/2
        d = computedist(xmid)
        if abs(d - sidelen) < 0.1 ** 15:
            return xmid
        if d > sidelen:
            xmax = xmid
        else:
            xmin = xmid


class Antiprism:
    def __init__(self, side_length=1):
        s = side_length
        n = 4

        antiprismcrratio = solvescale(s)
        # projection_dist(antiprismcr * s) / antiprismcr
        h = ((np.cos(np.pi / n) - np.cos(2 * np.pi / n)) / 2) ** (1 / 2)

        # 0, 2, 4, 6 going counterclockwise on top starting from x-axis
        # 1, 3, 5, 7 going counterclockwise on bottom starting at upper-right
        self.pointsnamed = {}
        for k in range(2 * n):
            self.pointsnamed[k] = Point(np.cos(k * np.pi / n) * np.sqrt(2) / 2 * antiprismcrratio,
                                        np.sin(k * np.pi / n) * np.sqrt(2) / 2 * antiprismcrratio,
                                        (-1) ** k * h * np.sqrt(2) / 2 * antiprismcrratio)
        self.points = self.pointsnamed.values()

        # names edges by endpoint names
        self.edgesnamed = set()
        self.edges = set()
        d = min({a.dist(b) for a in self.points for b in self.points if a != b})

        for i in range(len(self.points)):
            for j in range(len(self.points)):
                a = self.pointsnamed[i]
                b = self.pointsnamed[j]
                if abs(a.dist(b) - d) < delta:
                    self.edges.add(frozenset({a, b}))
                    self.edgesnamed.add(frozenset({i, j}))

        # top square is diamond shape
        # faces are frozen sets of names
        self.topface = frozenset([0, 2, 4, 6])  # positive z
        self.botface = frozenset([1, 3, 5, 7])  # negative z

    # internal update, no need to call outside
    def update(self, a=None, b=None):
        self.points = self.pointsnamed.values()
        self.edges = set()
        for i in range(len(self.points)):
            for j in range(len(self.points)):
                if frozenset([i, j]) in self.edgesnamed:
                    self.edges.add(frozenset({self.pointsnamed[i],
                                              self.pointsnamed[j]}))
        # self.edges.add(a)
        # self.edges.add(b)

    # modifies antiprism so that a to A, b to B, and c to C, preserves orientation
    # a, b, c are label names, A, B, C external points (and so are expressed as coordinates)
    def moveface(self, a, b, c, A, B, C):
        x, y, z = A.copy(), B.copy(), C.copy()
        center, radius = perpendicularbisector(self.pointsnamed[a], Point(0, 0, 0))
        for p in range(len(self.points)):
            self.pointsnamed[p] = self.pointsnamed[p].inverse(center, radius)
        center2, radius2 = perpendicularbisector(x, Point(0, 0, 0))
        x = x.inverse(center2, radius2)
        y = y.inverse(center2, radius2)
        z = z.inverse(center2, radius2)

        v1 = y.minus(x)
        v2 = z.minus(x)
        nx = v1.cross(v2).scale()
        v1 = self.pointsnamed[b].minus(self.pointsnamed[a])
        v2 = self.pointsnamed[c].minus(self.pointsnamed[a])
        na = v1.cross(v2).scale()
        axis = na.cross(nx)
        theta = np.arccos(nx.dot(na))
        rm = rotation_matrix([axis.x, axis.y, axis.z], theta)
        for p in range(len(self.points)):
            rp = np.dot(rm, [self.pointsnamed[p].x, self.pointsnamed[p].y,
                             self.pointsnamed[p].z])
            self.pointsnamed[p] = Point(rp[0], rp[1], rp[2])
        v1 = self.pointsnamed[b].scale()
        v2 = y.scale()
        axis = v1.cross(v2)
        theta = np.arccos(v1.dot(v2))
        rm = rotation_matrix([axis.x, axis.y, axis.z], theta)
        for p in range(len(self.points)):
            rp = np.dot(rm, [self.pointsnamed[p].x, self.pointsnamed[p].y,
                             self.pointsnamed[p].z])
            self.pointsnamed[p] = Point(rp[0], rp[1], rp[2])

        for p in range(len(self.points)):
            self.pointsnamed[p] = self.pointsnamed[p].inverse(center2, radius2)

        self.update()


'''
#pt1, pt2 = spherecircleintersect(center, radius, self.pointsnamed[b], B, B.inverse())
#midpt = pt1 if pt1.dist() < pt2.dist() else pt2
#print(midpt)
#center, radius = threepointcircle(midpt, A, midpt.inverse())
#self.update(frozenset([midpt, A]), frozenset([B, self.pointsnamed[b]]))'''
