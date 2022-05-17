from util import *
import numpy as np

def computedist(x):
    return 2 * np.log((x + np.sqrt(49 * x ** 4 / 144 - x ** 2 / 6 + 1)) /
                  np.sqrt((1 - 7 * x ** 2 / 12) * (1 - 7 * x ** 2 / 12)))

def solvescale(sidelen):
    xmin = 0
    antiprismcr = np.sqrt(4+1/np.sin(np.pi/8)**2)/4
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

class Prism:
    def __init__(self, side_length=1):
        s = side_length
        n = 3
        prismcr = np.sqrt(7 / 12)

        prismcrratio = solvescale(s)
        #projection_dist(prismcr * s) / prismcr

        # 0, 1, 2 going counterclockwise on the top face starting at x-axis
        # 3, 4, 5 going counterclockwise on the bottom face starting at x-axis
        self.pointsnamed = {}
        for k in range(n):
            self.pointsnamed[k] = Point(np.cos(k * 2 * np.pi / n) * np.sqrt(3) / 3 * prismcrratio,
                                        np.sin(k * 2 * np.pi / n) * np.sqrt(3) / 3 * prismcrratio,
                                        prismcrratio / 2)
            self.pointsnamed[k + 3] = Point(np.cos(k * 2 * np.pi / n) * np.sqrt(3) / 3 * prismcrratio,
                                            np.sin(k * 2 * np.pi / n) * np.sqrt(3) / 3 * prismcrratio,
                                            -prismcrratio / 2)
        self.points = self.pointsnamed.values()
        self.edges = set()
        self.edgesnamed = set()
        d = min({a.dist(b) for a in self.points for b in self.points if a != b})

        for i in range(len(self.points)):
            for j in range(len(self.points)):
                a = self.pointsnamed[i]
                b = self.pointsnamed[j]
                if abs(a.dist(b) - d) < delta:
                    self.edges.add(frozenset({a, b}))
                    self.edgesnamed.add(frozenset({i, j}))

        '''# triangle is arrow pointing in positive x direction
        # faces are frozen sets of points
        self.face1 = frozenset([0,1,3,4])  # square face 'above' x-axis birds eye view
        self.face2 = frozenset([1,2,4,5])  # square face 'perpendicular' to x-axis birds eye view
        self.face3 = frozenset([0,2,3,5])  # square face 'below' x-axis birds eye view
        self.topface = frozenset([0,1,2])  # positive z
        self.botface = frozenset([3,4,5])  # negative z'''

        self.face1 = [3, 0, 2]
        self.face2 = [5, 2, 4]
        self.face3 = [4, 1, 3]

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
