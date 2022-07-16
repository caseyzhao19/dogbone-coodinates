from util import *
import numpy as np


def computedist(x):
    return 2 * np.log((x + np.sqrt((9 + 4 * np.sqrt(2)) * x ** 4 / 32
                                   - np.sqrt(2) * x ** 2 / 4 + 1)) /
                      np.sqrt((1 - (4 + np.sqrt(2)) * x ** 2 / 8) ** 2))


def solvescale(sidelen):
    xmin = 0
    antiprismcr = 0.822664388008
    xmax = 1 / antiprismcr
    while True:
        xmid = (xmax + xmin) / 2
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
    def moveface(self, a, b, c, A, B, C, f=None, cam_angle=0):
        if f:
            triangles = []
            for j in range(8):
                triangles.append([self.pointsnamed[j], self.pointsnamed[(j + 1) % 8],
                                  self.pointsnamed[(j + 2) % 8]])
            triangles.append([self.pointsnamed[0], self.pointsnamed[2], self.pointsnamed[6]])

            triangles.append([self.pointsnamed[2], self.pointsnamed[4], self.pointsnamed[6]])
            triangles.append([self.pointsnamed[1], self.pointsnamed[5], self.pointsnamed[7]])
            triangles.append([self.pointsnamed[1], self.pointsnamed[3], self.pointsnamed[5]])
            for a1, b1, c1 in triangles:
                d1 = a1.inverse()
                center, radius = fourpointsphere(a1, b1, c1, d1, False)
                if radius != -1:
                    f.write("object {\n")
                    f.write(
                        "\t sphere { <" + str(center.x) + "," + str(center.y) + "," + str(center.z) + ">," + str(
                            radius) + "\n")
                    f.write("\t\t texture{ pigment{ color rgb <1, 0.24, 0.1> filter 0.5}}\n")
                    f.write("\t\t}\n")
                    f.write("\t clipped_by {mesh {\n")
                    O = Point(0, 0, 0)
                    for p, q, r in [(a1, b1, c1), (O, a1, b1), (O, a1, c1), (O, b1, c1)]:
                        f.write(
                            "\t\t triangle { <" + str(p.x) + "," + str(p.y) + "," + str(p.z) + ">, <" + str(q.x) + "," +
                            str(q.y) + "," + str(q.z) + ">, <" + str(r.x) + "," + str(r.y) + "," + str(r.z) + "> }\n")
                    f.write("\t\t inside_vector <0, 0, 1>}}\n")
                    f.write("rotate <" + str(cam_angle) + ", 0, 0>}\n")
                else:
                    f.write("object {\n")
                    f.write("\t mesh {\n")
                    f.write("\t\t triangle { <" + str(center[0].x) + "," + str(center[0].y) + "," + str(
                        center[0].z) + ">, <"
                            + str(center[1].x) + "," + str(center[1].y) + "," + str(center[1].z) + ">, <"
                            + str(center[2].x) + "," + str(center[2].y) + "," + str(center[2].z) +
                            ">}} texture{ pigment{ color rgb <1, 0.24, 0.1> filter 0.5}} rotate <" + str(
                        cam_angle) + ", 0, 0>}\n")
            print("mod:", c)
            for a1, b1 in self.edgesnamed:
                a1 = self.pointsnamed[a1]
                b1 = self.pointsnamed[b1]
                c1 = a1.inverse()
                center, radius = threepointcircle(a1, b1, c1)
                print(radius)
                if radius != -1:
                    n = a1.minus(c1).cross(b1.minus(c1)).scale()
                    default = Point(0, 1, 0)
                    axis = n.cross(default)
                    theta = np.arccos(n.dot(default))
                    rm = rotation_matrix([axis.x, axis.y, axis.z], theta)
                    f.write("intersection {")
                    f.write("object {")
                    f.write(
                        "torus {" + str(radius) + ", 0.005} texture{pigment{ color rgb <0.4, 0.05, 0.05>}} matrix <")
                    for row in [0, 1, 2]:
                        for col in [0, 1, 2]:
                            f.write(str(rm[row][col]))
                            f.write(", ")
                    f.write(str(center.x) + ", " + str(center.y) + ", " + str(center.z))
                    center, radius = fourpointsphere(self.pointsnamed[0], self.pointsnamed[1],
                                                     self.pointsnamed[2], self.pointsnamed[3])
                    f.write(">} sphere {<" + str(center.x) + ", " + str(center.y) + ", " + str(center.z) + "> "
                            + str(radius + 0) + " pigment{ color rgb <0.4, 0.05, 0.05>}} rotate <" + str(
                        cam_angle) + ", 0, 0>}\n")
                else:
                    f.write("#declare myspline = spline{linear_spline 0, <"
                            + str(a1.x) + ", " + str(a1.y) + ", " + str(a1.z) +
                            "> 1, <" + str(b1.x) + ", " + str(b1.y) + ", " + str(b1.z) +
                            ">}\n")
                    f.write("#declare ctr=0; #while (ctr < 1) sphere { myspline(ctr), "
                            "0.005 pigment { rgb <0.4, 0.05, 0.05> } rotate <" + str(
                        cam_angle) + ", 0, 0>} #declare ctr = ctr + 0.01; #end\n")
        x, y, z = A.copy(), B.copy(), C.copy()
        center, radius = perpendicularbisector(self.pointsnamed[a], Point(0, 0, 0))
        if f:
            # reflect center
            f.write("intersection {\n")
            f.write("\t sphere { <" + str(center.x) + ", " + str(center.y) + ", " + str(center.z) + ">," + str(
                radius) + "\n")
            f.write("\t\t texture{ pigment{ color rgb <1, 0.85, 0.8> filter 0.999}} ")
            f.write("\t\t} sphere{<0,0,0>, 1 texture{ pigment{ color rgb <1, 1, 1> filter 0.999}}} rotate <" + str(
                cam_angle) + ", 0, 0>}\n")
        for p in range(len(self.points)):
            self.pointsnamed[p] = self.pointsnamed[p].inverse(center, radius)
        if f:
            triangles = []
            for j in range(8):
                triangles.append([self.pointsnamed[j], self.pointsnamed[(j + 1) % 8],
                                  self.pointsnamed[(j + 2) % 8]])
            triangles.append([self.pointsnamed[0], self.pointsnamed[2], self.pointsnamed[6]])

            triangles.append([self.pointsnamed[2], self.pointsnamed[4], self.pointsnamed[6]])
            triangles.append([self.pointsnamed[1], self.pointsnamed[5], self.pointsnamed[7]])
            triangles.append([self.pointsnamed[1], self.pointsnamed[3], self.pointsnamed[5]])
            for a1, b1, c1 in triangles:
                d1 = a1.inverse()
                center, radius = fourpointsphere(a1, b1, c1, d1, False)
                if radius != -1:
                    f.write("object {\n")
                    f.write(
                        "\t sphere { <" + str(center.x) + "," + str(center.y) + "," + str(center.z) + ">," + str(
                            radius) + "\n")
                    f.write("\t\t texture{ pigment{ color rgb <1, 0.24, 0.1> filter 0.5}}\n")
                    f.write("\t\t}\n")
                    f.write("\t clipped_by {mesh {\n")
                    O = Point(0, 0, 0)
                    for p, q, r in [(a1, b1, c1), (O, a1, b1), (O, a1, c1), (O, b1, c1)]:
                        f.write(
                            "\t\t triangle { <" + str(p.x) + "," + str(p.y) + "," + str(p.z) + ">, <" + str(q.x) + "," +
                            str(q.y) + "," + str(q.z) + ">, <" + str(r.x) + "," + str(r.y) + "," + str(r.z) + "> }\n")
                    f.write("\t\t inside_vector <0, 0, 1>}}\n")
                    f.write("rotate <" + str(cam_angle) + ", 0, 0>}\n")
                else:
                    f.write("object {\n")
                    f.write("\t mesh {\n")
                    f.write("\t\t triangle { <" + str(center[0].x) + "," + str(center[0].y) + "," + str(
                        center[0].z) + ">, <"
                            + str(center[1].x) + "," + str(center[1].y) + "," + str(center[1].z) + ">, <"
                            + str(center[2].x) + "," + str(center[2].y) + "," + str(center[2].z) +
                            ">}} texture{ pigment{ color rgb <1, 0.24, 0.1> filter 0.5}} rotate <" + str(
                        cam_angle) + ", 0, 0>}\n")
            print("mod:", c)
            for a1, b1 in self.edgesnamed:
                a1 = self.pointsnamed[a1]
                b1 = self.pointsnamed[b1]
                c1 = a1.inverse()
                center, radius = threepointcircle(a1, b1, c1)
                print(radius)
                if radius != -1:
                    n = a1.minus(c1).cross(b1.minus(c1)).scale()
                    default = Point(0, 1, 0)
                    axis = n.cross(default)
                    theta = np.arccos(n.dot(default))
                    rm = rotation_matrix([axis.x, axis.y, axis.z], theta)
                    f.write("intersection {")
                    f.write("object {")
                    f.write(
                        "torus {" + str(radius) + ", 0.005} texture{pigment{ color rgb <0.4, 0.05, 0.05>}} matrix <")
                    for row in [0, 1, 2]:
                        for col in [0, 1, 2]:
                            f.write(str(rm[row][col]))
                            f.write(", ")
                    f.write(str(center.x) + ", " + str(center.y) + ", " + str(center.z))
                    center, radius = fourpointsphere(self.pointsnamed[0], self.pointsnamed[1],
                                                     self.pointsnamed[2], self.pointsnamed[3])
                    f.write(">} sphere {<" + str(center.x) + ", " + str(center.y) + ", " + str(center.z) + "> "
                            + str(radius + 0) + " pigment{ color rgb <0.4, 0.05, 0.05>}} rotate <" + str(
                        cam_angle) + ", 0, 0>}\n")
                else:
                    f.write("#declare myspline = spline{linear_spline 0, <"
                            + str(a1.x) + ", " + str(a1.y) + ", " + str(a1.z) +
                            "> 1, <" + str(b1.x) + ", " + str(b1.y) + ", " + str(b1.z) +
                            ">}\n")
                    f.write("#declare ctr=0; #while (ctr < 1) sphere { myspline(ctr), "
                            "0.005 pigment { rgb <0.4, 0.05, 0.05> } rotate <" + str(
                        cam_angle) + ", 0, 0>} #declare ctr = ctr + 0.01; #end\n")
        center2, radius2 = perpendicularbisector(x, Point(0, 0, 0))
        if f:
            # reflect dest
            f.write("intersection {\n")
            f.write("\t sphere { <" + str(center2.x) + ", " + str(center2.y) + ", " + str(center2.z) + ">," + str(
                radius2) + "\n")
            f.write("\t\t texture{ pigment{ color rgb <0.8, 0.92, 1> filter 0.999}} ")
            f.write("\t\t} sphere{<0,0,0>, 1 texture{ pigment{ color rgb <1, 1, 1> filter 0.999}}} rotate <" + str(
                cam_angle) + ", 0, 0>}\n")
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

        if f:
            triangles = []
            desttriangles = []
            for j in range(8):
                triangles.append([self.pointsnamed[j], self.pointsnamed[(j + 1) % 8],
                                  self.pointsnamed[(j + 2) % 8]])
            desttriangles.append([self.pointsnamed[0], self.pointsnamed[2], self.pointsnamed[6]])

            desttriangles.append([self.pointsnamed[2], self.pointsnamed[4], self.pointsnamed[6]])
            triangles.append([self.pointsnamed[1], self.pointsnamed[5], self.pointsnamed[7]])
            triangles.append([self.pointsnamed[1], self.pointsnamed[3], self.pointsnamed[5]])
            for a1, b1, c1 in desttriangles:
                d1 = a1.inverse()
                center, radius = fourpointsphere(a1, b1, c1, d1, False)
                if radius != -1:
                    f.write("object {\n")
                    f.write(
                        "\t sphere { <" + str(center.x) + "," + str(center.y) + "," + str(center.z) + ">," + str(
                            radius) + "\n")
                    f.write("\t\t texture{ pigment{ color rgb <1, 0.24, 0.1> filter 0.5}}\n")
                    f.write("\t\t}\n")
                    f.write("\t clipped_by {mesh {\n")
                    O = Point(0, 0, 0)
                    for p, q, r in [(a1, b1, c1), (O, a1, b1), (O, a1, c1), (O, b1, c1)]:
                        f.write(
                            "\t\t triangle { <" + str(p.x) + "," + str(p.y) + "," + str(p.z) + ">, <" + str(q.x) + "," +
                            str(q.y) + "," + str(q.z) + ">, <" + str(r.x) + "," + str(r.y) + "," + str(r.z) + "> }\n")
                    f.write("\t\t inside_vector <0, 0, 1>}}\n")
                    f.write("rotate <" + str(cam_angle) + ", 0, 0>}\n")
                else:
                    f.write("object {\n")
                    f.write("\t mesh {\n")
                    f.write("\t\t triangle { <" + str(center[0].x) + "," + str(center[0].y) + "," + str(
                        center[0].z) + ">, <"
                            + str(center[1].x) + "," + str(center[1].y) + "," + str(center[1].z) + ">, <"
                            + str(center[2].x) + "," + str(center[2].y) + "," + str(center[2].z) +
                            ">}} texture{ pigment{ color rgb <0.1, 0.9, 1> filter 0.5}} rotate <" + str(
                        cam_angle) + ", 0, 0>}\n")
            for a1, b1, c1 in triangles:
                d1 = a1.inverse()
                center, radius = fourpointsphere(a1, b1, c1, d1, False)
                if radius != -1:
                    f.write("object {\n")
                    f.write(
                        "\t sphere { <" + str(center.x) + "," + str(center.y) + "," + str(center.z) + ">," + str(
                            radius) + "\n")
                    f.write("\t\t texture{ pigment{ color rgb <1, 0.24, 0.1> filter 0.5}}\n")
                    f.write("\t\t}\n")
                    f.write("\t clipped_by {mesh {\n")
                    O = Point(0, 0, 0)
                    for p, q, r in [(a1, b1, c1), (O, a1, b1), (O, a1, c1), (O, b1, c1)]:
                        f.write(
                            "\t\t triangle { <" + str(p.x) + "," + str(p.y) + "," + str(p.z) + ">, <" + str(q.x) + "," +
                            str(q.y) + "," + str(q.z) + ">, <" + str(r.x) + "," + str(r.y) + "," + str(r.z) + "> }\n")
                    f.write("\t\t inside_vector <0, 0, 1>}}\n")
                    f.write("rotate <" + str(cam_angle) + ", 0, 0>}\n")
                else:
                    f.write("object {\n")
                    f.write("\t mesh {\n")
                    f.write("\t\t triangle { <" + str(center[0].x) + "," + str(center[0].y) + "," + str(
                        center[0].z) + ">, <"
                            + str(center[1].x) + "," + str(center[1].y) + "," + str(center[1].z) + ">, <"
                            + str(center[2].x) + "," + str(center[2].y) + "," + str(center[2].z) +
                            ">}} texture{ pigment{ color rgb <1, 0.24, 0.1> filter 0.5}} rotate <" + str(
                        cam_angle) + ", 0, 0>}\n")
            print("mod:", c)
            for a1, b1 in self.edgesnamed:
                a1 = self.pointsnamed[a1]
                b1 = self.pointsnamed[b1]
                c1 = a1.inverse()
                center, radius = threepointcircle(a1, b1, c1)
                print(radius)
                if radius != -1:
                    n = a1.minus(c1).cross(b1.minus(c1)).scale()
                    default = Point(0, 1, 0)
                    axis = n.cross(default)
                    theta = np.arccos(n.dot(default))
                    rm = rotation_matrix([axis.x, axis.y, axis.z], theta)
                    f.write("intersection {")
                    f.write("object {")
                    f.write(
                        "torus {" + str(radius) + ", 0.005} texture{pigment{ color rgb <0.4, 0.05, 0.05>}} matrix <")
                    for row in [0, 1, 2]:
                        for col in [0, 1, 2]:
                            f.write(str(rm[row][col]))
                            f.write(", ")
                    f.write(str(center.x) + ", " + str(center.y) + ", " + str(center.z))
                    center, radius = fourpointsphere(self.pointsnamed[0], self.pointsnamed[1],
                                                     self.pointsnamed[2], self.pointsnamed[3])
                    f.write(">} sphere {<" + str(center.x) + ", " + str(center.y) + ", " + str(center.z) + "> "
                            + str(radius + 0) + " pigment{ color rgb <0.4, 0.05, 0.05>}} rotate <" + str(
                        cam_angle) + ", 0, 0>}\n")
                else:
                    f.write("#declare myspline = spline{linear_spline 0, <"
                            + str(a1.x) + ", " + str(a1.y) + ", " + str(a1.z) +
                            "> 1, <" + str(b1.x) + ", " + str(b1.y) + ", " + str(b1.z) +
                            ">}\n")
                    f.write("#declare ctr=0; #while (ctr < 1) sphere { myspline(ctr), "
                            "0.005 pigment { rgb <0.4, 0.05, 0.05> } rotate <" + str(
                        cam_angle) + ", 0, 0>} #declare ctr = ctr + 0.01; #end\n")

        for p in range(len(self.points)):
            self.pointsnamed[p] = self.pointsnamed[p].inverse(center2, radius2)

        self.update()


'''
#pt1, pt2 = spherecircleintersect(center, radius, self.pointsnamed[b], B, B.inverse())
#midpt = pt1 if pt1.dist() < pt2.dist() else pt2
#print(midpt)
#center, radius = threepointcircle(midpt, A, midpt.inverse())
#self.update(frozenset([midpt, A]), frozenset([B, self.pointsnamed[b]]))'''
