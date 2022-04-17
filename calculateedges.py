# i did some comments, weeeee! still hard to read tho, is very un-good code

import numpy as np
import os

# cutoff for when doubles need to check for 'equal', otherwise finicky because rounding and stuff
delta = 0.0000001


# 3d point class with some functions
class Point:
    # x, y, z coordinates
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    # finds distance between self and other, if no other, then finds norm aka distance from origin
    def dist(self, other=None):
        if not other:
            other = Point(0, 0, 0)
        return ((self.x - other.x) ** 2 + (self.y - other.y) ** 2 + (self.z - other.z) ** 2) ** (1 / 2)

    # if circle is infinite radius, is the same as regular reflection
    # indicate this special case by r = -1 and c = (i,j,k) is three points on the plane
    # otherwise, takes circle inversion, where circle is inputted as center c and radius r
    # if no arguments, then takes inverse with unit circle centered at origin
    def inverse(self, c=None, r=1):
        if c is None:
            c = Point(0, 0, 0)
        # plane inversion
        if r == -1:
            (i, j, k) = c
            # ax + by + cz + d = 0 is equation of plane
            v1 = j.minus(i)
            v2 = k.minus(i)
            n = v1.cross(v2)
            a = n.x
            b = n.y
            c = n.z
            d = -n.x * i.x - n.y * i.y - n.z * i.z
            # reflects over plane we found
            k = (-a * self.x - b * self.y - c * self.z - d) / (a * a + b * b + c * c)
            x2 = a * k + self.x
            y2 = b * k + self.y
            z2 = c * k + self.z
            x3 = 2 * x2 - self.x
            y3 = 2 * y2 - self.y
            z3 = 2 * z2 - self.z
            return Point(x3, y3, z3)
        # regular inversion over circle
        d = self.dist(c)
        if d == 0:
            print("origin has no inverse")
            return None
        return Point(c.x + (r ** 2) * (self.x - c.x) / (d ** 2),
                     c.y + (r ** 2) * (self.y - c.y) / (d ** 2),
                     c.z + (r ** 2) * (self.z - c.z) / (d ** 2))

    # returns midpoint between self and other
    def midpoint(self, other):
        return Point((self.x + other.x) / 2, (self.y + other.y) / 2, (self.z + other.z) / 2)

    # returns self - other
    def minus(self, other):
        return Point(self.x-other.x, self.y-other.y, self.z-other.z)

    # returns self + other
    def plus(self, other):
        return Point(self.x+other.x, self.y+other.y, self.z+other.z)

    # returns dot product self * other
    def dot(self, other):
        return self.x*other.x + self.y*other.y + self.z*other.z

    # returns cross product self x other
    def cross(self, other):
        return Point(self.y*other.z-self.z*other.y, -(self.x*other.z-self.z*other.x),
                     self.x*other.y-self.y*other.x)

    # scales self by scalar alpha
    def scale(self, alpha):
        return Point(self.x * alpha, self.y * alpha, self.z * alpha)

    # equality comparison
    def __eq__(self, other):
        return abs(self.x - other.x) < delta and abs(self.y - other.y) < delta and abs(self.z - other.z) < delta

    # something something hash
    def __hash__(self):
        return round(1 / delta * ((self.x + self.y) * (self.x + self.z + 1) * (self.y + self.z + 2) / 2) + self.y)

    # checks for collinearity but i dont think i ever use it since i check if det is zero instead for methods below
    def isCollinear(self, o1, o2):
        dist1 = self.dist(o1)
        dist2 = self.dist(o2)
        dist3 = o1.dist(o2)
        maxdist = max(dist1, dist2, dist3)
        sumdist = dist1+dist2+dist3
        # basically triangle inequality but make sure its the equality
        if abs(sumdist - 2 * maxdist) < delta:
            return True
        return False

# takes euclidean distance and finds ball distance but i think its only valid for diameters probably
# i found some stuff on wikipedia and did some un-good math for some special case for something
def projection_dist(d):
    return (np.exp(d) - 1)/(np.exp(d) + 1)

# very definitely stolen am not this clevery
# if points are all on same plane, then sphere is undefined, use plane reflection instead
# return indicates by returning three points on plane and -1 for radius
def fourpointsphere(a, b, c, d):
    A = a.minus(d)
    B = b.minus(d)
    C = c.minus(d)
    # something something linear algebra
    det = A.x * (B.y * C.z - B.z * C.y) - A.y * (B.x * C.z - B.z * C.x) + A.z * (B.x * C.y - B.y * C.x)
    if abs(det) < delta:
        return (a,b,c), -1
    j = A.cross(B).scale(C.dot(C))
    k = C.cross(A).scale(B.dot(B))
    l = B.cross(C).scale(A.dot(A))
    m = j.plus(k).plus(l).scale(0.5 / det)
    radius = m.dist()
    center = d.plus(m)
    return center, radius

# also very definitely stolen, stackexchange people are very smart
# finds circle through three points, if they are collinear, returns -1 for both center, radius (shouldn't happen)
def threepointcircle(a,b,c):
    t = b.minus(a)
    u = c.minus(a)
    v = c.minus(b)
    w = t.cross(u)
    wsl = w.dist() ** 2
    if abs(wsl) < delta:
        return None, -1
    iwsl2 = 0.5 / wsl
    center = a.plus((u.scale(t.dot(t)*u.dot(v)).minus(t.scale(u.dot(u)*t.dot(v)))).scale(iwsl2))  # i am sorry
    radius = np.sqrt(t.dot(t) * u.dot(u) * v.dot(v) * iwsl2 * 0.5)
    return center, radius

# given circle with center c and radius r, and line passing through points p and q, find two intersections (hopefully not tangent)
# stolen from paul bourke's website
def linecircleintersections(p, q, c=None, r=1):
    if c is None:
        c = Point(0.0, 0.0, 0.0)
    x1, y1, z1 = p.x, p.y, p.z  # SORRY!
    x2, y2, z2 = q.x, q.y, q.z
    x3, y3, z3 = c.x, c.y, c.z
    a = (x2-x1) ** 2 + (y2-y1) ** 2 + (z2-z1) ** 2
    b = 2 * ((x2-x1)*(x1-x3) + (y2-y1)*(y1-y3) + (z2-z1)*(z1-z3))
    c = x3**2 + y3**2 + z3**2 + x1**2 + y1**2 + z1**2 - 2 * (x3*x1 + y3*y1 + z3*z1) - r**2
    u1 = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
    u2 = (-b - np.sqrt(b**2 - 4*a*c))/(2*a)
    i1 = p.plus(q.minus(p).scale(u1))
    i2 = p.plus(q.minus(p).scale(u2))
    return i1, i2

# schlafli symbol, p indicates p-gon, q is number of faces per vertex, r is number of polyhedra per edge
p = 4
q = 3
r = 5

# some magical formulas i pulled off some paper i found online
# finds circumradius and inradius of center polyhedron
hpq = np.arccos(np.sqrt(np.cos(np.pi/p) ** 2 + np.cos(np.pi/q) ** 2))
hqr = np.arccos(np.sqrt(np.cos(np.pi/q) ** 2 + np.cos(np.pi/r) ** 2))
ir = np.arccosh(np.sin(np.pi/p) * np.cos(np.pi/r) / np.sin(hpq))      # inradius
cr = np.arccosh(np.cos(np.pi/p) * np.cos(np.pi/q) * np.cos(np.pi/r) / # circumradius
               (np.sin(hpq) * np.sin(hqr)))
ir = projection_dist(ir)
cr = projection_dist(cr)

# ok ok so this is only good for cubes rn, will have to find the initial array for the center
# shape for other platonic solids, seems tedious, leaving as TODO!
d = cr / np.sqrt(3)
corners = set()
for i in [-d, d]:
    for j in [-d, d]:
        for k in [-d, d]:
            corners.add(Point(i, j, k))

corners = frozenset(corners)
edges = set()
for c1 in corners:
    for c2 in corners:
        if abs(c1.dist(c2) - 2*d) < delta:
            edges.add(frozenset([c1, c2]))

# this is lame vvv oh well
facesnotfrozen = [set() for i in range(6)]
for c in corners:
    if c.x > 0: facesnotfrozen[0].add(c)
    else: facesnotfrozen[1].add(c)
    if c.y > 0: facesnotfrozen[2].add(c)
    else: facesnotfrozen[3].add(c)
    if c.z > 0: facesnotfrozen[4].add(c)
    else: facesnotfrozen[5].add(c)

faces = set()
for i in range(6):
    s = facesnotfrozen[i]
    faces.add(frozenset(s))
faces = frozenset(faces)

# maps each polyhedron with its faces, index is frozenset of corners (kinda stupid but it works)
cubes = {corners: faces}
newcubes = dict()
newedges = set()

# loops through cubes, reflects cube over each face and updates the 'cubes' map
# then goes through all edges, checks if edge is in the current cube, it it is, reflect!
# very inefficient rn, but small enough dataset that runtime is negligible
for cube in cubes:
    for face in cubes[cube]:
        c1, c2, c3, c4 = face
        c, r = fourpointsphere(c1, c2, c3, c1.inverse())
        cori = frozenset([corner.inverse(c, r) for corner in cube])
        newfaces = set()
        for faceother in cubes[cube]:
            o1, o2, o3, o4 = faceother
            newfaces.add(frozenset([o1.inverse(c,r), o2.inverse(c,r), o3.inverse(c,r), o4.inverse(c,r)]))
        newfaces = frozenset(newfaces)
        newcubes[cori] = newfaces
        for edge in edges:
            e1, e2 = edge
            if e1 in cube and e2 in cube:
                newedges.add(frozenset([e1.inverse(c, r), e2.inverse(c, r)]))

# should probably refactor, since its basically the previous loop, but switching roles of newcubes and cubes
# wont let me add to the set while looping through the set so thats what i ended up doing :(
# since its just 3 iterations anyways, was too lazy to make recursion method
for cube in newcubes:
    for face in newcubes[cube]:
        c1, c2, c3, c4 = face
        c, r = fourpointsphere(c1, c2, c3, c1.inverse())
        cori = frozenset([corner.inverse(c, r) for corner in cube])
        newfaces = set()
        for faceother in newcubes[cube]:
            o1, o2, o3, o4 = faceother
            newfaces.add(frozenset([o1.inverse(c,r), o2.inverse(c,r), o3.inverse(c,r), o4.inverse(c,r)]))
        newfaces = frozenset(newfaces)
        cubes[cori] = newfaces
        for edge in newedges:
            e1, e2 = edge
            if e1 in cube and e2 in cube:
                edges.add(frozenset([e1.inverse(c, r), e2.inverse(c, r)]))

cutoff = 0.95
# third iteration! ctrl-c, ctrl-v of first iteration
# if you are tryna debug, should comment this iteration out, and change 'newedges' to 'edges' in the printout loop below
# takes forever to draw 1400 edges!!!
for cube in cubes:
    for face in cubes[cube]:
        c1, c2, c3, c4 = face
        c, r = fourpointsphere(c1, c2, c3, c1.inverse())
        cori = frozenset([corner.inverse(c, r) for corner in cube])
        maxdist = max([cc.dist() for cc in cori])
        if maxdist > cutoff: continue
        newfaces = set()
        for faceother in cubes[cube]:
            o1, o2, o3, o4 = faceother
            newfaces.add(frozenset([o1.inverse(c,r), o2.inverse(c,r), o3.inverse(c,r), o4.inverse(c,r)]))
        newfaces = frozenset(newfaces)
        newcubes[cori] = newfaces
        for edge in edges:
            e1, e2 = edge
            if e1 in cube and e2 in cube:
                newedges.add(frozenset([e1.inverse(c, r), e2.inverse(c, r)]))
print(len(newedges))

for cube in newcubes:
    for face in newcubes[cube]:
        c1, c2, c3, c4 = face
        c, r = fourpointsphere(c1, c2, c3, c1.inverse())
        cori = frozenset([corner.inverse(c, r) for corner in cube])
        maxdist = max([cc.dist() for cc in cori])
        if maxdist > cutoff: continue
        newfaces = set()
        for faceother in newcubes[cube]:
            o1, o2, o3, o4 = faceother
            newfaces.add(frozenset([o1.inverse(c,r), o2.inverse(c,r), o3.inverse(c,r), o4.inverse(c,r)]))
        newfaces = frozenset(newfaces)
        cubes[cori] = newfaces
        for edge in newedges:
            e1, e2 = edge
            if e1 in cube and e2 in cube:
                edges.add(frozenset([e1.inverse(c, r), e2.inverse(c, r)]))
print(len(edges))

for cube in cubes:
    for face in cubes[cube]:
        c1, c2, c3, c4 = face
        c, r = fourpointsphere(c1, c2, c3, c1.inverse())
        cori = frozenset([corner.inverse(c, r) for corner in cube])
        maxdist = max([cc.dist() for cc in cori])
        if maxdist > cutoff: continue
        newfaces = set()
        for faceother in cubes[cube]:
            o1, o2, o3, o4 = faceother
            newfaces.add(frozenset([o1.inverse(c,r), o2.inverse(c,r), o3.inverse(c,r), o4.inverse(c,r)]))
        newfaces = frozenset(newfaces)
        newcubes[cori] = newfaces
        for edge in edges:
            e1, e2 = edge
            if e1 in cube and e2 in cube:
                newedges.add(frozenset([e1.inverse(c, r), e2.inverse(c, r)]))
print(len(newedges))

# given hyperbolic center and radius, finds euclidean center and radius
def findeuclideanspherewithhyperboliccenterandradius(p, r):
    a, b = linecircleintersections(p, Point(0,0,0))
    ratio = r * p.dist(b) / p.dist(a)
    p1 = a.plus(b.minus(a).scale(1 / (ratio + 1)))
    ratio = r * p.dist(a) / p.dist(b)
    p2 = a.plus(b.minus(a).scale(ratio / (ratio + 1)))
    ec1 = p1.midpoint(p2)
    r1 = ec1.dist(p1)
    return ec1, r1

points = set()
piperadius = 0.20
piperad = np.exp(piperadius)

# printout into csv of coordinates that i feed to solidworks
# the 999999999 looks stupid, but basically indicates to solidworks that i want a line, not an arc
with open('edges.txt', 'w') as f:
    for edge in newedges:
        c1, c2 = edge
        c3 = c1.inverse()
        c, r = threepointcircle(c1, c2, c3)
        mid = c1.midpoint(c2)
        C1, R1 = findeuclideanspherewithhyperboliccenterandradius(c1, piperad)
        C2, R2 = findeuclideanspherewithhyperboliccenterandradius(c2, piperad)
        points.add((C1, R1))
        points.add((C2, R2))
        if r == -1:
            C3 = Point(999999999, 999999999, 999999999)
            _, R3 = findeuclideanspherewithhyperboliccenterandradius(mid, piperad)
        else:
            # finds midpoint of arc, so that threepointarc in solidworks has a third point to work with
            diff = mid.minus(c)
            sdiff = diff.scale(r/diff.dist())
            c3 = sdiff.plus(c)
            C3, R3 = findeuclideanspherewithhyperboliccenterandradius(c3, piperad)
            # writes to csv file
        c, r = threepointcircle(C1, C2, C3)
        C4, _ = linecircleintersections(C1, c, C1, R1)
        C5, _ = linecircleintersections(C2, c, C2, R2)
        C6, _ = linecircleintersections(C3, c, C3, R3)
        f.write(str(C1.x) + ", " + str(C1.y) + ", " + str(C1.z) + ", " +
                str(C2.x) + ", " + str(C2.y) + ", " + str(C2.z) + ", " +
                str(C3.x) + ", " + str(C3.y) + ", " + str(C3.z) + ", " +
                str(C4.x) + ", " + str(C4.y) + ", " + str(C4.z) + ", " +
                str(C5.x) + ", " + str(C5.y) + ", " + str(C5.z) + ", " +
                str(C6.x) + ", " + str(C6.y) + ", " + str(C6.z) + ", ")

with open('edges.txt', 'rb+') as f:
    f.seek(-2, os.SEEK_END)
    f.truncate()

biggerpts = dict()
for p, r in points:
    if p not in biggerpts.keys():
        biggerpts[p] = r
    else:
        biggerpts[p] = max(biggerpts[p], r)

with open('points.txt', 'w') as f:
    for p in biggerpts.keys():
        r = biggerpts[p]
        f.write(str(p.x) + ", " + str(p.y) + ", " + str(p.z) + ", " + str(r) + ", ")

print(len(biggerpts))

with open('points.txt', 'rb+') as f:
    f.seek(-2, os.SEEK_END)
    f.truncate()