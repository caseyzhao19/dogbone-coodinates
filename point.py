# i did some comments, weeeee! still hard to read tho, is very un-good code

import numpy as np
import os

# cutoff for when doubles need to check for 'equal', otherwise finicky because rounding and stuff
delta = 0.0000000001


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

    # scales self by scalar alpha, no parameter returns normalized
    def scale(self, alpha=None):
        if alpha is None:
            alpha = 1/self.dist()
        return Point(self.x * alpha, self.y * alpha, self.z * alpha)

    # equality comparison
    def __eq__(self, other):
        return self.dist(other) < delta

    # something something hash
    def __hash__(self):
        return round(1 / delta * ((self.x + self.y) * (self.x + self.z + 1) * (self.y + self.z + 2) / 2) + self.y)

    def slope(self, other):
        if abs(other.x - self.x) < delta:
            return None
        return (other.y - self.y) / (other.x - self.x)

    # checks for collinearity but i dont think i ever use it since i check if det is zero instead for methods below
    def isCollinear(self, o1, o2):
        dist1 = self.dist(o1)
        dist2 = self.dist(o2)
        dist3 = o1.dist(o2)
        maxdist = max(dist1, dist2, dist3)
        sumdist = dist1+dist2+dist3
        # basically triangle inequality but make sure its the equality
        if abs(sumdist - 2 * maxdist) < delta**5:
            return True
        return False

    def __str__(self):
        return str(self.x) + " " + str(self.y) + " " + str(self.z)

    def copy(self):
        return Point(self.x, self.y, self.z)
'''
def CalculateLineLineIntersection(line1Point1, line1Point2, line2Point1, line2Point2):
    resultSegmentPoint1 = Point(0,0,0)
    resultSegmentPoint2 = Point(0,0,0)

    p1 = line1Point1
    p2 = line1Point2
    p3 = line2Point1
    p4 = line2Point2
    p13 = p1.minus(p3)
    p43 = p4.minus(p3)

    if p43.dist() ** 2 < delta ** 3:
        print('43')
        return None, None

    p21 = p2.minus(p1)
    if p21.dist() ** 2 < delta ** 3:
        print('21')
        return None, None

    d1343 = p13.x * p43.x + p13.y * p43.y + p13.y * p43.y
    d4321 = p43.x * p21.x + p43.y * p21.y + p43.y * p21.y
    d1321 = p13.x * p21.x + p13.y * p21.y + p13.y * p21.y
    d4343 = p43.x * p43.x + p43.y * p43.y + p43.y * p43.y
    d2121 = p21.x * p21.x + p21.y * p21.y + p21.y * p21.y

    denom = d2121 * d4343 - d4321 * d4321
    if abs(denom) < delta ** 3:
        print('denom')
        return None, None

    numer = d1343 * d4321 - d1321 * d4343

    mua = numer / denom
    mub = (d1343 + d4321 * mua) / d4343

    resultSegmentPoint1.x = p1.x + mua * p21.x
    resultSegmentPoint1.y = p1.y + mua * p21.y
    resultSegmentPoint1.z = p1.z + mua * p21.z
    resultSegmentPoint2.x = p3.x + mub * p43.x
    resultSegmentPoint2.y = p3.y + mub * p43.y
    resultSegmentPoint2.z = p3.z + mub * p43.z

    return resultSegmentPoint1, resultSegmentPoint2

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

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])'''
'''
def rotatetopfacenormal(edges, x, y, z, faces=None, debug=False):
    toppts = set()
    botpts = set()
    leftpts = None
    rightpts = None
    if faces is not None:
        leftpts, rightpts = faces
    for e in edges:
        e1, e2 = e
        if e1.z > 0: toppts.add(e1)
        else: botpts.add(e1)
        if e2.z > 0: toppts.add(e2)
        else: botpts.add(e2)
    #print(len(botpts))
    v1 = y.minus(x)
    v2 = z.minus(x)
    n = v1.cross(v2).scale()
    if x.minus(n).dist() > x.plus(n).dist():
        n = n.scale(-1)
    axis = n.cross(Point(0,0,1))
    theta = np.arccos(n.dot(Point(0,0,1)))-np.pi
    rm = rotation_matrix([axis.x, axis.y, axis.z], theta)

    newedges = set()
    newtoppts = set()
    newbotpts = set()
    newleftpts = set()
    newrightpts = set()
    for e in edges:
        e1, e2 = e
        e1r = np.dot(rm, [e1.x, e1.y, e1.z])
        e2r = np.dot(rm, [e2.x, e2.y, e2.z])
        newedges.add(frozenset([Point(e1r[0], e1r[1], e1r[2]), Point(e2r[0], e2r[1], e2r[2])]))
        if e1 in toppts: newtoppts.add(Point(e1r[0], e1r[1], e1r[2]))
        if e2 in toppts: newtoppts.add(Point(e2r[0], e2r[1], e2r[2]))
        if e1 in botpts: newbotpts.add(Point(e1r[0], e1r[1], e1r[2]))
        if e2 in botpts: newbotpts.add(Point(e2r[0], e2r[1], e2r[2]))
        if faces is not None:
            if e1 in leftpts: newleftpts.add(Point(e1r[0], e1r[1], e1r[2]))
            if e2 in leftpts: newleftpts.add(Point(e2r[0], e2r[1], e2r[2]))
            if e1 in rightpts: newrightpts.add(Point(e1r[0], e1r[1], e1r[2]))
            if e2 in rightpts: newrightpts.add(Point(e2r[0], e2r[1], e2r[2]))

    #newedges.add(frozenset([Point(0,0,0), n.scale(2)]))
    perp1 = x.minus(n.scale(x.dot(n) / n.dist() ** 2))
    toppt = list(newtoppts)[0]
    perp2 = toppt.minus(n.scale(toppt.dot(n) / n.dist() ** 2))
    theta = np.arccos(perp1.dot(perp2)/perp1.dist()/perp2.dist())
    if perp1.cross(perp2).dot(n) > 0:
        theta = -theta
    #print('theta', theta)
    rm = rotation_matrix([n.x, n.y, n.z], theta)

    newedges2 = set()
    newtoppts2 = set()
    newbotpts2 = set()
    newleftpts2 = set()
    newrightpts2 = set()
    for e in newedges:
        e1, e2 = e
        e1r = np.dot(rm, [e1.x, e1.y, e1.z])
        e2r = np.dot(rm, [e2.x, e2.y, e2.z])
        newedges2.add(frozenset([Point(e1r[0], e1r[1], e1r[2]), Point(e2r[0], e2r[1], e2r[2])]))
        if e1 in newtoppts: newtoppts2.add(Point(e1r[0], e1r[1], e1r[2]))
        if e2 in newtoppts: newtoppts2.add(Point(e2r[0], e2r[1], e2r[2]))
        if e1 in newbotpts: newbotpts2.add(Point(e1r[0], e1r[1], e1r[2]))
        if e2 in newbotpts: newbotpts2.add(Point(e2r[0], e2r[1], e2r[2]))
        if faces is not None:
            if e1 in newleftpts: newleftpts2.add(Point(e1r[0], e1r[1], e1r[2]))
            if e2 in newleftpts: newleftpts2.add(Point(e2r[0], e2r[1], e2r[2]))
            if e1 in newrightpts: newrightpts2.add(Point(e1r[0], e1r[1], e1r[2]))
            if e2 in newrightpts: newrightpts2.add(Point(e2r[0], e2r[1], e2r[2]))
    ptmap = {}
    for p1 in [x, y, z]:
        for p2 in newtoppts2:
            if p1 not in ptmap.keys():
                ptmap[p1] = p2
            else:
                ptmap[p1] = p2 if p1.dist(p2) < p1.dist(ptmap[p1]) else ptmap[p1]
    #print(x.x, x.y, x.z)
    #print(ptmap[x].x, ptmap[x].y, ptmap[x].z)
    #print(z.x, z.y, z.z)
    #print(ptmap[z].x, ptmap[z].y, ptmap[z].z)

    p1, p2 = CalculateLineLineIntersection(x, ptmap[x], z, ptmap[z])
    if p1 is None:
        p1, p2 = CalculateLineLineIntersection(y, ptmap[y], z, ptmap[z])
    if p1 is None:
        p1, p2 = CalculateLineLineIntersection(x, ptmap[x], y, ptmap[y])
    c = p1.midpoint(p2)
    r = np.sqrt(c.dist(y) * c.dist(ptmap[y]))
    newedges3 = set()
    newbotface = set()
    newleftface = set()
    newrightface = set()
    for e in newedges2:
        e1, e2 = e
        e1i = e1.inverse(c, r)
        e2i = e2.inverse(c, r)
        newedges3.add(frozenset([e1i, e2i]))
        if e1 in newbotpts2 and e2 in newbotpts2:
            newbotface.add(frozenset([e1i, e2i]))
        if e1 in newleftpts2 and e2 in newleftpts2:
            newleftface.add(frozenset([e1i, e2i]))
        if e1 in newrightpts2 and e2 in newrightpts2:
            newrightface.add(frozenset([e1i, e2i]))
    if debug:
        return newedges2, frozenset(newbotface)
    if faces is None:
        return newedges3, frozenset(newbotface)
    else:
        return newedges3, frozenset(newleftface), frozenset(newrightface)'''

'''
s = 1
antiprismcr = 0.822664388
prismcr = np.sqrt(7/12)
antiprismcrratio = projection_dist(antiprismcr*s)/antiprismcr/np.sqrt(2)
prismcrratio = projection_dist(prismcr*s)/prismcr/np.sqrt(3)
n=4
m=3
h = ((np.cos(np.pi/n)-np.cos(2*np.pi/n))/2)**(1/2)
antiprismpts = {Point(np.cos(k*np.pi/n)*antiprismcrratio, np.sin(k*np.pi/n)*antiprismcrratio, (-1)**k*h*antiprismcrratio)
                for k in range(2*n)}

prismpts = set()
for k in range(m):
    prismpts.add(Point(np.cos(k*2*np.pi/m)*prismcrratio, np.sin(k*2*np.pi/m)*prismcrratio, prismcrratio*np.sqrt(3)/2))
    prismpts.add(Point(np.cos(k*2*np.pi/m)*prismcrratio, np.sin(k*2*np.pi/m)*prismcrratio, -prismcrratio*np.sqrt(3)/2))

antiprismedges = set()
d = min({a.dist(b) for a in antiprismpts for b in antiprismpts if a != b})
for a in antiprismpts:
    for b in antiprismpts:
        if abs(a.dist(b)-d) < delta:
            antiprismedges.add(frozenset({a, b}))
#antiprismedges.add(frozenset([Point(0,0,0), Point(0,0,1)]))
d = min({a.dist(b) for a in prismpts for b in prismpts if a != b})
prismedges = set()
for a in prismpts:
    for b in prismpts:
        if abs(a.dist(b)-d) < delta:
            prismedges.add(frozenset({a,b}))

face1 = [] #square faces
face2 = []
face3 = []
for e in prismedges:
    e1, e2 = e
    if e1.x < 0 and e2.x < 0:
        face1.append(e)
    if e1.y > -e1.x/np.sqrt(3) and e2.y > -e2.x/np.sqrt(3):
        face2.append(e)
    if e1.y < e1.x/np.sqrt(3) and e2.y < e2.x/np.sqrt(3):
        face3.append(e)
prismfaces = {frozenset(face1), frozenset(face2), frozenset(face3)}

edgesall = set(prismedges)
newfaces = set()
for face in prismfaces:
    x, y = list(face)[0]
    z, _ = list(face)[1]
    for e in face:
        if x in e and y not in e:
            e1, e2 = e
            z = e1 if x != e1 else e2
    newedges, newface = rotatetopfacenormal(antiprismedges, x, y, z)
    newfaces.add(newface)
    edgesall = edgesall.union(newedges)

rm = rotation_matrix([0,1,0], np.pi/2)
hprismedges = set()
for e in prismedges:
    e1, e2 = e
    e1r = np.dot(rm, [e1.x, e1.y, e1.z])
    e2r = np.dot(rm, [e2.x, e2.y, e2.z])
    hprismedges.add(frozenset([Point(e1r[0], e1r[1], e1r[2]), Point(e2r[0], e2r[1], e2r[2])]))

hprismfaces = set()
prismfaces.remove(frozenset(face1))
for f in prismfaces:
    newface = set()
    for e in f:
        e1, e2 = e
        e1r = np.dot(rm, [e1.x, e1.y, e1.z])
        e2r = np.dot(rm, [e2.x, e2.y, e2.z])
        newface.add(Point(e1r[0], e1r[1], e1r[2]))
        newface.add(Point(e2r[0], e2r[1], e2r[2]))
    hprismfaces.add(frozenset(newface))

newfaces2 = set()
for face in newfaces:
    x, y = list(face)[0]
    z, _ = list(face)[1]
    for e in face:
        if x in e and y not in e:
            e1, e2 = e
            z = e1 if x != e1 else e2

    #print(x.x, x.y, x.z)
    #print(y.x, y.y, y.z)
    #print(z.x, z.y, z.z)
    newedges, newleftface, newrightface = rotatetopfacenormal(hprismedges, x, y, z, hprismfaces)
    newfaces2.add(newleftface)
    newfaces2.add(newrightface)
    edgesall = edgesall.union(newedges)
#print(len(newfaces2))
#p = list(list(list(newfaces2)[0])[1])[0]
#print(p.x, p.y, p.z)
for face in [list(newfaces2)[0]]:
    x, y = list(face)[0]
    z, _ = list(face)[1]
    for e in face:
        if x in e and y not in e:
            e1, e2 = e
            z = e1 if x != e1 else e2
    print(x.x, x.y, x.z)
    print(y.x, y.y, y.z)
    print(z.x, z.y, z.z)
    newedges, newface = rotatetopfacenormal(antiprismedges, x, y, z, None, True)
    newfaces.add(newface)
    edgesall = edgesall.union(newedges)

print(len(edgesall))
for edge in edgesall:
    c1, c2 = edge
    c3 = c2.inverse()
    if c3 is None:
        c3 = c1.inverse()
    c, r = threepointcircle(c1, c2, c3)
    mid = c1.midpoint(c2)
    if r == -1:
        c3 = Point(999999999, 999999999, 999999999)
    else:
        # finds midpoint of arc, so that threepointarc in solidworks has a third point to work with
        diff = mid.minus(c)
        sdiff = diff.scale(r/diff.dist())
        c3 = sdiff.plus(c)
    print(str(c1.x) + ", " + str(c1.y) + ", " + str(c1.z) + ", " +
          str(c2.x) + ", " + str(c2.y) + ", " + str(c2.z) + ", " +
          str(c3.x) + ", " + str(c3.y) + ", " + str(c3.z) + ", ")'''