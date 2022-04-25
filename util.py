import numpy as np
from point import Point

delta = 0.0001

# takes euclidean distance and finds ball distance but i think its only valid for diameters probably
# i found some stuff on wikipedia and did some un-good math for some special case for something
def projection_dist(d):
    return (np.exp(d) - 1)/(np.exp(d) + 1)

def projection_dist_inv(d):
    return np.log((1+d)/(1-d))

# finds line intersection between two points, that is, smallest line segment between two lines, returns
# endpoints of that smallest line segment
def CalculateLineLineIntersection(line1Point1, line1Point2, line2Point1, line2Point2):
    resultSegmentPoint1 = Point(0,0,0)
    resultSegmentPoint2 = Point(0,0,0)

    p1 = line1Point1
    p2 = line1Point2
    p3 = line2Point1
    p4 = line2Point2
    p13 = p1.minus(p3)
    p43 = p4.minus(p3)

    if p43.dist() ** 2 < delta ** 5:
        print('43')
        return None, None

    p21 = p2.minus(p1)
    if p21.dist() ** 2 < delta ** 5:
        print('21')
        return None, None

    d1343 = p13.x * p43.x + p13.y * p43.y + p13.y * p43.y
    d4321 = p43.x * p21.x + p43.y * p21.y + p43.y * p21.y
    d1321 = p13.x * p21.x + p13.y * p21.y + p13.y * p21.y
    d4343 = p43.x * p43.x + p43.y * p43.y + p43.y * p43.y
    d2121 = p21.x * p21.x + p21.y * p21.y + p21.y * p21.y

    denom = d2121 * d4343 - d4321 * d4321
    if abs(denom) < delta ** 5:
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


# Return the rotation matrix associated with counterclockwise rotation about
# the given axis by theta radians.
def rotation_matrix(axis, theta):
    axis = np.asarray(axis)
    if np.sqrt(np.dot(axis, axis)) == 0:
        return np.array([[1, 0, 0],
                         [0, 1, 0],
                         [0, 0, 1]])
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

# very definitely stolen am not this clevery
# if points are all on same plane, then sphere is undefined, use plane reflection instead
# return indicates by returning three points on plane and -1 for radius
def fourpointsphere(a, b, c, d):
    A = a.minus(d)
    B = b.minus(d)
    C = c.minus(d)
    # something something linear algebra
    det = A.x * (B.y * C.z - B.z * C.y) - A.y * (B.x * C.z - B.z * C.x) + A.z * (B.x * C.y - B.y * C.x)
    if abs(det) < delta**5:
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
    if abs(wsl) < delta ** 5:
        return None, -1
    iwsl2 = 0.5 / wsl
    center = a.plus((u.scale(t.dot(t)*u.dot(v)).minus(t.scale(u.dot(u)*t.dot(v)))).scale(iwsl2))  # i am sorry
    radius = np.sqrt(t.dot(t) * u.dot(u) * v.dot(v) * iwsl2 * 0.5)
    return center, radius

# 2d
def circleIntersect(c1, r1, c2, r2):
    # stolen also :(
    # i was having weird issues so i gave up
    d = c1.dist(c2)
    a = (r1 ** 2 - r2 ** 2 + d ** 2) / (2 * d)
    h = (r1 ** 2 - a ** 2) ** (1 / 2)

    x2 = c1.x + a * (c2.x - c1.x) / d
    y2 = c1.y + a * (c2.y - c1.y) / d

    x3 = x2 + h * (c2.y - c1.y) / d
    y3 = y2 - h * (c2.x - c1.x) / d
    x4 = x2 - h * (c2.y - c1.y) / d
    y4 = y2 + h * (c2.x - c1.x) / d

    return Point(x3, y3, 0), Point(x4, y4, 0)

# 2d
def perpBisect(a, b):
    c1, r1 = hyperCircle(a, b)
    c2, r2 = hyperCircle(b, a)
    p1, p2 = circleIntersect(c1, r1, c2, r2)
    return p1, p2

# 2d
def tangent(c, pt, l=10):
    slope = c.slope(pt)
    if slope == 0:
        p1 = Point(pt.x, pt.y + l, 0)
        p2 = Point(pt.x, pt.y - l, 0)
    elif slope is None:
        p1 = Point(pt.x-l, pt.y, 0)
        p2 = Point(pt.x+l, pt.y, 0)
    else:
        p1 = Point(pt.x-l, pt.y+l/slope, 0)
        p2 = Point(pt.x+l, pt.y-l/slope, 0)
    return p1, p2

def intersect(a,b,c,d):
    px = ((a.x*b.y-a.y*b.x)*(c.x-d.x)-(a.x-b.x)*(c.x*d.y-c.y*d.x))/(
                (a.x-b.x)*(c.y-d.y)-(a.y-b.y)*(c.x-d.x))
    py = ((a.x*b.y-a.y*b.x)*(c.y-d.y)-(a.y-b.y)*(c.x*d.y-c.y*d.x))/(
                (a.x-b.x)*(c.y-d.y)-(a.y-b.y)*(c.x-d.x))
    return Point(px, py, 0)

def hyperCircle(c, pt):
    if c == Point(0, 0, 0):
        C = c
        R = C.dist(pt)
    elif c.isCollinear(pt, Point(0, 0, 0)):
        ci = c.inverse()
        c2 = c.midpoint(ci)
        ai = pt.inverse(c2, c2.dist(c))
        C = pt.midpoint(ai)
        R = C.dist(ai)
    else:
        ct, rt = threepointcircle(c, pt, pt.inverse())
        t = tangent(ct, pt)
        C = intersect(t[0], t[1], Point(0, 0, 0), c)
        R = C.dist(pt)
    return C, R

# given two pts in 3d space find sphere that is perpendicular bisector
# this is probably vry stupid but i wanted to recycle the 2d code because im lazy
# basically maps 3d case to xy plane does the thingy and takes the inverse
def perpendicularbisector(a, b):
    c = Point(0,0,0)
    v1 = b.minus(a)
    if c.isCollinear(a,b):
        n = Point(0,-v1.z,v1.y).scale()
    else:
        v2 = c.minus(a)
        n = v1.cross(v2).scale()
    axis = n.cross(Point(0, 0, 1))
    theta = np.arccos(n.dot(Point(0, 0, 1)))
    rm = rotation_matrix([axis.x, axis.y, axis.z], theta)
    ar = np.dot(rm, [a.x, a.y, a.z])
    ar = Point(ar[0], ar[1], ar[2])
    br = np.dot(rm, [b.x, b.y, b.z])
    br = Point(br[0], br[1], br[2])
    c1, r1 = hyperCircle(ar, br)
    c2, r2 = hyperCircle(br, ar)
    p1, p2 = circleIntersect(c1, r1, c2, r2)
    C, R = threepointcircle(p1, p2, p1.inverse())
    C = np.dot(np.linalg.inv(rm), [C.x, C.y, C.z])
    return Point(C[0], C[1], C[2]), R

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

# given sphere with center C and radius R, and a,b,c on a circle, find intersections
# also terrible because i am lazy
# i need help
# with life
# and stuff
def spherecircleintersect(C, R, a, b, c):
    c, r = threepointcircle(a, b, c)
    v1 = b.minus(a)
    v2 = c.minus(a)
    n = v1.cross(v2).scale()
    d = n.dot(c.minus(C))
    ci = C.plus(n.scale(d))
    ri = np.sqrt(R ** 2 - d ** 2)
    axis = n.cross(Point(0, 0, 1))
    theta = np.arccos(n.dot(Point(0, 0, 1)))
    rm = rotation_matrix([axis.x, axis.y, axis.z], theta)
    c = np.dot(rm, [c.x, c.y, c.z])
    c = Point(c[0], c[1], c[2])
    ci = np.dot(rm, [ci.x, ci.y, ci.z])
    ci = Point(ci[0], ci[1], ci[2])
    p1, p2 = circleIntersect(c, r, ci, ri)
    p1 = np.dot(np.linalg.inv(rm), [p1.x, p1.y, p1.z])
    p2 = np.dot(np.linalg.inv(rm), [p2.x, p2.y, p2.z])
    return Point(p1[0], p1[1], p1[2]), Point(p2[0], p2[1], p2[2])


