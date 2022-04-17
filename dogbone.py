from util import *
from prism import Prism
from point import Point
from antiprism import Antiprism

prism = Prism(0.2)
alledges = set(prism.edges)
#for p, q in alledges: print(p, q)
antiprism = Antiprism(0.2)
squarefaces = [prism.face1, prism.face2, prism.face3]
newsquarefaces = []
for face in squarefaces:
    antiprism.moveface(2, 0, 4, prism.pointsnamed[face[0]], prism.pointsnamed[face[1]], prism.pointsnamed[face[2]])
    newsquarefaces.append([antiprism.pointsnamed[3], antiprism.pointsnamed[1], antiprism.pointsnamed[5]])
    alledges = alledges.union(antiprism.edges)
for i in range(5):
    squarefaces = []
    for face in newsquarefaces:
        prism.moveface(0, 3, 2, face[0], face[1], face[2])
        squarefaces.append([prism.pointsnamed[2], prism.pointsnamed[1], prism.pointsnamed[5]])
        squarefaces.append([prism.pointsnamed[1], prism.pointsnamed[0], prism.pointsnamed[4]])
        alledges = alledges.union(prism.edges)

    newsquarefaces = []
    for face in squarefaces:
        antiprism.moveface(2, 0, 4, face[0], face[1], face[2])
        newsquarefaces.append([antiprism.pointsnamed[3], antiprism.pointsnamed[1], antiprism.pointsnamed[5]])
        alledges = alledges.union(antiprism.edges)


'''antiprism.moveface(2, 0, 4, prism.pointsnamed[3], prism.pointsnamed[0], prism.pointsnamed[5])
alledges = alledges.union(antiprism.edges)
prism.moveface(0, 3, 2, antiprism.pointsnamed[3], antiprism.pointsnamed[1], antiprism.pointsnamed[5])
alledges = alledges.union(prism.edges)

prism = Prism()
antiprism.moveface(2, 0, 4, prism.pointsnamed[5], prism.pointsnamed[2], prism.pointsnamed[4])
alledges = alledges.union(antiprism.edges)
prism.moveface(0, 3, 2, antiprism.pointsnamed[3], antiprism.pointsnamed[1], antiprism.pointsnamed[5])
alledges = alledges.union(prism.edges)

prism = Prism()
antiprism.moveface(2, 0, 4, prism.pointsnamed[4], prism.pointsnamed[1], prism.pointsnamed[3])
alledges = alledges.union(antiprism.edges)
prism.moveface(0, 3, 2, antiprism.pointsnamed[3], antiprism.pointsnamed[1], antiprism.pointsnamed[5])
alledges = alledges.union(prism.edges)'''


print(len(alledges))
for edge in alledges:
    c1, c2 = edge
    c3 = c2.inverse()
    if c3 is None:
        c3 = c1.inverse()
    c, r = threepointcircle(c1, c2, c3)
    mid = c1.midpoint(c2)
    if True:#r == -1 or r>2:
        c3 = Point(999999999, 999999999, 999999999)
    else:
        # finds midpoint of arc, so that threepointarc in solidworks has a third point to work with
        diff = mid.minus(c)
        sdiff = diff.scale(r / diff.dist())
        c3 = sdiff.plus(c)
    print(str(c1.x) + ", " + str(c1.y) + ", " + str(c1.z) + ", " +
          str(c2.x) + ", " + str(c2.y) + ", " + str(c2.z) + ", " +
          str(c3.x) + ", " + str(c3.y) + ", " + str(c3.z) + ", ")