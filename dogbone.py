from util import *
from prism import Prism
from point import Point
from antiprism import Antiprism

prism = Prism(2)
alledges = set(prism.edges)
alltrianglefaces = [[prism.pointsnamed[0], prism.pointsnamed[1], prism.pointsnamed[2]],
                    [prism.pointsnamed[3], prism.pointsnamed[4], prism.pointsnamed[5]]]
#for p, q in alledges: print(p, q)
antiprism = Antiprism(2)
squarefaces = [prism.face1, prism.face2, prism.face3]
newsquarefaces = []
for face in squarefaces:
    antiprism.moveface(2, 0, 4, prism.pointsnamed[face[0]], prism.pointsnamed[face[1]], prism.pointsnamed[face[2]])
    newsquarefaces.append([antiprism.pointsnamed[3], antiprism.pointsnamed[1], antiprism.pointsnamed[5]])
    alledges = alledges.union(antiprism.edges)
    for i in range(8):
        alltrianglefaces.append([antiprism.pointsnamed[i], antiprism.pointsnamed[(i+1) % 8],
                                 antiprism.pointsnamed[(i+2) % 8]])
for i in range(0):
    squarefaces = []
    for face in newsquarefaces:
        prism.moveface(0, 3, 2, face[0], face[1], face[2])
        squarefaces.append([prism.pointsnamed[2], prism.pointsnamed[1], prism.pointsnamed[5]])
        squarefaces.append([prism.pointsnamed[1], prism.pointsnamed[0], prism.pointsnamed[4]])
        alledges = alledges.union(prism.edges)
        alltrianglefaces.append([prism.pointsnamed[0], prism.pointsnamed[1], prism.pointsnamed[2]])
        alltrianglefaces.append([prism.pointsnamed[3], prism.pointsnamed[4], prism.pointsnamed[5]])

    newsquarefaces = []
    for face in squarefaces:
        antiprism.moveface(2, 0, 4, face[0], face[1], face[2])
        newsquarefaces.append([antiprism.pointsnamed[3], antiprism.pointsnamed[1], antiprism.pointsnamed[5]])
        alledges = alledges.union(antiprism.edges)
        for j in range(8):
            alltrianglefaces.append([antiprism.pointsnamed[j], antiprism.pointsnamed[(j+1) % 8],
                                     antiprism.pointsnamed[(j+2) % 8]])

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
print(len(alltrianglefaces))

for a, b, c in alltrianglefaces:
    d = a.inverse()
    center, radius = fourpointsphere(a,b,c,d)
    print("intersection {")
    print("\t sphere { <", center.x, ",", center.y, ",", center.z, ">,", radius)
    print("\t\t texture{ pigment{ color rgb<1,0.7,0.1> }")
    print("\t\t\t finish { specular .7 }}} ")
    print("\t clipped_by {mesh {")
    O = Point(0, 0, 0)
    for p,q,r in [(a,b,c), (O, a, b) , (O, a, c), (O, b, c)]:
        print("\t\t triangle { <", p.x, ",", p.y, ",", p.z, ">, <", q.x, ",", q.y, ",", q.z, ">, <",
              r.x, ",", r.y, ",", r.z, "> }")
    print("\t\t inside_vector <0, 0, 1>}}")
    print("}")

'''for edge in alledges:
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
          str(c3.x) + ", " + str(c3.y) + ", " + str(c3.z) + ", ")'''