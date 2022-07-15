from util import *
from prism import Prism
from point import Point
from antiprism import Antiprism

side_len = 0.00001
prism = Prism(side_len)
alledges = set(prism.edges)
alltrianglefaces = [[prism.pointsnamed[0], prism.pointsnamed[1], prism.pointsnamed[2]],
                    [prism.pointsnamed[3], prism.pointsnamed[4], prism.pointsnamed[5]]]
#for p, q in alledges: print(p, q)
antiprism = Antiprism(side_len)
antiprism = Antiprism(side_len)
squarefaces = [prism.face1, prism.face2, prism.face3]
newsquarefaces = []
for face in squarefaces:
    antiprism.moveface(2, 0, 4, prism.pointsnamed[face[0]], prism.pointsnamed[face[1]], prism.pointsnamed[face[2]])
    newsquarefaces.append([antiprism.pointsnamed[3], antiprism.pointsnamed[1], antiprism.pointsnamed[5]])
    alledges = alledges.union(antiprism.edges)
    for i in range(8):
        alltrianglefaces.append([antiprism.pointsnamed[i], antiprism.pointsnamed[(i+1) % 8],
                                 antiprism.pointsnamed[(i+2) % 8]])
for i in range(10):
    squarefaces = []
    for face in newsquarefaces:
        prism = Prism(side_len)
        prism.moveface(0, 3, 2, face[0], face[1], face[2])
        squarefaces.append([prism.pointsnamed[2], prism.pointsnamed[1], prism.pointsnamed[5]])
        squarefaces.append([prism.pointsnamed[1], prism.pointsnamed[0], prism.pointsnamed[4]])
        alledges = alledges.union(prism.edges)
        alltrianglefaces.append([prism.pointsnamed[0], prism.pointsnamed[1], prism.pointsnamed[2]])
        alltrianglefaces.append([prism.pointsnamed[3], prism.pointsnamed[4], prism.pointsnamed[5]])

    newsquarefaces = []
    for face in squarefaces:
        antiprism = Antiprism(side_len)
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
f = open("povcode.txt", "w")
for a, b, c in alltrianglefaces:
    d = a.inverse()
    center, radius = fourpointsphere(a,b,c,d)
    if radius == -1:
        print('hi')
        f.write("mesh {\n")
        f.write("\t triangle { <" + str(a.x) + "," + str(a.y) + "," + str(a.z) + ">, <" + str(b.x) + "," +
                str(b.y) + "," + str(b.z) + ">, <" + str(c.x) + "," + str(c.y) + "," + str(c.z) + ">}\n")
        f.write("\t\t texture{ pigment{ color rgb<1,0.7,0.1> }\n")
        f.write("\t\t\t finish { specular 1 }}\n")
        f.write("}\n")
    else:
        f.write("object {\n")
        f.write("\t sphere { <" + str(center.x) + "," + str(center.y) + "," + str(center.z) + ">," + str(radius) + "\n")
        f.write("\t\t texture{ pigment{ color rgb<1,0.7,0.1> }\n")
        f.write("\t\t\t finish { specular 1 }}}\n")
        f.write("\t clipped_by {mesh {\n")
        O = Point(0, 0, 0)
        for p,q,r in [(a,b,c), (O, a, b) , (O, a, c), (O, b, c)]:
            f.write("\t\t triangle { <" + str(p.x) + "," + str(p.y) + "," + str(p.z) + ">, <" + str(q.x) + "," +
                    str(q.y) + "," + str(q.z) + ">, <" + str(r.x) + "," + str(r.y) + "," + str(r.z) + "> }\n")
        f.write("\t\t inside_vector <0, 0, 1>}}\n")
        f.write("}\n")
f.close()
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