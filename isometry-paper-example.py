from util import *
from prism import Prism
from point import Point
from antiprismdemo import Antiprism
import matplotlib.pyplot as plt

f = open("./povlengthall/povisoex.pov", "w")
color = ["<0.95, 1, 0>", "<0.6, 0.6, 0.6>","<1, 0.25, 0>"]
f.write("#include \"colors.inc\"\n")
f.write("\n")
f.write("camera {\n")
f.write("\t location <0, 0, " + str(-3) + ">\n")
f.write("\t look_at <0, 0, 0>\n")
f.write("}\n")
f.write("\n")
f.write("light_source {\n")
f.write("\t <0,0,-10> color <1, 1, 1>\n")
f.write("}\n")
f.write("\n")
f.write("light_source {\n")
f.write("\t <10,0,-10> color <0.3, 0.3, 0.3>\n")
f.write("}\n")
f.write("\n")
f.write("background { rgb <1,1,1> }")
f.write("\n")
cam_angle = 158

# origin
f.write("object {\n")
f.write("\t sphere { <0, 0, -0.5>, 0.01\n")
f.write("\t\t texture{ pigment{ color rgb <0.2, 0.1, 0.1>} rotate <" + str(cam_angle) + ", 0, 0>}\n")
f.write("\t\t}}\n")


lastface = 0
lastprism = 0
side_length = 0.7#0.02 + len_idx/frames
prism = Prism(side_length)
alledges = set(prism.edges)
# for p, q in alledges: print(p, q)
antiprism = Antiprism(side_length)
squarefaces = [prism.face2]
newsquarefaces = []
for face in squarefaces:
    antiprism.moveface(2, 0, 4, prism.pointsnamed[face[0]], prism.pointsnamed[face[1]],
                       prism.pointsnamed[face[2]])
    newsquarefaces.append([antiprism.pointsnamed[3], antiprism.pointsnamed[1], antiprism.pointsnamed[5]])
    alledges = alledges.union(antiprism.edges)
    for face in newsquarefaces:
        prism = Prism(side_length)
        prism.moveface(0, 3, 2, face[0], face[1], face[2])
        lastface = [prism.pointsnamed[1], prism.pointsnamed[0], prism.pointsnamed[4], prism.pointsnamed[3]]
        lastprism = prism

desttriangle = []
desttriangle.append([lastface[0], lastface[1], lastface[2]])
desttriangle.append([lastface[1], lastface[2], lastface[3]])
'''for a, b in [(0, 1), (0, 2), (2, 3), (1, 3)]:
    a = lastface[a]
    b = lastface[b]
    c = a.inverse()
    print(a)
    center, radius = threepointcircle(a, b, c)
    n = a.minus(c).cross(b.minus(c)).scale()
    default = Point(0, 1, 0)
    axis = n.cross(default)
    theta = np.arccos(n.dot(default))
    rm = rotation_matrix([axis.x, axis.y, axis.z], theta)
    f.write("object {")
    f.write("intersection {")
    f.write("object {")
    f.write("torus {" + str(radius) + ", 0.005} texture{pigment{ color rgb <0.05, 0.35, 0.5>}} matrix <")
    for row in [0, 1, 2]:
        for col in [0, 1, 2]:
            f.write(str(rm[row][col]))
            f.write(", ")
    f.write(str(center.x) + ", " + str(center.y) + ", " + str(center.z))
    center, radius = fourpointsphere(lastprism.pointsnamed[0], lastprism.pointsnamed[1],
                                     lastprism.pointsnamed[2], lastprism.pointsnamed[3])
    f.write(">} sphere {<" + str(center.x) + ", " + str(center.y) + ", " + str(center.z) + "> "
            + str(radius + 0.0025) + " pigment{ color rgb <0.05, 0.35, 0.5>}}")
    f.write("rotate <" + str(cam_angle) + ", 0, 0>}}\n")'''

antiprism = Antiprism(0.7)
antiprism.moveface(2, 0, 4, lastface[0], lastface[1], lastface[2], f, cam_angle)
triangles = []
for j in range(8):
    triangles.append([antiprism.pointsnamed[j], antiprism.pointsnamed[(j + 1) % 8],
                     antiprism.pointsnamed[(j + 2) % 8]])
#triangles.append([antiprism.pointsnamed[0], antiprism.pointsnamed[2], antiprism.pointsnamed[6]])
#triangles.append([antiprism.pointsnamed[2], antiprism.pointsnamed[4], antiprism.pointsnamed[6]])
triangles.append([antiprism.pointsnamed[1], antiprism.pointsnamed[5], antiprism.pointsnamed[7]])
triangles.append([antiprism.pointsnamed[1], antiprism.pointsnamed[3], antiprism.pointsnamed[5]])

for a, b, c in triangles:
    d = a.inverse()
    center, radius = fourpointsphere(a, b, c, d)
    f.write("object {\n")
    f.write(
        "\t sphere { <" + str(center.x) + "," + str(center.y) + "," + str(center.z) + ">," + str(radius) + "\n")
    f.write("\t\t texture{ pigment{ color rgb <1, 0.24, 0.1> filter 0.5}}\n")
    f.write("\t\t}\n")
    f.write("\t clipped_by {mesh {\n")
    O = Point(0, 0, 0)
    for p, q, r in [(a, b, c), (O, a, b), (O, a, c), (O, b, c)]:
        f.write("\t\t triangle { <" + str(p.x) + "," + str(p.y) + "," + str(p.z) + ">, <" + str(q.x) + "," +
                str(q.y) + "," + str(q.z) + ">, <" + str(r.x) + "," + str(r.y) + "," + str(r.z) + "> }\n")
    f.write("\t\t inside_vector <0, 0, 1>}}\n")
    f.write("rotate <" + str(cam_angle) + ", 0, 0>}\n")


for a, b, c in desttriangle:
    d = a.inverse()
    center, radius = fourpointsphere(a, b, c, d)
    f.write("object {\n")
    f.write(
        "\t sphere { <" + str(center.x) + "," + str(center.y) + "," + str(center.z) + ">," + str(radius) + "\n")
    f.write("\t\t texture{ pigment{ color rgb <0.1, 0.9, 1> filter 0.5}}\n")
    f.write("\t\t}\n")
    f.write("\t clipped_by {mesh {\n")
    O = Point(0, 0, 0)
    for p, q, r in [(a, b, c), (O, a, b), (O, a, c), (O, b, c)]:
        f.write("\t\t triangle { <" + str(p.x) + "," + str(p.y) + "," + str(p.z) + ">, <" + str(q.x) + "," +
                str(q.y) + "," + str(q.z) + ">, <" + str(r.x) + "," + str(r.y) + "," + str(r.z) + "> }\n")
    f.write("\t\t inside_vector <0, 0, 1>}}\n")
    f.write("rotate <" + str(cam_angle) + ", 0, 0>}\n")

for a, b in antiprism.edgesnamed:
    a = antiprism.pointsnamed[a]
    b = antiprism.pointsnamed[b]
    c = a.inverse()
    center, radius = threepointcircle(a, b, c)
    n = a.minus(c).cross(b.minus(c)).scale()
    default = Point(0, 1, 0)
    axis = n.cross(default)
    theta = np.arccos(n.dot(default))
    rm = rotation_matrix([axis.x, axis.y, axis.z], theta)
    f.write("intersection {")
    f.write("object {")
    f.write("torus {" + str(radius) + ", 0.005} texture{pigment{ color rgb <0.4, 0.05, 0.05>}} matrix <")
    for row in [0, 1, 2]:
        for col in [0, 1, 2]:
            f.write(str(rm[row][col]))
            f.write(", ")
    f.write(str(center.x) + ", " + str(center.y) + ", " + str(center.z))
    center, radius = fourpointsphere(antiprism.pointsnamed[0], antiprism.pointsnamed[1],
                                     antiprism.pointsnamed[2], antiprism.pointsnamed[3])
    f.write(">} sphere {<" + str(center.x) + ", " + str(center.y) + ", " + str(center.z) + "> "
            + str(radius) + " pigment{ color rgb <0.4, 0.05, 0.05>}} rotate <" + str(cam_angle) + ", 0, 0>}\n")
# unit sphere
f.write("object {\n")
f.write("\t sphere { <0, 0, 0,>, 1.00001\n")
f.write("\t\t texture{ pigment{ color rgb <1, 1, 1> filter 0.9}}\n")
f.write("\t\t}}\n")
f.close()

