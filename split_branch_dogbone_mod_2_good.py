from util import *
from prism import Prism
from point import Point
from antiprism import Antiprism
import matplotlib.pyplot as plt

frames = 100
for len_idx in range(150):  # range(frames):
    if len_idx < 100:
        continue
    side_length = 0.02 + len_idx/frames
    prism = Prism(side_length)
    alledges = set(prism.edges)
    # for p, q in alledges: print(p, q)
    antiprism = Antiprism(side_length)
    squarefaces = [prism.face2]
    newsquarefaces = []
    # initialize square faces center prism
    branches = [[prism.pointsnamed[1], prism.pointsnamed[2], prism.pointsnamed[4]],
                [prism.pointsnamed[2], prism.pointsnamed[4], prism.pointsnamed[5]]]
    for face in squarefaces:
        antiprism.moveface(2, 0, 4, prism.pointsnamed[face[0]], prism.pointsnamed[face[1]],
                           prism.pointsnamed[face[2]])
        newsquarefaces.append([antiprism.pointsnamed[3], antiprism.pointsnamed[1], antiprism.pointsnamed[5]])
        alledges = alledges.union(antiprism.edges)
        for i in range(8):
            branches.append([antiprism.pointsnamed[i], antiprism.pointsnamed[(i + 1) % 8],
                                     antiprism.pointsnamed[(i + 2) % 8]])
    for i in range(1):
        squarefaces = []
        for face in newsquarefaces:
            prism = Prism(side_length)
            prism.moveface(0, 3, 2, face[0], face[1], face[2])
            squarefaces.append(
                [prism.pointsnamed[2], prism.pointsnamed[1], prism.pointsnamed[5], prism.pointsnamed[4]])
            squarefaces.append(
                [prism.pointsnamed[1], prism.pointsnamed[0], prism.pointsnamed[4], prism.pointsnamed[3]])
            alledges = alledges.union(prism.edges)
            branches.append([prism.pointsnamed[0], prism.pointsnamed[1], prism.pointsnamed[2]])
            branches.append([prism.pointsnamed[3], prism.pointsnamed[4], prism.pointsnamed[5]])

        newsquarefaces = []
        for face in squarefaces:
            antiprism = Antiprism(side_length)
            antiprism.moveface(2, 0, 4, face[0], face[1], face[2])
            newsquarefaces.append(
                [antiprism.pointsnamed[3], antiprism.pointsnamed[1], antiprism.pointsnamed[5],
                 antiprism.pointsnamed[7]])
            alledges = alledges.union(antiprism.edges)
            for j in range(8):
                branches.append([antiprism.pointsnamed[j], antiprism.pointsnamed[(j + 1) % 8],
                                         antiprism.pointsnamed[(j + 2) % 8]])
    disttowalls = []
    for i in range(12):
        curritertriangles = []
        squarefaces = []
        for face in newsquarefaces:
            prism = Prism(side_length)
            prism.moveface(2, 0, 5, face[0], face[1], face[2])
            squarefaces.append(
                [prism.pointsnamed[2], prism.pointsnamed[1], prism.pointsnamed[5], prism.pointsnamed[4]])
            squarefaces.append(
                [prism.pointsnamed[1], prism.pointsnamed[0], prism.pointsnamed[4], prism.pointsnamed[3]])
            alledges = alledges.union(prism.edges)
            curritertriangles.append([prism.pointsnamed[0], prism.pointsnamed[1], prism.pointsnamed[2]])
            curritertriangles.append([prism.pointsnamed[3], prism.pointsnamed[4], prism.pointsnamed[5]])

        newsquarefaces = []
        for face in squarefaces:
            antiprism = Antiprism(side_length)
            antiprism.moveface(2, 0, 4, face[0], face[1], face[2])
            newsquarefaces.append(
                [antiprism.pointsnamed[3], antiprism.pointsnamed[1], antiprism.pointsnamed[5],
                 antiprism.pointsnamed[7]])
            alledges = alledges.union(antiprism.edges)
            for j in range(8):
                curritertriangles.append([antiprism.pointsnamed[j], antiprism.pointsnamed[(j + 1) % 8],
                                 antiprism.pointsnamed[(j + 2) % 8]])
        branches.extend(curritertriangles)
        minwalldist = 1000
        for triangle in curritertriangles:
            for p in triangle:
                minwalldist = min(minwalldist, pointplanedistance(p.x, p.y, p.z, -1, -0.5, 0, 0))
                minwalldist = min(minwalldist, pointplanedistance(p.x, p.y, p.z, -1, 0.5, 0, 0))
        disttowalls.append(minwalldist)

    for face in newsquarefaces:
        branches.append([face[0], face[1], face[2]])
        branches.append([face[1], face[2], face[3]])
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
    max_dist = 0
    for a, b, c in branches:
        max_dist = max(max_dist, a.dist())
        max_dist = max(max_dist, b.dist())
        max_dist = max(max_dist, c.dist())

    print(disttowalls)
    figure, axes = plt.subplots()
    plt.figure(figsize=(4, 7))
    axes.set_ylabel('distance')
    axes.set_ylabel('iteration')
    plt.title("side length = " + "{:.2f}".format(side_length))
    plt.axhline(0, color='lightgray')
    plt.xlim([0.5, 12.5])
    plt.ylim([-0.8, 0.8])
    plt.xticks(range(1, 13))
    plt.yticks(np.arange(-0.8, 0.9, 0.1))
    plt.plot(range(1, 13), disttowalls)
    plt.savefig('./plots/plot' + str(len_idx) + '.jpg', dpi=600)


    print(len(alledges))
    str_len_idx = str(len_idx)
    if len_idx < 10:
        str_len_idx = "00" + str(len_idx)
    elif len_idx < 100:
        str_len_idx = "0" + str(len_idx)
    f = open("./povlengthall/povcodeall" + str_len_idx + ".pov", "w")
    color = ["<0.95, 1, 0>", "<0.6, 0.6, 0.6>","<1, 0.25, 0>"]
    f.write("#include \"colors.inc\"\n")
    f.write("\n")
    f.write("camera {\n")
    f.write("\t location <0, 0, " + str(-3 * max_dist) + ">\n")
    f.write("\t look_at <0, 0, 0>\n")
    f.write("}\n")
    f.write("\n")
    f.write("light_source {\n")
    f.write("\t <0,0,-10> color <0.3, 0.3, 0.3>\n")
    f.write("}\n")
    f.write("\n")
    f.write("light_source {\n")
    f.write("\t <10,0,-10> color <0.3, 0.3, 0.3>\n")
    f.write("}\n")
    f.write("\n")
    f.write("background { rgb <0.1,.3,.4> }")
    f.write("\n")
    f.write("union {")
    for branch in [0]:
        f.write("object {\n")
        f.write("\t mesh {\n")
        for a, b, c in branches:
            f.write("\t\t triangle { <" + str(a.x) + "," + str(a.y) + "," + str(a.z) + ">, <" + str(b.x) + "," +
                    str(b.y) + "," + str(b.z) + ">, <" + str(c.x) + "," + str(c.y) + "," + str(c.z) + "> }\n")
        f.write("\t\t inside_vector <0, 0, 1>\n")
        f.write("\t\t texture{ pigment{ color rgb" + color[branch] + "}\n")
        f.write("\t\t\t finish { specular 1 }}\n")
        f.write("\t}\n")
        f.write("rotate <0, 0," + str(120 * branch) + ">}\n")
    f.write("rotate <0, 0, 0>")
    f.write("}")
    f.write("plane { <0, 1, 0>, 4 }")
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