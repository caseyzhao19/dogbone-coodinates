import numpy as np
from prism import Prism
from point import Point
from antiprism import Antiprism

import sympy
from sympy.solvers import solve
from sympy import Symbol

prism = Prism(0.5)
antiprism = Antiprism(0.5)


# points p and q, find hyperbolic distance
def distance(u, v):
    if u.dist() >= 1 or v.dist() >= 1:
        print('undefined distance between vectors u, v outside unit sphere')
        return None
    return 2 * sympy.log((u.minus(v).dist() + sympy.sqrt(u.dist() ** 2 * v.dist() ** 2
                                                         - 2 * u.dot(v) + 1)) /
                         (sympy.sqrt((1 - u.dist() ** 2) * (1 - v.dist() ** 2))))


'''x = Symbol('x')
print(solve(2 * sympy.log((x+sympy.sqrt(49*x**4/144-x**2/6+1))/
                          sympy.sqrt((1-7*x**2/12)*(1-7*x**2/12)))-0.5, x))
'''
scale = 0.5
'''x = scale
z = x * np.sqrt(3)/2

u = Point(x, 0, z)
v = Point(-x / 2, -np.sqrt(3) * x / 2, z)
w = Point(x, 0, -z)
print(distance(u, v))
print(distance(u, w))
print(distance(u, v)/distance(u, w))'''

# conclusion: prism is the same ratio despite scale when regular

# is this true for antiprisms?
scale = 0.01
x = scale
n = 4
h = ((np.cos(np.pi / n) - np.cos(2 * np.pi / n)) / 2) ** (1 / 2)
z = h * scale * 2
u = Point(x, 0, z)
v = Point(0, x, z)
w = Point(np.sqrt(2) * x / 2, -np.sqrt(2) * x / 2, -z)
print(distance(u, v))
print(distance(u, w))
print(distance(u, v) / distance(u, w))
print(u.dist(v) / u.dist(w))

prismD = distance(prism.pointsnamed[0], prism.pointsnamed[1])
print(antiprism.pointsnamed[0])
antiprismD = distance(antiprism.pointsnamed[0], antiprism.pointsnamed[1])
print(prismD)
print(antiprismD)
#print(prismD / antiprismD)
