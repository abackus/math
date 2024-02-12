from make_linear_map import PL
import numpy as np
import matplotlib.pyplot as plt

# The boundary of this map stretches the quadrilateral to a very scalene triangle
crossing = PL.extension([
    np.array([0, 0]), # 0
    np.array([2, 0]), # 1
    np.array([4, 0]), # 2
    np.array([6, 0]), # 3
    np.array([8, 0]), # 4

    np.array([1, 5]), # 5
    np.array([2.5, 5]), # 6
    np.array([4, 5]), # 7
    np.array([5.5, 5]), # 8
    np.array([7, 5]), # 9

    np.array([2, 10]), # 10
    np.array([3, 10]), # 11
    np.array([4, 10]), # 12
    np.array([5, 10]), # 13
    np.array([6, 10]), # 14

    np.array([3, 15]), # 15
    np.array([3.5, 15]), # 16
    np.array([4, 15]), # 17
    np.array([4.5, 15]), # 18
    np.array([5, 15]), # 19

    np.array([4, 20]), # 20
    np.array([4, 20]), # 21
    np.array([4, 20]), # 22
    np.array([4, 20]), # 23
    np.array([4, 20]) # 24
])


crossing.plot()

g = crossing.deform_tight(0.1, 0.0005)
print("We improved the global Lipschitz constant from", max(crossing.localLip), "to", max(g.localLip))
print("Plotting...")
g.plot()

plt.show()