from make_linear_map import PL
import numpy as np
import matplotlib.pyplot as plt

# The boundary of this map is mapped to a self-crossing quadrilateral
crossing = PL.extension([
    np.array([0, 0]),
    np.array([2, 0]),
    np.array([4, 0]),
    np.array([6, 0]),
    np.array([8, 0]),
    np.array([3, 1.732]),
    np.array([4, 1.732]),
    np.array([5, 1.732]),
    np.array([6, 1.732]),
    np.array([6, 1.732]),
    np.array([6, 3.464]),
    np.array([5, 3.464]),
    np.array([5, 3.464]),
    np.array([5, 3.464]),
    np.array([4, 3.464]),
    np.array([9, 5.196]),
    np.array([7, 5.196]),
    np.array([5, 5.196]),
    np.array([3, 5.196]),
    np.array([2, 5.196]),
    np.array([12, 6.928]),
    np.array([9, 6.928]),
    np.array([6, 6.928]),
    np.array([3, 6.928]),
    np.array([0, 6.928])
])

from make_linear_map import VERTICES

crossing.plot()
print(["v: " + str(crossing.values[v]) for v in range(len(VERTICES))])

g = crossing.deform_tight(1, 0.01)
print("We improved the global Lipschitz constant from", max(crossing.localLip), "to", max(g.localLip))
print("Plotting...")
g.plot()

plt.show()