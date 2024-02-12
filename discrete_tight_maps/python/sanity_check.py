from make_linear_map import PL
import numpy as np
import matplotlib.pyplot as plt

# Sanity check:
# This initial data oscillates wildly on the interior, but agrees with the identity on the boundary.
# It should converge to a map close to the identity.
crazy = PL.extension([
    np.array([0, 0]),
    np.array([2, 0]),
    np.array([4, 0]),
    np.array([6, 0]),
    np.array([8, 0]),
    np.array([1, 1.732]),
    np.array([10, 9]),
    np.array([6, 0]),
    np.array([5, 0]),
    np.array([9, 1.732]),
    np.array([2, 3.464]),
    np.array([4, 0]),
    np.array([15, -2]),
    np.array([5, -10]),
    np.array([10, 3.464]),
    np.array([3, 5.196]),
    np.array([10, 19]),
    np.array([12, -10]),
    np.array([20, -2]),
    np.array([11, 5.196]),
    np.array([4, 6.928]),
    np.array([6, 6.928]),
    np.array([8, 6.928]),
    np.array([10, 6.928]),
    np.array([12, 6.928])
])

crazy.plot()

g = crazy.deform_tight(1, 0.01)
print("We improved the global Lipschitz constant from", max(crazy.localLip), "to", max(g.localLip))
print("Plotting...")
g.plot()

plt.show()