import numpy as np
gauss = np.random.normal
from scipy.linalg import svd

# We start by defining the simplicial complex of the domain.
# Here are the vertices:
VERTICES = [
    np.array([0, 0]),
    np.array([2, 0]),
    np.array([4, 0]),
    np.array([6, 0]),
    np.array([8, 0]),
    np.array([1, 1.732]),
    np.array([3, 1.732]),
    np.array([5, 1.732]),
    np.array([7, 1.732]),
    np.array([9, 1.732]),
    np.array([2, 3.464]),
    np.array([4, 3.464]),
    np.array([6, 3.464]),
    np.array([8, 3.464]),
    np.array([10, 3.464]),
    np.array([3, 5.196]),
    np.array([5, 5.196]),
    np.array([7, 5.196]),
    np.array([9, 5.196]),
    np.array([11, 5.196]),
    np.array([4, 6.928]),
    np.array([6, 6.928]),
    np.array([8, 6.928]),
    np.array([10, 6.928]),
    np.array([12, 6.928])
]

# A simplex consists of a triple of integers, which are the indices of the vertices.
SIMPLICES = [
    [0, 5, 1],
    [5, 1, 6],
    [1, 6, 2],
    [6, 2, 7],
    [2, 7, 3],
    [7, 3, 8],
    [8, 4, 9],
    [5, 10, 6],
    [10, 6, 11],
    [6, 11, 7],
    [11, 7, 12],
    [7, 12, 8],
    [12, 8, 13],
    [8, 13, 9],
    [13, 9, 14],
    [10, 15, 11],
    [15, 11, 16],
    [11, 16, 12],
    [16, 12, 17],
    [12, 17, 13],
    [17, 13, 18],
    [18, 14, 19],
    [15, 20, 16],
    [20, 16, 21],
    [16, 21, 17],
    [21, 17, 22],
    [17, 22, 18],
    [22, 18, 23],
    [18, 23, 19],
    [23, 19, 24]
]

# A list of indices of vertices which are "fixed".
BOUNDARY = [
    0,
    1,
    2,
    3,
    4,
    5,
    9,
    10,
    14,
    15,
    19,
    20,
    21,
    22,
    23,
    24
]

def link(v):
    # Returns points in the link of a vertex of index v 
    pts = []
    for simplex in SIMPLICES:
        if v in simplex:
            w = simplex.copy()
            w.remove(v)
            for i in range(100):
                pts.append((VERTICES[w[0]] * (i/100)) + (VERTICES[w[1]] * ((100 - i)/100)))
    return pts

def linear_extn(x1, x2, y1, y2):
    # Returns the linear map (in two dimensions) which sends xi to yi (and 0 to 0).
    # We assume that x1 and x2 are linearly independent.
    A = np.array([[x1[0], x1[1], 0, 0], [0, 0, x1[0], x1[1]], [x2[0], x2[1], 0, 0], [0, 0, x2[0], x2[1]]])
    b = np.array([y1[0], y1[1], y2[0], y2[1]])
    x = np.linalg.solve(A, b)
    return np.array([[x[0], x[1]], [x[2], x[3]]])

def affine_extn(x0, x1, x2, y0, y1, y2):
    # Returns the affine map (in two dimensions) which sends xi to yi.
    # We assume that x0, x1, x2 are in general position.
    M = linear_extn(x1 - x0, x2 - x0, y1 - y0, y2 - y0)
    b = y0 - np.matmul(M, x0)
    return M, b

class PL:
    def __init__(self, linear, constant):
        self.linear = linear
        self.constant = constant
        self.localLip = [svd(M)[1][0] for M in self.linear]
        self.localLipOnVertices = [
            max([
                self.localLip[s] for s in range(len(SIMPLICES)) if v in SIMPLICES[s]
            ]) for v in range(len(VERTICES))
        ]
        self.globalLipOnVertices = max(self.localLipOnVertices)

        self.values = []
        for v in range(len(VERTICES)):
            for s in range(len(SIMPLICES)):
                if v in SIMPLICES[s]:
                    self.values.append(np.matmul(self.linear[s], VERTICES[v]) + self.constant[s])
                    break
                
    def extension(y):
        # Returns the PL map (in two dimensions) for the simplices, which sends the ith vertex to y[i].
        Ms = []
        bs = []
        for s in SIMPLICES:
            M, b = affine_extn(VERTICES[s[0]], VERTICES[s[1]], VERTICES[s[2]], y[s[0]], y[s[1]], y[s[2]])
            Ms.append(M)
            bs.append(b)
        return PL(Ms, bs)

    def tighter(self, other):
        if self.globalLipOnVertices == 0:
            return True
        if (self.globalLipOnVertices - other.globalLipOnVertices)/self.globalLipOnVertices < -0.1:
            return True
        lhs = [self.localLip[s] if self.localLip[s] > other.localLip[s] else 0 for s in range(len(SIMPLICES))]
        rhs = [other.localLip[s] if other.localLip[s] > self.localLip[s] else 0 for s in range(len(SIMPLICES))]
        return (max(lhs) < max(rhs))
    
    def perturb(self, temperature):
        y = self.values.copy()
        for v in range(len(VERTICES)):
            if v not in BOUNDARY:
                local_temp = self.localLipOnVertices[v] * temperature / self.globalLipOnVertices
                y[v] = gauss(y[v], local_temp)
        return PL.extension(y)

def find_tight_map(f, initial_temperature, final_temperature):
    # Returns the tight map which agrees with f on BOUNDARY
    # The temperature is the amount we perturb f at a "typical" step.
    # It exponentially decays, so we try smaller perturbations later on.
    candidate = f
    temperature = initial_temperature
    i = 0
    while temperature > final_temperature:
        new_candidate = candidate.perturb(temperature)
        if new_candidate.tighter(candidate):
            candidate = new_candidate 
            print("Improvement at step", i, ": Temperature", temperature, ": Lip", candidate.globalLipOnVertices)
        temperature = 0.999 * temperature
        i = i + 1
    return candidate

identity = PL.extension(VERTICES)

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

g = find_tight_map(crazy, 1, 0.01)
print("Lipschitz constant of the tight competitor:", max(g.localLip))
print("Lipschitz constant of the wildly oscillating map:", max(crazy.localLip))