import numpy as np
gauss = np.random.normal
from scipy.linalg import svd

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# We start by defining the simplicial complex of the domain.
# Here are the vertices:
VERTICES = []
BOUNDARY = []

VERTICES_SORTED = []

ind = 0
for i in range(101):
    row = []
    for j in range(101):
        VERTICES.append([i/100, j/100])
        if any([i == 0, j == 0, i == 100, j == 100]):
            BOUNDARY.append(ind)
        ind = ind + 1
        row.append(ind)

SIMPLICES = []

for i in range(100):
    for j in range(100):
        SIMPLICES.append([
            VERTICES_SORTED[i][j],
            VERTICES_SORTED[i + 1][j],
            VERTICES_SORTED[i][j + 1]
        ])
        SIMPLICES.append([
            VERTICES_SORTED[i + 1][j],
            VERTICES_SORTED[i][j + 1],
            VERTICES_SORTED[i + 1][j + 1]
        ])

# TODO Everything under here

def linear_extn(x1, x2, y1, y2):
    # Returns the linear map (in two dimensions) which sends xi to yi (and 0 to 0).
    # We assume that x1 and x2 are linearly independent.
    A = np.array([[x1[0], x1[1]], [x2[0], x2[1]]])
    b = np.array([y1, y2])
    return(np.linalg.solve(A, b))

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
        self.localLip = [abs(M) for M in self.linear]
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
        if (self.globalLipOnVertices - other.globalLipOnVertices)/self.globalLipOnVertices < -0.0001:
            print("Found a map of smaller global Lip! Now", other.globalLipOnVertices)
            return True
        lhs = [self.localLip[s] if self.localLip[s] > other.localLip[s] else 0 for s in range(len(SIMPLICES))]
        rhs = [other.localLip[s] if other.localLip[s] > self.localLip[s] else 0 for s in range(len(SIMPLICES))]
        if max(lhs) < max(rhs) - 0.0001:
            print("Found a tighter map!")
            return True
        return False
    
    def perturb(self, temperature, localize):
        y = self.values.copy()
        for v in range(len(VERTICES)):
            if v not in BOUNDARY:
                local_temp = temperature
                if localize:
                    local_temp = self.localLipOnVertices[v] * local_temp / self.globalLipOnVertices
                y[v] = gauss(y[v], local_temp)
        return PL.extension(y)
    
    def plot(self):
        # We plot the map.
        x0 = min([self.values[v][0] for v in range(len(VERTICES))])
        x1 = max([self.values[v][0] for v in range(len(VERTICES))])
        y0 = min([self.values[v][1] for v in range(len(VERTICES))])
        y1 = max([self.values[v][1] for v in range(len(VERTICES))])

        fig, ax = plt.subplots()
        ax.set_xlim([x0, x1])
        ax.set_ylim([y0, y1])

        for s in range(len(SIMPLICES)):
            y = []
            for v in range(len(VERTICES)):
                if v in SIMPLICES[s]:
                    y.append(self.values[v])
            p = Polygon(np.array(y), facecolor=colors[s])
            ax.add_patch(p)

    def deform_tight(self, initial_temperature, final_temperature, intermediate_plot=False):
    # Returns the tight map which agrees with f on BOUNDARY
    # The temperature is the amount we perturb f at a "typical" step.
    # It exponentially decays, so we try smaller perturbations later on.
        candidate = self
        
        # print("First pass: trying to decrease the Lipschitz constant")
        time = 0
        most_recent_improvement = 0
        temperature = initial_temperature
        while temperature > final_temperature:
            if time % 2500 == 0:
                print("Current time:", time, " | Current temperature:", temperature)
                if intermediate_plot:
                    candidate.plot()
            if most_recent_improvement < time - 100:
                temperature = temperature * 0.999

            new_candidate = candidate.perturb(temperature, True)
            if new_candidate.tighter(candidate):
                candidate = new_candidate 
                most_recent_improvement = time
            time = time + 1
            
        # print("Second pass: trying to tighten")
        # time = 0
        # most_recent_improvement = 0
        # temperature = initial_temperature
        # while temperature > final_temperature:
        #     if time % 5000 == 0:
        #         print("Current time:", time, " | Current temperature:", temperature)
        #     if most_recent_improvement < time - 100:
        #         temperature = temperature * 0.999
        #     new_candidate = candidate.perturb(temperature, False)
        #     if new_candidate.tighter(candidate):
        #         candidate = new_candidate 
        #         most_recent_improvement = time
        #     time = time + 1

        return candidate