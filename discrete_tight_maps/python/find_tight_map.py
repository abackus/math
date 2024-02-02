import numpy
gauss = numpy.random.normal
array = numpy.array
from scipy.linalg import svd

TRIALS = 10**2

VERTICES = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.5, 0.5]]

# We assume that, if the BOUNDARY_FUNCTION has a domain of size n, then the first n VERTICES are on the boundary.
BOUNDARY_FUNCTION = [[0.0, 0.0], [2.0, 0.0], [0.0, 1.0], [3.0, 2.0]]

class Simplex:
    # A Simplex is a list of vertices.
    # Instead of working with the vertices themselves, we just work with their indices.
    def __init__(self, vertices):
        self.vertices = vertices

SIMPLICES = [
    Simplex([0, 1, 4]),
    Simplex([0, 2, 4]),
    Simplex([3, 1, 4]),
    Simplex([3, 2, 4])
]

class Function:
    # A function is a point assigned to each vertex.
    def __init__(self, values):
        self.values = values 
        self.lex = [self.Lip(s) for s in SIMPLICES].sort().reverse()

    def Lip(self, s):
        # Computes the Lipschitz constant of the PL extension of self to the simplex s.

        # We write the PL extension in point-slope form: u = y1 + A(x - x1)
        x1 = VERTICES[s.vertices[0]]
        y1 = self.values[self.vertices[0]]
        for i in range(1, len(s)):
            # Let A be the matrix...
            # TODO
            A = 0

        return svd(A)[1][0]

    def perturb(self, sigma):
        # Create a copy of self with Gaussian noise of standard deviation sigma added to each internal vertex.
        new_values = self.values.copy()
        for i in range(len(BOUNDARY_FUNCTION), len(VERTICES)):
            new_values[i] = new_values[i] + gauss(0, sigma)
        return Function(new_values)

    def __lt__ (self, other):
        for i in range(len(SIMPLICES)):
            if self.lex[i] < other.lex[i]:
                return True
        return False
    
def tight_map():
    initial_values = BOUNDARY_FUNCTION.copy()
    for _ in range(len(VERTICES) - len(BOUNDARY_FUNCTION)):
        initial_values.append(0.0)
    candidate = Function(initial_values)

    for i in range(TRIALS):
        new_candidate = candidate.perturb(2 ** (0 - i(TRIALS/10)))
        if new_candidate < candidate:
            candidate = new_candidate 
    
    return candidate

print("Lipschitz constant of the tight map is", tight_map().lex[0])