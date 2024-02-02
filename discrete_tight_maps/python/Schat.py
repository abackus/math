import numpy
gauss = numpy.random.normal
array = numpy.array
from scipy.linalg import svd

def Schat(x, y):
    m1 = array([[2, 2 * (x - 1)], [0, 2 * y]])
    m2 = array([[2 * x, 0], [2 * x - 1, 1]])
    m3 = array([[3, 3 - 2 * x], [3 * (y - 1), 1]])
    m4 = array([[2 * x, 5 - 2 * x], [2 * (1 - y), 2]])    
    Q = max(svd(m1)[1][0], svd(m2)[1][0], svd(m3)[1][0], svd(m4)[1][0])
    return Q


minx = 0
miny = 0
min_Schat = Schat(minx, miny)

trials = 10**2

for i in range(trials):
    sigma = 2**(0 - i/(trials/10))
    candidatex = gauss(minx, sigma)
    candidatey = gauss(miny, sigma)
    candidate_Schat = Schat(candidatex, candidatey)
    if candidate_Schat < min_Schat:
        minx = candidatex 
        miny = candidatey
        min_Schat = candidate_Schat
    if i % (trials/10) == 0:
        print("Current best Lipschitz constant is", min_Schat, ". Current standard variation is", sigma, ".")

print("Minimizing vertex is (", minx, ",", miny, "). Minimizing Lipschitz constant is", min_Schat, ".")