def J(p, t):
    '''Computes the pairing function J.'''
    return (int((p+t) * (p + t + 1) * 0.5) + p) # (p+t)*(p+t+1) is even so this makes sense

def pairing_projection(n):
    '''Computes the projections p, t such that J(p, t) = n'''
    i = 0 # plays the role of p + t
    while True: # will halt because J is surjective
        for t in range(0, i + 1):
            p = i - t
            if J(p, t) == n:
                return p, t
            i = i + 1
