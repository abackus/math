import math

def sharp(y):
  x = math.log(1 - y)
  return 1 - math.exp(x + 1)

def flat(y):
  x = math.log(1 - y)
  return 1 - math.exp(x - 1)

def value(n, y, beta):
  if n == 0:
    return max(y, 0)
  return max(y, 0, beta * (value(n - 1, sharp(y), beta) + value(n - 1, flat(y), beta)) / 2)
