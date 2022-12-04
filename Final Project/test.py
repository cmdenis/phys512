import numpy as np
def f(x):
    x[2] = 3

w = np.ones(10)

f(w)

print(w)