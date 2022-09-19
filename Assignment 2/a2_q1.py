import numpy as np
import matplotlib.pyplot as plt
import scipy as sci

def ef_ring(z, Q = 1, R = 1):
    return Q*z/(4*np.pi*sci.constants.epsilon_0)

print(sci.constants.epsilon_0)

