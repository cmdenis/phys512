import numpy as np
import os
import corner
import matplotlib.pyplot as plt
import random
from use_func import *

# Here we use importance sampling with our original chain

# Load data
chain = np.loadtxt("runs3/planck_chain.txt")

# Calculating weights according to gaussian distribution with mean and std given in Q4 
weights = np.exp(-(((0.054-chain[:, 3])/0.0074)**2)/2)/(0.0074*np.sqrt(2*np.pi))

importance_avg = np.average(chain, weights=weights, axis = 0)

print("The new parameters with importance sampling are:", importance_avg)
