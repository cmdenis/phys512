import numpy as np 
from matplotlib import pyplot as plt
import camb
from tqdm import tqdm
from use_func import *
import time
import datetime
import corner


# Loading the data
#planck = np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt', skiprows=1)
planck = np.loadtxt('COM_PowerSpect_CMB-TT-binned_R3.01.txt', skiprows=1)

# Storing info from raw data
ell = planck[:,0]
data_l = len(ell)
spec = planck[:,1]
errs = 0.5 * (planck[:,2] + planck[:,3])

# Guess parameters:
pars = np.asarray([69, 0.022, 0.12, 0.06, 2.10e-9, 0.95])
pars = np.asarray([6.68836115e+01, 2.20632912e-02, 1.19980709e-01, 2.11551525e+00, 1.89577396e-09, 9.66833369e-01])
#pars = np.array([6.88823057e+01,  3.26467430e-02,  1.23149023e-01,  5.57313513e-02, 8.99017735e-09,  9.54287992e-01])
#pars = np.array([6.89998762e+01, 3.00e-02, 1.23, 5.8e-02, 1.28348123e-09, 9.55e-01])

pred = get_spectrum(pars)[:data_l]


#plt.plot(ell, spec, label = "data")
plt.plot(ell, pred, label = "fit")
plt.show()
