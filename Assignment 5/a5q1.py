import numpy as np 
from scipy import integrate
from matplotlib import pyplot as plt
import camb
from use_func import *


# Loading the data
planck = np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt',skiprows=1)

# Storing info from raw data
ell = planck[:,0]
spec = planck[:,1]
errs = 0.5 * (planck[:,2] + planck[:,3])


# Guess parameters:
pars = np.asarray([69, 0.022, 0.12, 0.06, 2.10e-9, 0.95])


# Generating predicted data
model = get_spectrum(pars)
model = model[:len(spec)]

# Calculating chi^2
chisq = chi2(model, spec, errs)

# Degrees of freedom
dof = len(spec)-len(pars)

print("The Chi^2 is:", chisq)
print("The number of degrees of freedom are:", dof)
print("Mean Chi^2 with error:", dof, "±", np.sqrt(2*dof))
