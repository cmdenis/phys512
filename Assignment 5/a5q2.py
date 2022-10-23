from black import out
import numpy as np 
from scipy import integrate
from matplotlib import pyplot as plt
import camb
from use_func import *
import time
import datetime
import corner


# Guess parameters
pars = np.asarray([69, 0.022, 0.12, 0.06, 2.10e-9, 0.95])

# Loading the data
planck = np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt', skiprows=1) # Importing raw data
x_data = planck[:,0]    # X axis data
y_data = planck[:,1]    # Y axis data
errs = 0.5 * (planck[:,2] + planck[:,3])    # Errors on data

data_l = len(x_data)    # Length of data

# Finding the length of synthsized data
sim_len = len(get_spectrum(pars))


# Defining gradient taking function
def num_derivs(fun,pars,dp):
    A = np.empty([sim_len, len(pars)])  # Initializing the grad
    for i in range(len(pars)):
        pp = pars.copy()
        pp[i] = pars[i]+dp[i]   # Parameters on right side
        y_right = fun(pp)       # Right side of derivative
        pp[i] = pars[i]-dp[i]   # Parameters on left side
        y_left = fun(pp)        # Left side of derivative
        
        A[:,i]=((y_right-y_left)/2/dp[i])   # Computing rows of gradient
    return A

# Defining Newton method
def num_newton(fun, pars, dp,x , y, sigma = 1, niter = 5):
    
    chi2prev = 1E12 # Initial Large Chi^2

    inv_N = np.eye(len(sigma))/sigma**2 # Inverse of error matrix

    # Main Newton Method loop
    for i in range(niter):
    
        pred = fun(pars)[0:data_l]              # Predicted data
        r = y - pred                            # Residuals
        A = num_derivs(fun, pars, dp)[0:data_l] # Gradient

        lhs = A.T@inv_N@A                       # Finding lin alg rhs
        rhs = A.T@inv_N@r                       # Finding lin alg lhs
        step = np.linalg.inv(lhs)@rhs           # Finding size and direction   
        
        pars = pars + step                      # New position after step
 
        
        chi2 =  np.sum(r**2/sigma**2)           # Computing Chi^2

        # Printing stuff
        print("Chi^2 is:", chi2)              
        print("Chi^2 Diff is:", chi2prev - chi2)
        print("Parameters are:", pars)
        print('Step is:', step)

        chi2prev = chi2     # Updating Chi^2
        
    print("Done Newton Method\n")
    return pars, np.linalg.inv(lhs)

# Chi^2 function
def chisq(y, pred, err):
    r=y-pred
    return np.sum(r**2/err**2)


# Doing Newton Method
fun = get_spectrum
print("Starting Newton Method...\n")
pars, cov_mat = num_newton(fun, pars, pars*1e-8, x_data, y_data, errs, 1)
np.savetxt("planck_fit_params.txt", pars) # Storing params on cpu
np.savetxt("planck_fit_cov.txt", cov_mat) # Storing errors on cpu

print("The best fit parameters from Newton are:", pars)
print("Their error is:", np.sqrt(np.diagonal(cov_mat)))

