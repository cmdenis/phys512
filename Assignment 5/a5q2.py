import numpy as np 
from scipy import integrate
from matplotlib import pyplot as plt
import camb
from tqdm import tqdm
from use_func import *
import time
import datetime


# Loading the data
planck = np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt', skiprows=1)

# Storing info from raw data
ell = planck[:,0]
data_l = len(ell)
spec = planck[:,1]
errs = 0.5 * (planck[:,2] + planck[:,3])

# Guess parameters:
pars = np.asarray([69, 0.022, 0.12, 0.06, 2.10e-9, 0.95])



# Starting the initial loop conditions
looper = True
pred = get_spectrum(pars)[:data_l]
old_chi = chi2(pred, spec, errs)
print("Initial Chi^2:", old_chi)
trunc = 0
'''
# We do the procedure until chi^2 condition is satisfied
while looper:
    
    # Calculating predicted fit and gradient
    grad = grad_f(get_spectrum, pars, length = data_l)

    # Residuals
    r = spec - pred

    # Fragmenting linear algera stuff
    lhs = grad.T@grad
    rhs = grad.T@r 

    # Applying a step
    dp = np.linalg.inv(lhs)@rhs
    pars = pars + dp

    pred = get_spectrum(pars)[:data_l]
    new_chi = chi2(pred, spec, errs)

    # Check if if Chi^2 is minimized enough
    if old_chi > new_chi:
        if (old_chi - new_chi) < 0.01:
            looper = False
        else:
            old_chi = new_chi
    else:
        print("\nChi square not minimized... Step made it bigger...\n")
        looper = False
    print("\nParameters:", pars, "\nDisplacement", dp)
    print("Chi^2:", new_chi)
'''

'''Now the MCMC part'''

print("\n Now MCMC part \n")

newton_res = np.asarray([69, 0.022, 0.12, 0.06, 2.10e-9, 0.95])

step_n = 5                                  # Number of steps to take

p_init = newton_res                         # Initial Parameters

step_size = np.asarray([1.51210183e-05, 2.52357154e-06, 8.80062373e-07, 1.95227730e-06, 4.06549121e-10, 5.96586242e-07]) #np.sqrt(np.diagonal(cov_mat))   # Size of steps

chain = np.zeros([step_n, len(p_init)+1])   # Initializing the chain
chain[0, 0:-1] = p_init

cur_pos = chain[0, 0:-1]

pred = get_spectrum(pars)[:data_l]
cur_chi = chi2(pred, spec, errs)

chain[0, -1] = cur_chi

#test = step_size*np.random.randn(len(step_size))
#print(test)
#print(chi2(chain[0, 0:-1]+test, t, d))

file_name = str(datetime.datetime.now())

#assert(0==1)
for i in tqdm(range(1, step_n)):

    # Finding new position
    new_pos = chain[i-1, 0:-1] + step_size*np.random.randn(len(step_size))

    # Finding new chi square
    pred = get_spectrum(pars)[:data_l]
    new_chi = chi2(pred, spec, errs)

    #print("The new Chi^2 is:", new_chi)
    #print("The old Chi^2 is:", cur_chi)

    if new_chi < cur_chi:
        #print("Then we accept the change.")
        accept = True
    else:
        delt = new_chi - cur_chi
        #print("Then we don't immediately accept the change. Their difference is", delt)
        
        prob = np.exp(-0.5*delt)
        if np.random.rand() < prob:
            accept = True
        else:
            accept = False
    if accept:
        cur_chi = new_chi
        cur_pos = new_pos

    chain[i, 0:-1] = cur_pos
    chain[i, -1] = cur_chi

np.savetxt(file_name+".txt", chain)

# Printing values
best_param = np.mean(chain[10:-1, :], axis = 0)
param_std = np.std(chain[10:-1, :], axis = 0)

print("The best fit parameters are:", best_param[0:-1])
print("Their error is:", param_std[0:-1])

# Plotting 
plt.plot(chain[:, 0:-1], label = ["a", "t_0", "w", "b", "c", "dt"])
plt.title("Parameters Trajectory in MCMC")
plt.legend()
plt.xlabel("MCMC Steps")
plt.ylabel("Parameter Value")
plt.savefig("figs/a5q2_mcmc_parameters.jpg")
plt.show()
plt.clf()
