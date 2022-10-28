from dataclasses import replace
import numpy as np 
from scipy import integrate
from matplotlib import pyplot as plt
import camb
from tqdm import tqdm
from use_func import *
import time
import datetime
import corner
import os



# Loading the data
planck = np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt', skiprows=1) # Importing raw data
x_data = planck[:,0]    # X axis data
y_data = planck[:,1]    # Y axis data
errs = 0.5 * (planck[:,2] + planck[:,3])    # Errors on data

data_l = len(x_data)    # Length of data



def chisq(y, pred, err):
    r=y-pred
    return np.sum(r**2/err**2)

file_list = os.listdir("runs3")
file_list.sort()
prev_file = file_list[-1]


#assert(0==1)

print("Loading Newton Method Covariance Matrix and last step of simulation...\n")
pars = np.loadtxt("runs3/"+prev_file)[-1, 0:-1]
print(pars)
cov_mat = np.loadtxt("planck_fit_cov.txt")

# Storing info from raw data
ell = planck[:,0]
data_l = len(ell)
spec = planck[:,1]
errs = 0.5 * (planck[:,2] + planck[:,3])

pred = get_spectrum(pars)[:data_l]
cur_chi = chisq(pred, spec, errs)

print(cur_chi)


print("Starting Params:", pars, "with chi^2:", cur_chi)
print("Covariance Matrix:", cov_mat)

print("\n Entering MCMC hyperspace. Wavefunction collapsing in progress... \n")



step_n = 16000                                  # Number of steps to take

p_init = pars                         # Initial Parameters


chain = np.zeros([step_n, len(p_init)+1])   # Initializing the chain
chain[0, 0:-1] = p_init

cur_pos = chain[0, 0:-1]

#pred = get_spectrum(pars)[:data_l]
#cur_chi = chi2(pred, spec, errs)

chain[0, -1] = cur_chi

#test = step_size*np.random.randn(len(step_size))
#print(test)
#print(chi2(chain[0, 0:-1]+test, t, d))

file_name = str(datetime.datetime.now()).replace("/", "-").replace(":", "-")

#assert(0==1)
for i in tqdm(range(1, step_n)):

    # Finding new position
    new_pos = chain[i-1, 0:-1] + np.random.multivariate_normal(np.zeros(len(np.diagonal(cov_mat))), cov_mat)*0.3

    # Finding new chi square
    pred = get_spectrum(new_pos)[:data_l]
    new_chi = chisq(pred, spec, errs)

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
            #print("Ended up going through...")
            accept = True
        else:
            #print("Ended up NOT going through...")
            accept = False
    if accept:
        cur_chi = new_chi
        cur_pos = new_pos

    chain[i, 0:-1] = cur_pos
    chain[i, -1] = cur_chi

    #if np.mod(i, 100) == 0 :
    #    corner.corner(chain[:, 0:-1])
    #    plt.savefig("figs/"+file_name+"a5q2_mcmc_corner.jpg")

    np.savetxt("runs3/"+file_name+".txt", chain)

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
#plt.show()
plt.clf()

corner.corner(chain[:, 0:-1])
plt.savefig("figs/a5q2_mcmc_corner.jpg")
