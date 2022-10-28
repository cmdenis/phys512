import numpy as np
import os
import corner
import matplotlib.pyplot as plt
import random
from use_func import *

file_list = os.listdir("runs4")
file_list.sort()
#file_list.remove(".DS_Store")
file_list.remove("compiled_runs.txt")


comp_dat = []


for file in file_list:
    comp_dat += list(np.loadtxt("runs4/" + file))


comp_dat = np.array(comp_dat)

np.savetxt("runs4/compiled_runs.txt", comp_dat)

# Printing values
best_param = np.mean(comp_dat[:, :], axis = 0)
param_std = np.std(comp_dat[:, :], axis = 0)

print("The best fit parameters are:", best_param[0:-1])
print("Their error is:", param_std[0:-1])

# Calculating optimized Chi^2
def chisq(y, pred, err):
    r=y-pred
    return np.sum(r**2/err**2)

planck = np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt', skiprows=1) # Importing raw data
x_data = planck[:,0]    # X axis data
y_data = planck[:,1]    # Y axis data
errs = 0.5 * (planck[:,2] + planck[:,3])    # Errors on data
data_l = len(x_data)    # Length of data
pred = get_spectrum(best_param)[:data_l]
cur_chi = chisq(pred, y_data, errs)
print("Optimized Chi^2 is:", cur_chi)


for i in range(len(best_param)):
    plt.plot(comp_dat[:, i])
    plt.title("Parameters over steps")
    plt.xlabel("Steps")
    plt.ylabel("Parameter value")
    plt.savefig("figs/q4param"+str(i)+"run.jpg")
    plt.clf()

    plt.clf()
    plt.loglog(np.abs(np.fft.rfft(comp_dat[:, i]))**2)
    plt.title("Power Spectrum")
    plt.xlabel("Frequency Space")
    plt.ylabel("Intensity")
    plt.savefig("figs/q4param"+str(i)+"power.jpg")
    plt.clf()


plt.clf()
corner.corner(comp_dat[:, 0:-1])
plt.savefig("figs/a5q3_mcmc_corner_compiled.jpg")
plt.show()

p_const = 6.62e-34

de_mean_val = 1 -  best_param[2]/p_const**2 - best_param[3]/p_const**2

print("The mean value of dark energy is", de_mean_val)
print("With uncertainty:", param_std[2]/p_const**2 + param_std[3]/p_const**2)