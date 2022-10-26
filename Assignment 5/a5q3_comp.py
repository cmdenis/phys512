import numpy as np
import os
import corner
import matplotlib.pyplot as plt
import random

file_list = os.listdir("runs")
file_list.sort()
file_list.remove(".DS_Store")


comp_dat = []


for file in file_list:
    comp_dat += list(np.loadtxt("runs/" + file))


comp_dat = np.array(comp_dat)

np.savetxt("compiled_runs.txt", comp_dat)

# Printing values
best_param = np.mean(comp_dat[10:-1, :], axis = 0)
param_std = np.std(comp_dat[10:-1, :], axis = 0)

print("The best fit parameters are:", best_param[0:-1])
print("Their error is:", param_std[0:-1])


plt.plot(comp_dat[:, 1])
#plt.plot(np.fft.irfft(np.fft.fft(comp_dat[1000:, 1])**2))
#plt.plot(np.fft.irfft(abs(np.fft.fft(nse))**2))
plt.show()


corner.corner(comp_dat[:, 0:-1])
plt.savefig("figs/a5q3_mcmc_corner_compiled.jpg")
#plt.show()