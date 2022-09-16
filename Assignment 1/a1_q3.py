import numpy as np
import matplotlib.pyplot as plt
import scipy as sci



dat = np.loadtxt("lakeshore.txt") # Importing raw data


# Taking a look at the data
temperatures = np.array([i[0] for i in dat[::-1]]) # Making array for temperatures (from raw data)
voltages = np.array([i[1] for i in dat[::-1]]) # Making array for volatages (from raw data)
slopes = np.array([i[2] for i in dat])

print(min(voltages))

# Plotting the data
plt.scatter(voltages, temperatures)
plt.title("Temperature vs Voltages for Lakeshore 670 diode")
plt.ylabel("Temperature (˚C)")
plt.xlabel("Voltage (V)")
plt.show()




# Plotting the cubic spline to the data
cs = sci.interpolate.CubicSpline(voltages, temperatures)

# Finding a linear spline
cs1 = sci.interpolate.interp1d(voltages, temperatures)

# Plotting the data
v = np.linspace(min(voltages), max(voltages), 1000)
plt.scatter(voltages, temperatures, label = "Data Points")
plt.plot(v, cs(v), label = "Cubic Spline", color = 'r')
plt.title("Temperature vs Voltages for Lakeshore 670 diode")
plt.ylabel("Temperature (˚C)")
plt.xlabel("Voltage (V)")
plt.legend()
plt.savefig("figs/q3_interp_plot.jpg")
plt.show()

# Error between linear and cubic splines
print("The error is about:", np.std(cs(v) - cs1(v)))










# Defining lakeshore function

def lakeshore(V, data):
    '''Function to interpolate the data with a cubic spline'''
    temperatures = np.array([i[0] for i in data[::-1]]) # Making array for temperatures (from raw data)
    voltages = np.array([i[1] for i in data[::-1]]) # Making array for volatages (from raw data)
    approx_vol_size = abs(voltages[2]-voltages[1])
    
    # First we find the interpolation
    cs = sci.interpolate.CubicSpline(voltages, temperatures) # spline function
    inter_val = cs(V)

    # Now we roughly estimate the error
    lin_spline = sci.interpolate.interp1d(voltages, temperatures) # Linear interpolation
    err_range = np.array(np.linspace(V-approx_vol_size, V+approx_vol_size, 1000))

    approx_error = np.std(abs(lin_spline(err_range) - cs(err_range)), axis = 0)



    return inter_val, approx_error

print(lakeshore(np.array([1.11, 0.2, 0.3]), dat))