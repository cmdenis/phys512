import numpy as np 
from scipy import integrate
from matplotlib import pyplot as plt


# PART 1) Simulate the decay of U238

# Setting up our ode
def fun(x, y, half_life = [1, 0.1]):
    # This is for decay of U238
    dydx = np.zeros(len(half_life) + 1)
    dydx[0] = -y[0]/half_life[0]
    dydx[1] = y[0]/half_life[0] - y[1]/half_life[1]
    dydx[2] = y[1]/half_life[1]
    return dydx


# units for conversion to seconds
hours = 3.6e3
bil_years = 3.15e19
day = 8.64e4
year = 3.15e7
minute = 60



comp_name = [
    'U238',
    'Th234', 
    'Pr234', 
    'U234', 
    'Th230', 
    'Ra226',
    'Rn222',
    'Po218',
    'Pb214',
    'Bi214',
    'Po214',
    'Pb210',
    'Bi210',
    'Po210',
    'Pb206'
    ]

def fun(
    x, 
    y, 
    half_life = [
        4.468*bil_years,    # U238
        24.1*day,           # Th234
        6.70*hours,         # Pr234
        2.455e5*year,       # U234
        7.538e4 * year,     # Th230
        1.6e3 * year,       # Ra226
        3.8235 * day,       # Rn222
        3.1 * minute,       # Po218
        26.8 * minute,      # Pb214
        19.9 * minute,      # Bi214
        164.3 * 1e-6,       # Po214
        22.3 * year,        # Pb210
        5.015 * year,       # Bi210
        138376 * day        # Po210
    ]
    ):


    # This is for decay of U238

    half_life = np.array(half_life)/year # In the end I ended up switching to years
    dydx = np.zeros(len(half_life) + 1) # Initializing the array for the derivatives


    dydx[0] = -y[0]/half_life[0]        # Setting up the first transition

    for i in range(1, len(dydx) - 1):   # Going through a loop to create the vector
        dydx[i] = y[i-1]/half_life[i-1] - y[i]/half_life[i]
    
    dydx[-1] = y[-2]/half_life[-1]      # Setting up the last component

    return dydx



# Init simulation
y0 = np.insert(np.zeros(len(comp_name)-1), 0, 1)
x0 = 0
x1 = 1e12

sol = integrate.solve_ivp(fun, [x0, x1], y0, method = 'Radau')


for y_dat, name in zip(sol.y, comp_name):
    plt.loglog(sol.t, y_dat, label = name)
    
plt.title("U238 Radioactive decay")
plt.legend(prop={'size': 8})
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.savefig("figs/a3q2_decay_plot.jpg")
plt.show()




# PART 2) Plotting the ratio of two elements

plt.loglog(sol.t, sol.y[-1]/sol.y[0], label = "Ratio of Pb206/U238")

plt.title("Pb206/U238 ratio")
plt.legend()
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.savefig("figs/a3q2_decay_ratio_plot.jpg")
plt.show()


# Defining the analytic estimation for the ratio of U238 and Pb206
def rat_Pb_U(t):
    U_hl = 4.468*bil_years/year
    return np.exp(t/U_hl) - 1

v = np.linspace(x0, x1, 1000)
ana_y = rat_Pb_U(v)

plt.loglog(sol.t, sol.y[-1]/sol.y[0], label = "Numerical")
plt.loglog(v, ana_y, label = "Analytic")

plt.title("Pb206/U238 ratio")
plt.legend()
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.savefig("figs/a3q2_analytic_comp.jpg")
plt.show()


# PART 3) Plotting the ratio of U238 and Th230

plt.loglog(sol.t, sol.y[4]/sol.y[3], label = "Ratio of Th230/U238")

plt.title("Th230/U238 ratio")
plt.legend()
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.savefig("figs/a3q2_decay_ratio_plot2.jpg")
plt.show()

'''
# Defining the analytic estimation for the ratio of the two elements
def rat_Pb_U(t):
    U_hl = 4.468*bil_years/year
    return np.exp(t/U_hl) - 1

v = np.linspace(x0, x1, 1000)
ana_y = rat_Pb_U(v)

plt.loglog(sol.t, sol.y[-1]/sol.y[0], label = "Numerical")
plt.loglog(v, ana_y, label = "Analytic")

plt.title("Pb206/U238 ratio")
plt.legend()
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.savefig("figs/a3q2_analytic_comp.jpg")
plt.show()
'''
