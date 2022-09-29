import numpy as np 
from scipy import integrate
from matplotlib import pyplot as plt

def fun(x, y, half_life = [1, 0.1]):
    # This is for decay of U238
    dydx = np.zeros(len(half_life) + 1)
    dydx[0] = -y[0]/half_life[0]
    dydx[1] = y[0]/half_life[0] - y[1]/half_life[1]
    dydx[2] = y[1]/half_life[1]
    return dydx

hours = 3.6e3
bil_years = 3.15e19
day = 8.64e4
year = 3.15e7
minute = 60

half_lives = [
    4.468*bil_years,    # U238
    24*day,             # Th234
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

def fun(x, y, half_life = [4468, , 0.1]):
    # This is for decay of U238
    dydx = np.zeros(len(half_life) + 1)
    dydx[0] = -y[0]/half_life[0]
    dydx[1] = y[0]/half_life[0] - y[1]/half_life[1]
    dydx[2] = y[1]/half_life[1]
    return dydx

y0 = np.asarray([1, 0, 0])
x0 = 0
x1 = 1

ans_rk4 = integrate.solve_ivp(fun, [x0, x1], y0, method = 'RK45')
ans_stiff = integrate.solve_ivp(fun, [x0, x1], y0, method = 'Radau')
#print(ans_rk4.nfev, 'evaluations with rk4')
#print(ans_stiff.nfev, 'evaluations when implicitly')
#print('fincal values are:', ans_rk4.y[0, -1], ' and ', ans_stiff.y[0, -1])


for y_dat, labl in zip(ans_rk4.y, [1, 2, 3]):
    plt.plot(ans_rk4.t, y_dat, label = labl )
    
plt.title("2-State radioactive decay")
plt.legend()   
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.savefig("figs/a3q2_decay_plot.jpg")
plt.show()
