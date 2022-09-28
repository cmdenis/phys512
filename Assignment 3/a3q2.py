import numpy as np 
from scipy import integrate

def fun(x, y, half_life = [1, 1e-5]):
    # This is for a 2-state radio active decay
    dydx = np.zeros(len(half_life) + 1)
    dydx[0] = -y[0]/half_life[0]
    dydx[1] = y[0]/half_life[0] - y[1]/half_life[1]
    dydx[2] = y[1]/half_life[1]
    return dydx

y0 = np.asarray([1, 0, 0])
x0 = 0
x1 = 1

ans_rk4 = integrate.solve_ivp(fun, [x0, x1], y0)
ans_stiff = integrate.solve_ivp(fun, [x0, x1], y0, method = 'Radau')
print(ans_rk4.nfev, 'evaluations with rk4')
print(ans_stiff.nfev, 'evaluations when implicitly')
print('fincal values are:', ans_rk4.y[0, -1], ' and ', ans_stiff.y[0, -1])