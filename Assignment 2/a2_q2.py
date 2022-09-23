import numpy as np
import matplotlib.pyplot as plt
import scipy as sci


def offset_gauss(x):
    return 1+10*np.exp(-0.5*x**2/(0.1)**2)


call_count_norm = 0
def integrate(fun,a,b,tol):
    global call_count_norm 
    call_count_norm += 5

    x = np.linspace(a, b, 5)
    dx = x[1] - x[0]
    y = fun(x)
    # 3-points integral
    i1 = (y[0]+4*y[2]+y[4])/3*(2*dx)
    i2 = (y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/3*dx
    myerr = np.abs(i1-i2)
    if myerr<tol:
        return i2
    else:
        mid=(a+b)/2
        int1=integrate(fun,a,mid,tol/2)
        int2=integrate(fun,mid,b,tol/2)
        return int1+int2


call_count_adapt = 0

def integrate_adaptive(fun, a, b, tol, extra = None, call_count = 0):
    
    global call_count_adapt

    sub_div = 5
    if extra == None:
        print('Starting Adaptive Integrator')
        x = np.linspace(a, b, sub_div)
        dx = x[1] - x[0]
        call_count_adapt += len(x) # To count the number of 
        y = fun(x)
    else:
        x = np.linspace(a, b, sub_div)
        dx = x[1] - x[0]
        call_count_adapt += 2
        y = np.array([extra[0], fun(x[1]), extra[1], fun(x[3]), extra[2]])

    # 3-points integral
    i1 = (y[0]+4*y[2]+y[4])/3*(2*dx)
    i2 = (y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/3*dx
    myerr = np.abs(i1-i2)
    if myerr < tol:
        return i2
    else:
        mid = (a+b)/2
        int1 = integrate_adaptive(fun, a, mid, tol/2, extra = [y[0], y[1], y[2]])
        int2 = integrate_adaptive(fun, mid, b, tol/2, extra = [y[2], y[3], y[4]])
        return int1+int2


# Gaussian test
# Testing integrate
print("Starting Normal Integrator:")
ans = integrate(offset_gauss, -5, 5, 1e-7)
print("Counts are:", call_count_norm)
print('Normal Integrator:', ans, ans-(10+np.sqrt(2*np.pi)))

# Testing integrate_adaptive
ans = integrate_adaptive(offset_gauss, -5, 5, 1e-7)
print("Counts are:", call_count_adapt)
print('Adaptive Integrator:', ans, ans-(10+np.sqrt(2*np.pi)))

print("The ratio of counts is:", call_count_adapt/call_count_norm)


# Exponential test
# Testing integrate
print("Starting Normal Integrator:")
ans = integrate(np.exp, 0, 2, 1e-7)
print("Counts are:", call_count_norm)
print('Normal Integrator:', ans, ans-(np.exp(2)-np.exp(0)))

# Testing integrate_adaptive
ans = integrate_adaptive(np.exp, 0, 2, 1e-7)
print("Counts are:", call_count_adapt)
print('Adaptive Integrator:', ans, ans-(np.exp(2)-np.exp(0)))

print("The ratio of counts is:", call_count_adapt/call_count_norm)


# Sine test
# Testing integrate
print("Starting Normal Integrator:")
ans = integrate(np.sin, -5, 5, 1e-7)
print("Counts are:", call_count_norm)
print('Normal Integrator:', ans, ans-(np.cos(5)-np.cos(-5)))

# Testing integrate_adaptive
ans = integrate_adaptive(np.sin, -5, 5, 1e-7)
print("Counts are:", call_count_adapt)
print('Adaptive Integrator:', ans, ans-(np.cos(5)-np.cos(-5)))

print("The ratio of counts is:", call_count_adapt/call_count_norm)

