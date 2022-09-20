import numpy as np
import matplotlib.pyplot as plt
import scipy as sci


def offset_gauss(x):
    return 1+10*np.exp(-0.5*x**2/(0.1)**2)

def integrate(fun,a,b,tol):
    print('calling function from ', a, b)
    x = np.linspace(a, b, 5)
    dx = x[1] - x[0]
    y = fun(x)
    #do the 3-point integral
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

def integrate_adaptive(fun, a, b, tol, extra = None):
    print('calling adaptive function from ', a, b)
    nb_points = 5
    if extra == None:
        x = np.linspace(a, b, nb_points)
        dx = x[1] - x[0]
        y = fun(x)
    else:
        x = np.linspace(a, b, nb_points)
        dx = x[1] - x[0]
        y = np.array([extra[0], fun(x[1]), extra[1], fun(x[3]), extra[2]])

    #do the 3-point integral
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

# To answer the question, we should basically pass the x vector in the extra argument and then we will only need to compute two of the 5 values of the vector

ans = integrate(offset_gauss, -4, 6, 1e-6)
ans2 = integrate(offset_gauss, -4, 0, 1e-6) + integrate(offset_gauss, 0, 6,1e-6)
print('Normal Integrator:', ans, ans2, ans-(10+np.sqrt(2*np.pi)))


ans = integrate_adaptive(offset_gauss, -4, 6, 1e-6)
ans2 = integrate_adaptive(offset_gauss, -4, 0, 1e-6) + integrate_adaptive(offset_gauss, 0, 6,1e-6)
print('Adaptive Integrator:', ans, ans2, ans-(10+np.sqrt(2*np.pi)))
