import numpy as np
import matplotlib.pyplot as plt
import scipy as sci

def ndiff(fun, x, full = False):
    '''Function to take a derivate. Optional: can output the estimated error on the result.'''
    
    ini_step = 10**-5 # Initial step size (for third derivative)
    
    # Anonymous function to take derivatives
    diff_op = lambda func, x: (func(x + ini_step) - func(x - ini_step))/(2*ini_step)

    deriv_1 = lambda x: diff_op(fun, x)     # Calculating the first derivative
    deriv_2 = lambda x: diff_op(deriv_1, x) # Calculating the second derivative
    deriv_3 = lambda x: diff_op(deriv_1, x) # Calculating the thirs derivative
    
    # Optimal step size
    third_deriv = deriv_3(x)
    opt_step = (abs(3*np.finfo(float).eps/third_deriv))**(1/3)
    
    # Finding the derivative
    deriv = (fun(x + opt_step) - fun(x - opt_step))/(2*opt_step)

    # Conditional system for optional argument
    if full:
        # Returns derivative and estimated error
        est_err = np.finfo(float).eps/opt_step + third_deriv*opt_step**2/6
        return np.array([deriv, est_err])
    
    elif not full:
        # Returns only derivative
        return deriv

print(ndiff(np.exp, 1, full = True))