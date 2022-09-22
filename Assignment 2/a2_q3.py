import numpy as np
import matplotlib.pyplot as plt
import scipy as sci



# Part a)

fun = np.log2 # Function we want to fit

n = 20 # Number of sample points
x_s = np.linspace(0.5, 1, n) # Linspace for sample points
y_s = fun(x_s)   # Sample points
x_new = (x_s-0.75)*4  # Shifting and stretching the range to -1 to 1

cheb_fit = np.polynomial.chebyshev.chebfit(x_new, y_s, n-1)[:9] # Coefficient of fit

x = np.linspace(0.5, 1, 1001) # Linspace for result
y = np.polynomial.chebyshev.Chebyshev(cheb_fit)((x-0.75)*4) # Image of result

print("The coefficients for the Chebyshev polynomials are:", cheb_fit)
print("The error on that range is:", np.std(y-fun(x)))

def cheb_log2(x):
    return np.polynomial.chebyshev.Chebyshev(
        [-4.56893394e-01, 4.95054673e-01, -4.24689768e-02, 4.85768297e-03, -6.25084976e-04, 8.57981013e-05, -1.22671891e-05, 1.80404306e-06, -2.70834251e-07]
    )((x - 0.75) * 4)




plt.plot(x, cheb_log2(x), label = "cheb_log2")
plt.plot(x, fun(x), label = "np.log2")
plt.scatter(x_s, y_s, label = "Sample Points")
plt.legend()
plt.savefig("figs/a2q2_logfunc_comparison.jpg")
plt.show()

plt.plot(x, cheb_log2(x)-fun(x), label = "cheb_log2")
plt.plot(x_s, y_s*0)
plt.legend()
plt.savefig("figs/a2q2_residuals.jpg")
plt.show()


# Part b)

def mylog2(x):
    mantissa, expo = np.frexp(x)
    return cheb_log2(mantissa) + expo


# Check if log works
x = np.linspace(0.001, 10, 1001)
y1 = mylog2(x)
y2 = np.log2(x)

plt.clf()
plt.plot(x, y1, label = "mylog2")
plt.plot(x, y2, label = "np.log2")
plt.legend()
plt.savefig("figs/a2q3_partb.jpg")
plt.show()

plt.clf()
plt.plot(x, y1 - y2, label = "mylog2")
plt.plot(x, y2 - y2, label = "np.log2")
plt.legend()
plt.savefig("figs/a2q3_partb_residuals.jpg")
plt.show()
