import numpy as np
import matplotlib.pyplot as plt
import scipy as sci



# First plots
# First we create our sampling points
fun = np.cos
nb_points = 5
x_data = np.linspace(-np.pi/2, np.pi/2, nb_points)
y_data = fun(x_data)
v = np.linspace(-np.pi/2, np.pi/2, 100)


# Polynomial fitting
pp = np.polyfit(x_data, y_data, nb_points-1)

# Cubic spline fitting
cs = sci.interpolate.CubicSpline(x_data, y_data)


# Rational function fitting
n = 2
m = nb_points - n - 1
pcols = [x_data**k for k in range(n+1)]
pmat = np.vstack(pcols)

qcols = [-x_data**k*y_data for k in range(1, m+1)]
qmat = np.vstack(qcols)
mat = np.hstack([pmat.T, qmat.T])
coeffs = np.linalg.inv(mat)@y_data

p=0 
for i in range(n+1):
    p = p + coeffs[i]*v**i
qq = 1
for i in range(m):
    qq = qq + coeffs[n+1+i]*v**(i+1)
pred = p/qq

# We plot the resulting curves
plt.scatter(x_data, fun(x_data), label = "Sample Points")
plt.plot(v, fun(v), label = "Analytical Result")
plt.plot(v, np.polyval(pp, v), label = "Polynomial")
plt.plot(v, cs(v), label = "Cubic Spline")
plt.plot(v, pred, label = "Rational Funtion")
plt.title("Interpolations")
plt.legend()
plt.savefig("figs/q4_interp1.jpg")
plt.show()

# Now the residuals
plt.scatter(x_data, fun(x_data)*0, label = "Sample Points")
plt.plot(v, fun(v)*0, label = "Analytical Result")
plt.plot(v, np.polyval(pp, v)-fun(v), label = "Polynomial")
plt.plot(v, cs(v)-fun(v), label = "Cubic Spline")
plt.plot(v, pred-fun(v), label = "Rational Funtion")
plt.title("Residuals for all interpolations")
plt.legend()
plt.savefig("figs/q4_resi1.jpg")
plt.show()

# Getting the error for each interpolation
print("The error for the polynomial interpolation is", np.std(np.polyval(pp, v)-fun(v)))
print("The error for the cubic spline interpolation is", np.std(cs(v)-fun(v)))
print("The error for the rational function interpolation is", np.std(pred-fun(v)))



# Second Plots

# First we create our sampling points
fun = lambda x : 1/(1+x**2)
nb_points = 5
x_data = np.linspace(-np.pi/2, np.pi/2, nb_points)
y_data = fun(x_data)
v = np.linspace(-np.pi/2, np.pi/2, 100)


# Polynomial fitting
pp = np.polyfit(x_data, y_data, nb_points-1)

# Cubic spline fitting
cs = sci.interpolate.CubicSpline(x_data, y_data)


# Rational function fitting
n = 2
m = nb_points - n - 1
pcols = [x_data**k for k in range(n+1)]
pmat = np.vstack(pcols)

qcols = [-x_data**k*y_data for k in range(1, m+1)]
qmat = np.vstack(qcols)
mat = np.hstack([pmat.T, qmat.T])
coeffs = np.linalg.inv(mat)@y_data

p=0 
for i in range(n+1):
    p = p + coeffs[i]*v**i
qq = 1
for i in range(m):
    qq = qq + coeffs[n+1+i]*v**(i+1)
pred = p/qq




# We plot the resulting curves
plt.scatter(x_data, fun(x_data), label = "Sample Points")
plt.plot(v, fun(v), label = "Analytical Result")
plt.plot(v, np.polyval(pp, v), label = "Polynomial")
plt.plot(v, cs(v), label = "Cubic Spline")
plt.plot(v, pred, label = "Rational Funtion")
plt.title("Interpolations")
plt.legend()
plt.savefig("figs/q4_interp2.jpg")
plt.show()

# Now the residuals
plt.scatter(x_data, fun(x_data)*0, label = "Sample Points")
plt.plot(v, fun(v)*0, label = "Analytical Result")
plt.plot(v, np.polyval(pp, v)-fun(v), label = "Polynomial")
plt.plot(v, cs(v)-fun(v), label = "Cubic Spline")
plt.plot(v, pred-fun(v), label = "Rational Funtion")
plt.title("Residuals for all interpolations")
plt.legend()
plt.savefig("figs/q4_resi2.jpg")
plt.show()

# Getting the error for each interpolation
print("The error for the polynomial interpolation is", np.std(np.polyval(pp, v)-fun(v)))
print("The error for the cubic spline interpolation is", np.std(cs(v)-fun(v)))
print("The error for the rational function interpolation is", np.std(pred-fun(v)))


# Third plots

# First we create our sampling points
fun = lambda x : 1/(1+x**2)
nb_points = 10
x_data = np.linspace(-np.pi/2, np.pi/2, nb_points)
y_data = fun(x_data)
v = np.linspace(-np.pi/2, np.pi/2, 100)


# Polynomial fitting
pp = np.polyfit(x_data, y_data, nb_points-1)

# Cubic spline fitting
cs = sci.interpolate.CubicSpline(x_data, y_data)


# Rational function fitting
n = 4
m = nb_points - n - 1
pcols = [x_data**k for k in range(n+1)]
pmat = np.vstack(pcols)

qcols = [-x_data**k*y_data for k in range(1, m+1)]
qmat = np.vstack(qcols)
mat = np.hstack([pmat.T, qmat.T])
coeffs = np.linalg.inv(mat)@y_data

p=0 
for i in range(n+1):
    p = p + coeffs[i]*v**i
qq = 1
for i in range(m):
    qq = qq + coeffs[n+1+i]*v**(i+1)
pred = p/qq




# We plot the resulting curves
plt.scatter(x_data, fun(x_data), label = "Sample Points")
plt.plot(v, fun(v), label = "Analytical Result")
plt.plot(v, np.polyval(pp, v), label = "Polynomial")
plt.plot(v, cs(v), label = "Cubic Spline")
plt.plot(v, pred, label = "Rational Funtion")
plt.title("Interpolations")
plt.legend()
plt.savefig("figs/q4_interp3.jpg")
plt.show()

# Now the residuals
plt.scatter(x_data, fun(x_data)*0, label = "Sample Points")
plt.plot(v, fun(v)*0, label = "Analytical Result")
plt.plot(v, np.polyval(pp, v)-fun(v), label = "Polynomial")
plt.plot(v, cs(v)-fun(v), label = "Cubic Spline")
plt.plot(v, pred-fun(v), label = "Rational Funtion")
plt.title("Residuals for all interpolations")
plt.legend()
plt.savefig("figs/q4_resi3.jpg")
plt.show()

# Getting the error for each interpolation
print("The error for the polynomial interpolation is", np.std(np.polyval(pp, v)-fun(v)))
print("The error for the cubic spline interpolation is", np.std(cs(v)-fun(v)))
print("The error for the rational function interpolation is", np.std(pred-fun(v)))







# fourth plots

# First we create our sampling points
fun = lambda x : 1/(1+x**2)
nb_points = 10
x_data = np.linspace(-np.pi/2, np.pi/2, nb_points)
y_data = fun(x_data)
v = np.linspace(-np.pi/2, np.pi/2, 100)


# Polynomial fitting
pp = np.polyfit(x_data, y_data, nb_points-1)

# Cubic spline fitting
cs = sci.interpolate.CubicSpline(x_data, y_data)


# Rational function fitting
n = 4
m = nb_points - n - 1
pcols = [x_data**k for k in range(n+1)]
pmat = np.vstack(pcols)

qcols = [-x_data**k*y_data for k in range(1, m+1)]
qmat = np.vstack(qcols)
mat = np.hstack([pmat.T, qmat.T])
coeffs = np.linalg.pinv(mat)@y_data

p=0 
for i in range(n+1):
    p = p + coeffs[i]*v**i
qq = 1
for i in range(m):
    qq = qq + coeffs[n+1+i]*v**(i+1)
pred = p/qq

# We plot the resulting curves
plt.scatter(x_data, fun(x_data), label = "Sample Points")
plt.plot(v, fun(v), label = "Analytical Result")
plt.plot(v, np.polyval(pp, v), label = "Polynomial")
plt.plot(v, cs(v), label = "Cubic Spline")
plt.plot(v, pred, label = "Rational Funtion")
plt.title("Interpolations")
plt.legend()
plt.savefig("figs/q4_interp4.jpg")
plt.show()

# Now the residuals
plt.scatter(x_data, fun(x_data)*0, label = "Sample Points")
plt.plot(v, fun(v)*0, label = "Analytical Result")
plt.plot(v, np.polyval(pp, v)-fun(v), label = "Polynomial")
plt.plot(v, cs(v)-fun(v), label = "Cubic Spline")
plt.plot(v, pred-fun(v), label = "Rational Funtion")
plt.title("Residuals for all interpolations")
plt.legend()
plt.savefig("figs/q4_resi4.jpg")
plt.show()

# Getting the error for each interpolation
print("The error for the polynomial interpolation is", np.std(np.polyval(pp, v)-fun(v)))
print("The error for the cubic spline interpolation is", np.std(cs(v)-fun(v)))
print("The error for the rational function interpolation is", np.std(pred-fun(v)))