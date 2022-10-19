from use_func import get_spectrum
import numpy as np 
from scipy import integrate
from matplotlib import pyplot as plt


# Loading the data
planck = np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt', skiprows=1)
#planck = np.loadtxt('COM_PowerSpect_CMB-TT-binned_R3.01.txt', skiprows=1)

# Storing info from raw data
x_data = planck[:,0]
data_l = len(x_data)
y_data = planck[:,1]
errs = 0.5 * (planck[:,2] + planck[:,3])

# Guess parameters:
pars = np.asarray([69, 0.025, 0.12, 0.06, 2.10e-9, 0.95])

def num_derivs(fun,pars,dp):
    A=np.empty([data_l, len(pars)])
    for i in range(len(pars)):
        pp=pars.copy()
        pp[i]=pars[i]+dp[i]
        y_right=fun(pp)
        pp[i]=pars[i]#-dp[i]
        y_left=fun(pp)
        #A[:,i]=((y_right-y_left)/(2*dp[i]))[0:data_l]
        A[:,i]=((y_right-y_left)/(dp[i]))[0:data_l]
    return A

def num_newton(fun,pars,dp,x,y,sigma=1,niter=5):
    chi2prev = 1E12
    inv_N = np.eye(len(sigma))/sigma**2

    for i in range(niter):
        pred=fun(pars)[0:data_l]
        r=y-pred
        A=num_derivs(fun,pars,dp)
        lhs=A.T@inv_N@A
        rhs=A.T@inv_N@r
        u, s, vh = np.linalg.svd(lhs)

        inv_mat = vh.T@np.linalg.inv(np.diag(s))@u.T
        step=inv_mat@rhs
        #step=np.linalg.pinv(lhs)@rhs
        pars=pars+step
        chi2 =  np.sum(r**2/sigma**2)
        print(f"\nChi^2 = {chi2}")
        print(f"Chi^2 Diff = {chi2prev - chi2}")
        chi2prev = chi2
        print("Parameters are:", pars)
        print('Step is:', step)
    print("Done Newton Method\n")
    return pars,np.linalg.inv(lhs)


def chisq(y, pred, err):
    r=y-pred
    return np.sum(r**2/err**2)


fun = get_spectrum
p0 = np.asarray([69, 0.022, 0.12, 0.06, 2.10e-9, 0.95])
pars, cov_mat = num_newton(fun, p0, pars*1e-8, x_data, y_data, errs, 3)

pred = get_spectrum(pars)[0:data_l]
og_chi = chisq(y_data, pred, errs)

print("We found chi^2:", og_chi)

# Generating variations of samples
var_p = []
for i in range(5):
    var_p.append(pars + np.random.multivariate_normal(np.zeros(len(np.diagonal(cov_mat))), cov_mat))

sample_diff_chi = []
for i in range(len(var_p)):
    pred = get_spectrum(var_p[i])[0:data_l]
    sample_diff_chi.append(chisq(y_data, pred, errs) - og_chi)

print(pars)

print("Mean difference in Chi^2:", np.mean(sample_diff_chi))
print("Variations:", sample_diff_chi)

plt.plot(x_data, y_data, label = "Data")
plt.plot(x_data, pred, label = "Guess Fit")
plt.title("Data and fit")
plt.show()