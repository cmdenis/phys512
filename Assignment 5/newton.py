from use_func import get_spectrum
import numpy as np 
from scipy import integrate
from matplotlib import pyplot as plt




# Loading the data
planck = np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt', skiprows=1)

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
        pp[i]=pars[i]-dp[i]
        y_left=fun(pp)
        A[:,i]=((y_right-y_left)/(2*dp[i]))[0:data_l]
    return A

def num_newton(fun,pars,dp,x,y,sigma=1,niter=5):
    chi2prev = 1E12
    for i in range(niter):
        pred=fun(pars)[0:data_l]
        r=y-pred
        A=num_derivs(fun,pars,dp)
        sigma = np.linalg(sigma)**2
        lhs=A.T@sigma@A
        rhs=A.T@sigma@r
        step=np.linalg.svd(lhs)@rhs
        pars=pars+step
        chi2 =  np.sum(r**2@sigma**2)
        print(f"chi2 = {chi2}")
        print(f"chi2diff = {chi2prev - chi2}")
        chi2prev = chi2
        #print('step is ',step)
    return pars,np.linalg.inv(lhs)






fun = get_spectrum
p0 = np.asarray([69, 0.022, 0.12, 0.06, 2.10e-9, 0.95])
p0 = np.asarray([60,0.02,0.1,0.05,2.00e-9,1.0])
#p0 = np.array([6.84081297e+01, 2.28896704e-02, 1.19637761e-01, -5.63652856e-02, 1.68178216e-09, 9.78754015e-01])

plt.plot(x_data, get_spectrum(p0)[0:data_l], label = "Guess Fit")
plt.plot(x_data, y_data, label = "Data")
plt.title("Data and fit")
#plt.show()



triple_num_pars,triple_num_cov = num_newton(fun, p0, np.repeat(1e-10, len(p0)), x_data, y_data, 1, 5)

print(triple_num_pars)
