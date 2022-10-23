import numpy as np
import camb

def chi2(pred, data, err):
    return np.sum(((pred - data)/err)**2)


def get_spectrum(pars, lmax = 3000):
    H0 = pars[0]    # Hubble constant
    ombh2 = pars[1] # Baryon density
    omch2 = pars[2] # Dark matter density
    tau = pars[3]   # Optical depth
    As = pars[4]    # Primordial amp of spectrum
    ns = pars[5]    # Primordial tilt of spectrum
    pars = camb.CAMBparams()
    pars.set_cosmology(H0 = H0, ombh2=ombh2, omch2=omch2, mnu=0.06, omk=0, tau=tau)
    pars.InitPower.set_params(As=As, ns=ns, r=0)
    pars.set_for_lmax(lmax, lens_potential_accuracy=0)
    results = camb.get_results(pars)
    powers = results.get_cmb_power_spectra(pars, CMB_unit='muK')
    cmb = powers['total']
    tt = cmb[:, 0]    #you could return the full power spectrum here if you wanted to do say EE
    return tt[2:]

# First we create a derivative taking function
def p_deriv(func, p_ind, p):
    shift = 1e-12
    # Creating array
    dp = np.zeros(len(p))
    # Check if derivative index is an integer
    if isinstance(p_ind, int):
        dp[p_ind] = shift 
    else:
        raise ValueError("Derivative index must be an integer.")
    return (func(p + dp) - func(p - dp))/(2*shift)    # Two sided derivative

# Define numerical gradient taker
def grad_f(f, p, length = None):
    # Initializing the gradient
    grad = np.zeros([2507, p.size])
    # Finding the derivatives to build the gradient
    if length == None:
        for param in range(len(p)):
            grad[:, param] = p_deriv(f, param, p)
    else:
        for param in range(len(p)):
            grad[:, param] = p_deriv(f, param, p)[:length]
    return grad