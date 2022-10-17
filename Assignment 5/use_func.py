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