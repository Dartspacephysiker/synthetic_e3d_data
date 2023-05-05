import numpy as np

def func_gammaT_div_gamma0(gamma1_div_gamma0,integT,fwhmRange=3,dutyCycle=0.25):
    """
    Calculate relative noise γ_T / γ_0. See Slack chat with Ilkka in ISSI - Understanding mesoscale ionospheric electrodynamics using regional data assimilation

    gamma1_div_gamma0 : Relative noise for one bit in the modulation cycle (i.e., "a ridiculously short time period," in Ilkka's words). I use gamma0 = 0.01 in get_noise_estimates.R
    integT            : Integration time in seconds
    fwhmRange         : Range ambiguity [km]
    dutyCycle         : Transmitter duty cycle
    """
    return gamma1_div_gamma0 * np.sqrt( fwhmRange/ (1.5e5 * integT * dutyCycle) )

def vError(gammaT_div_gamma0, errVi_div_Vi0, Vi0=1):
    """
    gammaT_div_gamma0 : Relative noise for integration time T.          [s  ] Output from func_gammaT_div_gamma0.
    errVi_div_Vi0     : Relative error, err(Vi/Vi0), from lookup table. [   ]
    Vi0               : Ilkka says that they use Vi0 = 1 m/s.           [m/s]
    """
    return gammaT_div_gamma0 * errVi_div_Vi0 * Vi0


def isoError(gammaT_div_gamma0, paramGEMINI, erriso_div_iso0):
    return gammaT_div_gamma0 * erriso_div_iso0 * paramGEMINI

