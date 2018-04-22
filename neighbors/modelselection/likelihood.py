# Third-party
import numpy as np
from scipy.linalg import block_diag
from scipy.special import logsumexp

# Project
from .linalg import get_y_Cinv, get_M, get_Ainv_nu_Delta

__all__ = ['FML_helper', 'ln_marg_v_likelihood']


def FML_helper(d, data, Cov, Vinv, v_scatter):
    """

    Parameters
    ----------
    d : `numpy.ndarray`
        Should have shape ``(nstars,)`` distance for each star.
    data : `numpy.ndarray`
        Should have shape ``(nstars, 6)``:
        [ra, dec, parallax, pmra, pmdec, radial_velocity] for each star.
        Units should be: [radian, radian, mas, mas/yr, mas/yr, km/s].
    Cov : `numpy.ndarray`
        Should have shape ``(nstars, 6, 6)``: full covariance matrix for
        each star. Units should match units of ``data``.
    Vinv : `numpy.ndarray`
        Should have shape ``(3, 3)`` and units [1/(km**2/s**2)]: this is the
        inverse variance matrix for one component of the mixture-of-Gaussians
        prior for the velocity distribution.
    v_scatter : numeric [km/s]
        Velocity tolerance in what we call "comoving."

    """
    nstars = len(d)

    y_Cinvs = [get_y_Cinv(d[i], data[i, 3:], Cov[i, 3:, 3:], v_scatter)
               for i in range(nstars)]
    y = np.concatenate([x[0] for x in y_Cinvs])
    Cinv = block_diag(*[x[1] for x in y_Cinvs])

    M = get_M(data[:, 0], data[:, 1])

    Ainv, nu, Delta = get_Ainv_nu_Delta(d, M, Cinv, y, Vinv)
    sgn, log_detAinv = np.linalg.slogdet(Ainv/(2*np.pi))
    log_detA = -log_detAinv

    assert sgn > 0
    return 0.5*log_detA - Delta


def ln_marg_v_likelihood(d, data, Cov, prior, v_scatter=0.,):
    """
    Marginalized over true v.

    Parameters
    ----------
    d : `numpy.ndarray`
        Should have shape ``(nstars,)`` distance for each star.
    data : `numpy.ndarray`
        Should have shape ``(nstars, 6)``:
        [ra, dec, parallax, pmra, pmdec, radial_velocity] for each star.
        Units should be: [radian, radian, mas, mas/yr, mas/yr, km/s].
    Cov : `numpy.ndarray`
        Should have shape ``(nstars, 6, 6)``: full covariance matrix for
        each star. Units should match units of ``data``.
    prior : dict
        For now, this is just a dictionary that contains information needed for
        the velocity prior. This should contain keys: 'Vinvs' and 'weights'.
        These are the inverse variance matrix for the mixture-of-Gaussians prior
        for the velocity distribution, and the weights for each mixture
        component (these should sum to 1).
    v_scatter : numeric [km/s]
        Velocity tolerance in what we call "comoving."
    """
    return logsumexp([FML_helper(d, data, Cov, Vinv, v_scatter)
                      for Vinv in prior['Vinvs']],
                     b=prior['weights'])
