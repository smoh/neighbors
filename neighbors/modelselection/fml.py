# Third-party
import numpy as np
from scipy.special import logsumexp
from scipy.integrate import simps

# Project
from .likelihood import ln_marg_v_likelihood

__all__ = ['ln_H1_FML', 'ln_H2_FML']


rlim = 2000. # pc; HACK: hard-coded limit. Deal With It!
A = np.log(3) - 3*np.log(rlim) - 0.5 * np.log(2 * np.pi)
def ln_prob(r, plx, err):
    """
    Evaluate log(likelihood * prior). This assumes a uniform space density
    prior.
    """
    B = -np.log(err)
    C = -0.5 * (1000 / r - plx)**2 / err**2
    D = 2 * np.log(r)
    return A + B + C + D


def get_mode(plx, err):
    """Mode of the distance posterior with a uniform space density prior."""
    snr = plx / err
    return 1000 * snr**2 / (4 * plx) * (1 - np.sqrt(1 - 8 / snr**2))


def get_rgrid(plx, err, ngrid=25):
    return np.linspace(1000/(plx + 4*err),
                       min(rlim, 1000/(plx - 4*err)),
                       ngrid)


def ln_H1_FML(data1, Cov1, data2, Cov2, Vinvs, prior_weights,
              n_dist_grid=32, v_scatter=0.):
    # TODO: need to construct a 2D d1, d2 grid here!
    ll_H1_at_samples = np.array([ln_H1_marg_v_likelihood(d1, d2, star1, star2, Vinv, v_scatter, prior_weights=prior_weights)
                                 for d1, d2 in zip(dist1, dist2)])
    return logsumexp(ll_H1_at_samples) - np.log(float(ll_H1_at_samples.size))


##############################################################################
#  Likelihood 2:
#
def ln_integrand2(r, data, Cov, Vinvs, prior_weights):
    ln_fr = ln_marg_v_likelihood([r], data, Cov, Vinvs,
                                 prior_weights=prior_weights)
    ln_p = ln_prob(r, data[0, 2], np.sqrt(Cov[0, 2, 2]), rlim)
    return ln_fr + ln_p


def ln_Q(data, Cov, Vinvs, prior_weights, n_dist_grid):
    rgrid = get_rgrid(data[0, 2], np.sqrt(Cov[0, 2, 2]), n_dist_grid)
    func_vals = [ln_integrand2(r, data, Cov, Vinvs, prior_weights)
                 for r in rgrid]
    return simps(np.exp(func_vals), x=rgrid)


def ln_H2_FML(star1, star2, Vinv, n_dist_samples=128, v_scatter=0.,
              prior_weights=None):
    dist1 = get_posterior_distance_samples(star1, size=n_dist_samples)
    dist2 = get_posterior_distance_samples(star2, size=n_dist_samples)
    ll_H2_at_samples = np.array([ln_H2_marg_v_likelihood(d1, d2, star1, star2, Vinv, v_scatter, prior_weights=prior_weights)
                                 for d1,d2 in zip(dist1, dist2)])
    return logsumexp(ll_H2_at_samples) - np.log(float(ll_H2_at_samples.size))
