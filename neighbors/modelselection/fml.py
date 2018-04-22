# Standard library
import math

# Third-party
import astropy.units as u
from astropy.utils.misc import isiterable
import numpy as np
from scipy.integrate import simps

# Project
from .data import data_to_vec_cov
from .likelihood import ln_marg_v_likelihood

__all__ = ['ComovingFML']


def get_mode(plx, err):
    """Mode of the distance posterior with a uniform space density prior."""
    snr = plx / err
    return 1000 * snr**2 / (4 * plx) * (1 - np.sqrt(1 - 8 / snr**2))


def get_rgrid(plx, err, rlim, ngrid=25):
    """
    Get a uniform grid of true distance to compute the fully marginalized
    likelihood for a star or set of stars.
    """
    return np.linspace(1000 / (plx + 4*err),
                       min(rlim, 1000 / (plx - 4*err)),
                       ngrid)


# First term in the log-likelihood + log-prior evaluation
_const = math.log(3) - 0.5 * math.log(2 * np.pi)

def ln_dist_prob(r, plx, err, rlim):
    """
    Evaluate log(likelihood * prior) for true distance, ``r``. This assumes
    a uniform space density prior.

    Parameters
    ----------
    TODO

    """
    A = _const - 3 * math.log(rlim)
    B = -math.log(err)
    C = -0.5 * ((1000 / r - plx) / err)**2
    D = 2 * math.log(r)
    return A + B + C + D


def ln_dist_integrand(rs, data, Cov, v_prior, rlim, v_scatter):
    """
    Integrand for the distance marginalization for one or two stars.

    Parameters
    ----------
    TODO

    """
    ln_f = ln_marg_v_likelihood(rs, data, Cov, prior=v_prior,
                                v_scatter=v_scatter)
    ln_p = [ln_dist_prob(rs[i], data[i, 2], np.sqrt(Cov[i, 2, 2]), rlim)
            for i in range(len(rs))]
    return ln_f + sum(ln_p)

class ComovingFML:
    """
    Helper class for computing the two fully marginalized likelihoods used for
    selecting candidate comoving stars.

    Parameters
    ----------
    gaia_data : ``GaiaData`` instance
        TODO
    v_prior : dict
        TODO
    v_scatter : `~astropy.units.Quantity`, optional [km/s]
        TODO
    rlim : `~astropy.units.Quantity`, optional [kpc]
        Maximum true distance to use for the assumed uniform space density prior
        over distance.
    n_dist_grid : int, iterable, optional
        The number of grid points to used when computing the marginalization
        over true distance. This can either be an integer (same number of grid
        points for both stars), or an iterable (number of grid points for each
        star). If not specified, this is estimated automatically for each star.
    """
    @u.quantity_input(v_scatter=u.km/u.s, rlim=u.kpc)
    def __init__(self, gaia_data, v_prior,
                 v_scatter=0*u.km/u.s, rlim=2*u.kpc, n_dist_grid=None):
        # Retrieve and convert data for each star
        data, Cov = data_to_vec_cov(gaia_data)
        if data.shape[0] != 2:
            raise ValueError("This class currently only supports pairs of "
                             "stars. You passed in data for {0} stars."
                             .format(data.shape[0]))

        if data.shape != Cov.shape[:2]:
            raise ValueError("Covariance matrix has invalid shape.")

        self.data = data
        self.Cov = Cov

        self.rlim = rlim
        self._rlim = self.rlim.to(u.pc).value

        self.v_scatter = v_scatter
        self._v_scatter = self.v_scatter.to(u.km/u.s).value

        self.v_prior = v_prior

        if n_dist_grid is None:
            # HACK: TODO: estimate this automatically using some heuristic
            n_dist_grid = (25, 25)

        elif isiterable(n_dist_grid):
            n_dist_grid = tuple(map(int, n_dist_grid))

        else:
            try:
                n_dist_grid = (int(n_dist_grid), int(n_dist_grid))
            except Exception as e:
                raise TypeError("Invalid specification of number of distance "
                                "grid points: {0}".format(n_dist_grid))

        self.n_dist_grid = n_dist_grid

        # parallax at index 2, parallax error: note that we drop the cross
        # terms!!!
        rgrid1 = get_rgrid(self.data[0, 2], np.sqrt(self.Cov[0, 2, 2]),
                           self._rlim, self.n_dist_grid[0])
        rgrid2 = get_rgrid(self.data[1, 2], np.sqrt(self.Cov[1, 2, 2]),
                           self._rlim, self.n_dist_grid[1])
        self.rgrids = [rgrid1, rgrid2]

    ###########################################################################
    # Likelihood 1:
    #
    def ln_H1_FML(self):
        ln_H = np.array([ln_dist_integrand([r1, r2], self.data, self.Cov,
                                           self.v_prior, self._rlim,
                                           self._v_scatter)
                         for r1 in self.rgrids[0] for r2 in self.rgrids[1]])
        ln_H = ln_H.reshape(self.n_dist_grid)

        # Simpson's rule, twice to do 2D integral:
        return np.log(simps(simps(np.exp(ln_H), x=self.rgrids[0], axis=0),
                            x=self.rgrids[1]))

    ###########################################################################
    # Likelihood 2:
    #
    def ln_Q(self, index):
        """
        Compute the log(Q) integral (Eq. 10) for the star with the specified
        index in the internal data array.
        """
        data = self.data[index]
        Cov = self.Cov[index]
        rgrid = self.rgrids[index]
        func_vals = [ln_dist_integrand([r], data[None], Cov[None],
                                       self.v_prior, self._rlim,
                                       self._v_scatter)
                     for r in rgrid]
        return np.log(simps(np.exp(func_vals), x=rgrid))

    def ln_H2_FML(self):
        ln_Q1 = self.ln_Q(index=0)
        ln_Q2 = self.ln_Q(index=1)
        return ln_Q1 + ln_Q2
