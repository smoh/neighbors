# Third-party
import numpy as np

# Project
from .units import dist_pm_to_vel, vel_to_dist_pm

__all__ = ['get_u_vec', 'get_tangent_basis', 'get_y_Cinv', 'get_M',
           'get_Ainv_nu_Delta']


def get_u_vec(lon, lat):
    """
    Given two sky coordinates as a longitude and latitude (RA, Dec),
    return a unit vector that points in the direction of the sky position.

    Sky positions should be in radians!

    Parameters
    ----------
    lon : numeric [rad]
        Longitude in radians.
    lat : numeric [rad]
        Latitude in radians.

    Returns
    -------
    u_hat : `numpy.ndarray`
        Unit 3-vector.

    """
    u_hat = np.array([np.cos(lon) * np.cos(lat),
                      np.sin(lon) * np.cos(lat),
                      np.sin(lat)])
    return u_hat


def get_tangent_basis(ra, dec):
    """
    Row vectors are the tangent-space basis at (ra, dec) on the unit sphere.

    Parameters
    ----------
    lon : numeric, array_like [rad]
        Longitude in radians.
    lat : numeric, array_like [rad]
        Latitude in radians.
    """
    ra = np.array(ra)
    dec = np.array(dec)

    M = np.zeros(ra.shape + (3, 3))
    M[..., 0, 0] = -np.sin(ra)
    M[..., 0, 1] = np.cos(ra)

    M[..., 1, 0] = -np.sin(dec) * np.cos(ra)
    M[..., 1, 1] = -np.sin(dec) * np.sin(ra)
    M[..., 1, 2] = np.cos(dec)

    M[..., 2, 0] = np.cos(dec) * np.cos(ra)
    M[..., 2, 1] = np.cos(dec) * np.sin(ra)
    M[..., 2, 2] = np.sin(dec)

    return M.squeeze()


def get_y_Cinv(d, pm_rv, Cov, jitter=0.):
    """
    Construct the vector `y`, and the inverse covariance matrix for the same
    scaled parameters, given true distance d.

    Parameters
    ----------
    d : numeric [pc]
        Distance in parsecs.
    pm_rv : array_like
        Iterable with pmra, pmdec, rv in units [mas/yr, mas/yr, km/s].
    Cov : array_like
        3x3 covariance matrix for just proper motion and radial velocity.
        Here we break the covariance between parallax and proper motion.
        We also assume that the radial velocity is uncorrelated with the proper
        motions.
    jitter : numeric [km/s]
        Extra velocity dispersion to add to the covariance matrix. This acts as
        a "tolerance" for what we consider to be "comoving."

    """
    y = np.array([d * pm_rv[0] * dist_pm_to_vel,
                  d * pm_rv[1] * dist_pm_to_vel,
                  pm_rv[2]])

    Cov = np.array(Cov, copy=True)

    Cov[0, :] *= d * dist_pm_to_vel # pmra
    Cov[1, :] *= d * dist_pm_to_vel # pmdec
    Cov[:, 0] *= d * dist_pm_to_vel # pmra
    Cov[:, 1] *= d * dist_pm_to_vel # pmdec

    # Add extra variance
    Cov += np.eye(Cov.shape[0]) * jitter**2

    Cinv = np.zeros_like(Cov)
    Cinv[:2, :2] = np.linalg.inv(Cov[:2, :2])
    Cinv[2, 2] = 1 / Cov[2, 2]

    return y, Cinv


def get_M(ra, dec):
    """
    Construct the projection matrix M, which should have shape ``(3*nstars,
    3)``, where ``nstars = len(ra)``.

    Parameters
    ----------
    ra : array_like [rad]
    dec : array_like [rad]
    """
    return np.vstack(get_tangent_basis(ra, dec))


def get_Ainv_nu_Delta(d, M_dirty, Cinv_dirty, y_dirty, Vinv):
    """
    Dirty because they may contain missing data that we mask out.

    Parameters
    ----------
    d : numeric [pc]
        Distance.
    M_dirty : array_like
        Transformation matrix.
    Cinv_dirty : array_like [1/(km/s)^2]
    y_dirty : array_like [km/s]
    Vinv : array_like
        1/(km/s)^2

    """
    d = np.atleast_1d(d)

    # If inverse variances are zero, some data is missing
    mask = np.diag(Cinv_dirty) != 0
    Cinv = Cinv_dirty[mask]
    Cinv = Cinv[:, mask]
    _, log_detCinv = np.linalg.slogdet(Cinv / (2*np.pi))

    M = M_dirty[mask]
    y = y_dirty[mask]

    # using ji vs. ij does the transpose of M
    Ainv = np.einsum('ji,jk,ks->is', M, Cinv, M) + Vinv

    # using ji vs. ij does the transpose
    Bb = np.einsum('ji,jk,k->i', M, Cinv, y)
    nu = -np.linalg.solve(Ainv, Bb)

    sgn, log_detVinv = np.linalg.slogdet(Vinv / (2*np.pi))

    yT_Cinv_y = np.einsum('i,ji,j->', y, Cinv, y)
    nuT_Ainv_nu = np.einsum('i,ji,j->', nu, Ainv, nu)
    Delta = (-sum([2*np.log(dd) for dd in d]) - 0.5*log_detCinv -
             0.5*log_detVinv + 0.5*yT_Cinv_y - 0.5*nuT_Ainv_nu)

    return Ainv, nu, Delta
