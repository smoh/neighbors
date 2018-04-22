# Third-party
import astropy.units as u
import numpy as np

__all__ = ['data_to_vec_cov']


def data_to_vec_cov(gaia_data):
    """
    Convert the input ``GaiaData`` object into raw data and covariance
    arrays to be used in the comoving star analysis.
    """
    if not gaia_data._has_rv:
        rv = np.zeros(len(gaia_data))
    else:
        rv = gaia_data.radial_velocity.to(u.km/u.s).value

    vec = np.array([gaia_data.ra.to(u.radian).value,
                    gaia_data.dec.to(u.radian).value,
                    gaia_data.parallax.to(u.mas).value,
                    gaia_data.pmra.to(u.mas/u.yr).value,
                    gaia_data.pmdec.to(u.mas/u.yr).value,
                    rv]).T

    return vec, gaia_data.get_cov()
