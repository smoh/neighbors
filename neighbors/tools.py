import numpy as np
from numpy import pi, log, deg2rad, rad2deg

__all__ = [
    "sample_uniform_sphere",
    "sample_on_sphere",
    "xyz_to_rradec",
    "get_projection_matrix",
    "project_at_position",
    "deproject_at_position"
]

def sample_uniform_sphere(rmax, size=1):
    """
    Randomly sample points inside a uniform density sphere

    rmax : float
        size of the sphere
    size : int, optional, default: 1
        number of samples

    Returns (size, 3) array of (X,Y,Z) of random points inside rmax
    """
    U = np.random.uniform(size=size)
    x = np.random.normal(size=(size,3))
    return rmax*U[:,None]**(1./3.)*x/np.linalg.norm(x,axis=1)[:,None]


def sample_on_sphere(size=1):
    """
    Randomly sample points on a unit sphere (i.e., random direction)

    size : int, optional, default: 1
        number of samples

    Returns (X,Y,Z) of vectors. The size of all vectors is unity.

    See also:
    http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    """
    phi = np.random.uniform(0,np.pi*2, size=size)
    costheta = np.random.uniform(-1,1, size=size)
    theta = np.arccos(costheta)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return x, y, z


def xyz_to_rradec(x, y, z):
    """
    Transform rectangular coordinates to spherical coordinates (r, ra, dec).
    The angles (ra, dec) are in degrees.

    Note that the only difference with usual spherical coordinates is that
    the angle *dec* is 0 at z=0, and lies between -pi/2 to pi/2.

    Returns (r, ra, dec); ra, dec in radians.
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    dec = np.rad2deg(np.arcsin(z/r))
    ra = np.rad2deg((np.arctan2(y, x) + 2*np.pi) % (2*np.pi))
    return r, ra, dec


def get_projection_matrix(ra, dec):
    """
    Returns the projection matrix x, y, z -> ra, dec, r.
    If ra, dec is scalar, the returned matrix has shape (3,3).
    IF ra, dec is n-size array, the returned matrix has shape (n, 3, 3).

    ra, dec : scalar or array-like
        in radians
    """
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)
    if ra.size != dec.size:
        raise ValueError("ra and dec must have same size")
    s = []
    for cra, cdec in zip(ra, dec):
        s.append(
            np.array([
                [-np.sin(cra), np.cos(cra), 0.],
                [-np.sin(cdec)*np.cos(cra), -np.sin(cdec)*np.sin(cra), np.cos(cdec)],
                [np.cos(cdec)*np.cos(cra), np.cos(cdec)*np.sin(cra), np.sin(cdec)]
            ]))
    return np.array(s).squeeze()


def project_at_position(vrec, ra, dec):
    """
    project [vx, vy, vz] at (ra, dec)

    ra, dec : scalar or array-like
        in radians

    Returns [v_ra, v_dec, v_r]
    """
    return get_projection_matrix(ra, dec).dot(vrec)


def deproject_at_position(vrec, ra, dec):
    """
    de-project [v_ra, v_dec, v_r] at (ra, dec)

    Returns [vx, vy, vz]
    """
    return get_projection_matrix(ra, dec).T.dot(vrec)


def make_cov(d):
    """
    Generate (parallax, pmra, pmdec) covariance matrix
    from Gaia data for a single star
    """
    cov = np.zeros([3,3])
    cov[0,0] = d['parallax_error']**2
    cov[1,1] = d['pmra_error']**2
    cov[2,2] = d['pmdec_error']**2
    cov[0,1] = d['parallax_pmra_corr'] * d['parallax_error'] * d['pmra_error']
    cov[0,2] = d['parallax_pmdec_corr'] * d['parallax_error'] * d['pmdec_error']
    cov[1,2] = d['pmra_pmdec_corr'] * d['pmra_error'] * d['pmdec_error']
    cov[1,0] = d['parallax_pmra_corr'] * d['parallax_error'] * d['pmra_error']
    cov[2,0] = d['parallax_pmdec_corr'] * d['parallax_error'] * d['pmdec_error']
    cov[2,1] = d['pmra_pmdec_corr'] * d['pmra_error'] * d['pmdec_error']
    return cov
