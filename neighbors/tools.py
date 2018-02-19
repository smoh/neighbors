import numpy as np
from numpy import pi, log, deg2rad, rad2deg

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

    Returns [v_ra, v_dec, v_r]
    """
    return get_projection_matrix(ra, dec).dot(vrec)


def deproject_at_position(vrec, ra, dec):
    """
    de-project [v_ra, v_dec, v_r] at (ra, dec)

    Returns [vx, vy, vz]
    """
    return get_projection_matrix(ra, dec).T.dot(vrec)
