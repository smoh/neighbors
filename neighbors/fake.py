"""
Module for generating fake binary samples
"""
import numpy as np
import astropy.units as u
import astropy.constants as c

from .tools import get_projection_matrix

__all__ = [
    'sample_powerlaw',
    'sample_uniform_sphere',
    'sample_on_sphere',
    'xyz_to_rradec',
    'BinarySampler'
]

def lognuniform(low=0, high=1, size=None, base=np.exp(1)):
    return np.power(base, np.random.uniform(low, high, size))

def sample_powerlaw(n, a, b, size=1):
    """
    Randomly sample a variable x with power-law pdf
        pdf(x) \propto x^n; x in [a, b)

    n : float
        power of the pdf
    a : float
        lower limit of the random variable
    b : float
        upper limit of the random variable
    size : int, optional, default: 1
        size of the sample

    Returns 1d array of size *size*

    Note: This is basically inverse transform sampling; only for power-law
    cases the inverse transform is analytic so it's straightforward to
    implement without numerical integration or interpolation.
    """
    r = np.random.random(size=size)
    if n == -1:
        ag, bg = np.log(a), np.log(b)
        return np.exp((ag + (bg-ag)*r))
    else:
        ag, bg = a**(n+1), b**(n+1)
        return (ag + (bg-ag)*r)**(1./(n+1))


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
    Transform rectangular coordinates to spherical coordinates (r, ra, dec)

    Note that the only difference with usual spherical coordinates is that
    the angle *dec* is 0 at z=0, and lies between -pi/2 to pi/2.

    Returns (r, ra, dec); ra, dec in radians.
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    dec = np.arcsin(z/r)
    ra = (np.arctan2(y, x) + 2*np.pi) % (2*np.pi)
    return r, ra, dec


class BinarySampler(object):
    """class for sampling random samples of binaries"""

    @u.quantity_input(mina=u.pc, maxa=u.pc)
    def __init__(self, powera=-1., mina=1e-3*u.pc, maxa=1*u.pc, size=1):
        """
        Make a fake sample of binaries with the following assumptions:
        - assumes equal solar mass
        - power-law distribution of semi-major axis
        - thermal distribution of eccentriciy, pdf(e) \propto e^2
        - random in inclination and true anomaly
        - distributed uniformly within a sphere of *rmax*
        - Gaussian velocity distribution with dispersion *sigma*

        powera : float, optional, default: -1
            power of the semi-major axis distribution
        mina, maxa : astropy Quantity, optional
            minimum and maximum semi-major axis of binaries
            The default is from 0.001 pc to 1 pc.
        size : int, optional, default: 1
            size of the sample
        """
        # NOTE: all angles are in radians
        self.size = size
        # semimajor axis pdf(a) \propto a^*powera*
        self.a = sample_powerlaw(powera, mina.value, maxa.value, size=size)*mina.unit
        # self.a = lognuniform(np.log10(10*u.R_sun.to(u.pc)), 1, size=size, base=10)*u.pc
        # eccentricity pdf(e) \propto e^2
        self.e = np.sqrt(np.random.uniform(size=size))
        # random true anomaly
        self.f = np.random.uniform(size=size)*np.pi*2.
        # period
        self.T = (np.sqrt(4*np.pi**2/(c.G * 2*u.solMass) * self.a**3)).to(u.yr)

        # random inclination
        self.inc = np.arccos(1-2*np.random.uniform(size=size))

        # velocity of relative vector in orbital plane
        # x-axis directs from focus to pericenter
        self.vx = (- 2*np.pi/self.T * self.a / np.sqrt(1 - self.e**2) * np.sin(self.f)).to(u.km/u.s)
        self.vy = (- 2*np.pi/self.T * self.a / np.sqrt(1 - self.e**2) * (self.e + np.sin(self.f))).to(u.km/u.s)
        # size of relative vector
        self.r = self.a * (1-self.e**2) / (1 + self.e*np.cos(self.f))

    def project_relative_vector(self):
        """
        """
        projected_separation = np.abs(self.r*np.cos(self.inc))
        projected_orbital_velocity = np.abs(np.hypot(self.vx, self.vy)*np.cos(self.inc))
        return projected_separation, projected_orbital_velocity

    def project(self, dmax=200, sigv=20):
        """"""
        # velocity of star 1 and 2 in the center of mass frame
        vx1, vy1 = self.vx*0.5, self.vy*0.5
        vx2, vy2 = -self.vx*0.5, -self.vy*0.5

        # position of to star 1
        x1, y1, z1 = sample_uniform_sphere(dmax, size=self.size).T*u.pc
        r1, ra1, dec1 = xyz_to_rradec(x1.value, y1.value, z1.value)
        # vector from star 1 to star 2 (random orientation)
        rx, ry, rz = sample_on_sphere(size=self.size) * self.r
        x2, y2, z2 = x1+rx, y1+ry, z1+rz
        r2, ra2, dec2 = xyz_to_rradec(x2.value, y2.value, z2.value)

        vcm_x, vcm_y, vcm_z = np.random.normal(0, sigv, size=(3, self.size))*u.km/u.s
        vx1_tot, vy1_tot, vz1_tot = vx1*np.cos(self.inc)+vcm_x, vy1*np.cos(self.inc)+vcm_y, vcm_z
        vx2_tot, vy2_tot, vz2_tot = vx2*np.cos(self.inc)+vcm_x, vy2*np.cos(self.inc)+vcm_y, vcm_z

        Vx1, Vy1, Vz1 = vx1_tot, vy1_tot, vz1_tot
        Vx2, Vy2, Vz2 = vx2_tot, vy2_tot, vz2_tot

        A1 = get_projection_matrix(ra1, dec1)
        A2 = get_projection_matrix(ra2, dec2)

        vradecr1 = np.einsum('ijk,ki->ij', A1,
                             np.hstack([Vx1[:,None], Vy1[:,None], Vz1[:,None]]).value.T)
        vradecr2 = np.einsum('ijk,ki->ij', A2,
                             np.hstack([Vx2[:,None], Vy2[:,None], Vz2[:,None]]).value.T)
        dvradecr = vradecr1-vradecr2
        dvperp = np.hypot(dvradecr[:,0], dvradecr[:,1])

        vradecr1 = np.einsum(
            'ijk,ki->ij', A1,
            np.hstack([
                (vx1*np.cos(self.inc))[:,None],
                (vy1*np.cos(self.inc))[:,None],
                (vx1*0)[:,None]
            ]).value.T)
        vradecr2 = np.einsum(
            'ijk,ki->ij', A2,
            np.hstack([
                (vx2*np.cos(self.inc))[:,None],
                (vy2*np.cos(self.inc))[:,None],
                (vx1*0)[:,None]
            ]).value.T)
        dvradecr = vradecr1-vradecr2
        dvperp_wrong = np.hypot(dvradecr[:,0], dvradecr[:,1])

        projected_separation = np.abs(self.r*np.cos(self.inc))
        return projected_separation, dvperp, vars()
