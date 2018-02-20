# -*- coding: utf-8 -*-
"""
module for generating fake data
"""
import numpy as np

from .. import tools


class FakeStar(object):
    def __init__(self, distance, ra, dec, vx, vy, vz):
        self.distance = distance
        self.ra = ra
        self.dec = dec
        self.vx = vx
        self.vy = vy
        self.vz = vz

        self.true_parallax = 1./distance
        self.true_vra, self.true_vdec, self.true_vr = \
            tools.project_at_position((vx,vy,vz), np.deg2rad(ra), np.deg2rad(dec))

        self.true_pmra = self.true_vra * self.true_parallax / 4.74
        self.true_pmdec = self.true_vdec * self.true_parallax / 4.74

    def observe(self, cov=None, snr=None):
        """
        cov (parallax, pmra, pmdec)
        """
        if snr:
            if cov:
                raise ValueError("cov and snr cannot be both used")
            cov = np.diag([(self.true_parallax/snr)**2,
                           (self.true_pmra/snr)**2,
                           (self.true_pmdec/snr)**2])
        parallax, pmra, pmdec = np.random.multivariate_normal(
            [self.true_parallax, self.true_pmra, self.true_pmdec], cov)
        return ObservedStar(self.ra, self.dec, parallax, pmra, pmdec, cov=cov)

    def angdist_to(self, star):
        """ Returns angular distance to star in degrees """
        return angdist_to(self.ra, self.dec, star.ra, star.dec)

    def __str__(self):
        return "\n".join([
            "FakeStar ra={0.ra:.3f}, dec={0.dec:.3f}, d={0.distance:.2f}, "
            "(vx, vy, vz) = ({0.vx:.2f},{0.vy:.2f},{0.vz:.2f})",
            "truth: parallax={0.true_parallax:.3f}, vra={0.true_vra:.2f} "
            "vdec={0.true_vdec:.2f}, pmra={0.true_pmra:.3f}, pmdec={0.true_pmdec:.3f}"
            ]).format(self)



class FakeEnsemble(object):
    """ class representing a fake ensemble of comoving stars """
    def __init__(self, d=0.1, ux=10., uy=10., uz=10., sigmau=0.1,
                 ra=45., dec=30., dangle=1., N=3):
        """
        Initialize a fake ensemble of comoving stars

        d : float
            distance
        ux, uy, uz : float
            comoving velocity
        sigmau : float
            internal velocity dispersion
        ra, dec : float
            sky coordinates in degrees
        dangle : float
            angular size in degrees
            star coordinates are sampled uniform random from (-dangle, dangle)
        N : int
            number of stars in ensemble
        """
        self.d = d
        self.ux = ux
        self.uy = uy
        self.uz = uz
        self.sigmau = sigmau
        self.ra = ra
        self.dec = dec
        self.dangle = dangle
        self.N = N
        self.stars = []

        for i in range(N):
            self.stars.append(self.sample())

    def __str__(self):
        return " ".join([
            "Ensemble N = {0.N:2d}",
            "ra, dec = {0.ra:.3f}, {0.dec:.3f}",
            "d = {0.d:.2f}",
            "(ux, uy, uz) = ({0.ux:.2f}, {0.uy:.2f}, {0.uz:.2f})",
            ]).format(self)

    def sample(self):
        """Sample a star from this fake Ensemble"""
        ra, dec = np.random.uniform(
            (self.ra-self.dangle,self.dec-self.dangle),
            (self.ra+self.dangle,self.dec+self.dangle))
        if self.sigmau >0:
            vx, vy, vz = np.random.normal((self.ux, self.uy, self.uz), self.sigmau)
        else:
            vx, vy, vz = self.ux, self.uy, self.uz
        d = np.random.uniform(0.95*self.d, 1.05*self.d)
        return FakeStar(d, ra, dec, vx, vy, vz)

    def angdist_to(self, star):
        return angdist_to(self.ra, self.dec, star.ra, star.dec)


class ObservedStar(object):
    """class representing astrometric source
    assumed units:
        ra, dec : deg
        parallax : mas
        pmra, pmdec : mas/yr
        rv : km/s
    """
    def __init__(self, ra, dec, parallax, pmra, pmdec, cov, rv=None, rv_error=None):
        self.ra = ra
        self.dec = dec
        self.parallax = parallax
        self.pmra = pmra
        self.pmdec = pmdec

        self.ra_rad, self.dec_rad = np.deg2rad(self.ra), np.deg2rad(self.dec)

        if cov.shape != (3,3):
            raise ValueError("cov must be 3x3 matrix")
        if (cov.T != cov).all():
            raise ValueError("cov must be symmetric")
        if (cov<0).any():
            raise ValueError("cov contains one or more negative values")
        self.cov = cov

        self.rv = rv
        self.rv_error = rv_error

        self.parallax_error = np.sqrt(cov[0,0])
        self.pmra_error = np.sqrt(cov[1,1])
        self.pmdec_error = np.sqrt(cov[2,2])

        self.parallax_snr = self.parallax / self.parallax_error
        self.pmra_snr = self.pmra / self.pmra_error
        self.pmdec_snr = self.pmdec / self.pmdec_error

    def project_at_position(self, vrec):
        """project velocity vector at star's position

        vvec : [vx, vy, vz]
        """
        return tools.project_at_position(vrec, self.ra_rad, self.dec_rad)

    def __str__(self):
        return " ".join([
            "ObservedStar ra={0.ra:.3f}, dec={0.dec:.3f}, parallax={0.parallax:.5f}, "
            "pmra={0.pmra:.3f} pmdec={0.pmdec:.3f}"
            ]).format(self)

    @property
    def distance(self):
        return 1./self.parallax #* 1e3

    @property
    def vra(self):
        return self.distance * self.pmra# / _KMSPC_TO_MASYR

    @property
    def vdec(self):
        return self.distance * self.pmdec# / _KMSPC_TO_MASYR

    @property
    def vr(self):
        return self.rv
