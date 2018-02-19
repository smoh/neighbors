# -*- coding: utf-8 -*-
"""
module for loading real and fake data
"""
import os
import numpy as np
import pandas as pd
import six

import fitsio
from astropy.io import fits
from astropy.table import Table
import astropy.coordinates as coord
import astropy.units as u

from . import tools
from .utils import cached_property

__all__ = [
    'load_mwsc', 'load_oh17', 'load_malo13',
    'TGASData', 'TGASStar'
    'AstrometricSource',
    'MWSC'
]

# TODO: auto-download if data does not exist
_DATADIR = os.path.dirname(os.path.dirname(__file__)[:-1])+'/data'

datasets = {
    'mwsc': 'J/A+A/585/A101',
    'malo13': 'J/ApJ/762/88',
    'gaia17cl': 'J/A+A/601/A19'
}

def load_mwsc():
    """Load Milky Way Star Cluster catalog in astropy Table"""
    dirname = _DATADIR + '/mwsc'
    table = Table.read(
        dirname+'/catalog.dat',
        readme=dirname+'/ReadMe',
        format='ascii.cds')
    return table


class MWSC(object):
    """class representing Milky Way Star Cluster catalog"""
    def __init__(self):
        self.table = load_mwsc()

    @cached_property
    def coords(self):
        return coord.Galactic(np.array(self.table['GLON'])*u.deg,
                              np.array(self.table['GLAT'])*u.deg,
                              np.array(self.table['d'])*u.pc)


def load_oh17():
    """Returns (star, pair, group) tables of Oh+17 catalog in pandas DataFrame
    """
    dirname = _DATADIR + '/oh17'
    star = pd.read_csv(dirname+"/star.csv")
    pair = pd.read_csv(dirname+"/pair.csv")
    group = pd.read_csv(dirname+"/group.csv")
    return star, pair, group


def load_tgas(columns='default'):
    """Returns TGAS in an array

    columns : list
      list of columns to read (passed to `fitsio.read`)
      Set to `None` to get all columns
      Set to 'default' to get astrometry-related columns only
    """
    if columns == 'default':
        columns = ['hip', 'tycho2_id', 'source_id',
                   'ra', 'dec', 'parallax', 'pmra', 'pmdec',
                   'phot_g_mean_mag', 'l', 'b']
    return fitsio.read(_DATADIR+'/stacked_tgas.fits', columns=columns)


def load_galah():
    """Returns GALAH+TGAS catalog in astropy table"""
    return Table.read(_DATADIR+"/galah/catalog.dat",
                      format='ascii.cds',
                      readme=_DATADIR+"/galah/ReadMe")


def load_malo13():
    """
    Returns table of Malo+2013 of bona-fide nearby moving group members
    """
    dirname = _DATADIR + '/malo13'
    malo = Table.read(
        dirname+'/table3.dat',
        readme=dirname+'/ReadMe',
        format='ascii.cds')
    return malo

def load_gaia17cl():
    dirname = _DATADIR + '/gaia17cl'
    table = Table.read(
        dirname+'/tabled.dat',
        readme=dirname+'/ReadMe',
        format='ascii.cds')
    return table


class TGAS(object):
    """ class representing TGAS data"""
    def __init__(self):
        self.data = load_tgas()

    @cached_property
    def skycoords(self):
        return coord.ICRS(self.data['ra']*u.deg,
                          self.data['dec']*u.deg,
                          1000./self.data['parallax']*u.pc)
    @cached_property
    def galcoords(self):
        return self.skycoords.transform_to('galactic')

    @cached_property
    def vra(self):
        return self.data['pmra']/self.data['parallax']*4.74

    @cached_property
    def vdec(self):
        return self.data['pmdec']/self.data['parallax']*4.74




class AstrometricSource(object):
    """ class representing an astrometric source """

    _unit_map = {
        'ra': u.degree,
        'dec': u.degree,
        'parallax': u.milliarcsecond,
        'pmra': u.milliarcsecond/u.year,
        'pmdec': u.milliarcsecond/u.year,
        'ra_error': u.degree,
        'dec_error': u.degree,
        'parallax_error': u.milliarcsecond,
        'pmra_error': u.milliarcsecond/u.year,
        'pmdec_error': u.milliarcsecond/u.year,
    }

    def __init__(self, data):
        self._data = data

    def __getattr__(self, name):
        # to prevent recursion errors:
        #   http://nedbatchelder.com/blog/201010/surprising_getattr_recursion.html
        if name == '_data':
            raise AttributeError()
        if name in AstrometricSource._unit_map:
            return self._data[name] * AstrometricSource._unit_map[name]
        elif name in self._data.dtype.names:
            return self._data[name]
        else:
            raise AttributeError("Object {} has no attribute '{}' and source data "
                                 "table has no column with that name.".format(self, name))


class TGASData(object):
    _unit_map = {
        'ra': u.degree,
        'dec': u.degree,
        'parallax': u.milliarcsecond,
        'pmra': u.milliarcsecond/u.year,
        'pmdec': u.milliarcsecond/u.year,
        'ra_error': u.degree,
        'dec_error': u.degree,
        'parallax_error': u.milliarcsecond,
        'pmra_error': u.milliarcsecond/u.year,
        'pmdec_error': u.milliarcsecond/u.year,
    }

    def __init__(self, filename_or_data, rv=None, rv_err=None):

        # radial velocities
        if rv is not None:
            if rv_err is None:
                raise ValueError("If radial velocity is provided, you must also "
                                 "provide an error.")

            if not hasattr(rv, 'unit') or not hasattr(rv_err, 'unit'):
                raise TypeError("Input radial velocity and error must be an Astropy "
                                "Quantity object.")

            elif not rv.unit.is_equivalent(u.km/u.s) or not rv_err.unit.is_equivalent(u.km/u.s):
               raise u.UnitsError("Radial velocity unit is not convertible to km/s!")

            self._rv = rv.to(u.km/u.s).value
            self._rv_err = rv_err.to(u.km/u.s).value

        else:
            self._rv = 0. # need to do this so y can have float type
            self._rv_err = None

        # TODO: maybe support memory-mapping here?
        if isinstance(filename_or_data, six.string_types):
            self._data = np.array(fits.getdata(filename_or_data, 1))

        else:
            self._data = np.array(filename_or_data)

    def __getattr__(self, name):
        # to prevent recursion errors:
        #   http://nedbatchelder.com/blog/201010/surprising_getattr_recursion.html
        if name == '_data':
            raise AttributeError()

        if name in TGASData._unit_map:
            return self._data[name] * TGASData._unit_map[name]

        elif name in self._data.dtype.names:
            return self._data[name]

        else:
            raise AttributeError("Object {} has no attribute '{}' and source data "
                                 "table has no column with that name.".format(self, name))

    def __getitem__(self, slc):
        sliced = self._data[slc]

        if self._rv_err is not None:
            rv = self._rv[slc]
            rv_err = self._rv_err[slc]
        else:
            rv = None
            rv_err = None

        if sliced.ndim == 0: # this is only one row
            return TGASStar(row=sliced, rv=rv, rv_err=rv_err)

        else: # many rows
            return TGASData(sliced, rv=rv, rv_err=rv_err)

    def __len__(self):
        return len(self._data)

    @property
    def rv(self):
        return self._rv*u.km/u.s

    # Other convenience methods
    def get_distance(self, lutz_kelker=True):
        """
        Return the distance with or without the Lutz-Kelker correction.
        """

        if lutz_kelker:
            snr = self._data['parallax'] / self._data['parallax_error']
            tmp = self._data['parallax'] * (0.5 + 0.5*np.sqrt(1 - 16/snr**2))

        else:
            tmp = self._data['parallax']

        return 1000./tmp * u.pc

    def get_vtan(self, lutz_kelker=True):
        """
        Return the tangential velocity computed using the proper motion
        and distance.
        """
        d = self.get_distance(lutz_kelker=lutz_kelker)
        vra = (self.pmra * d).to(u.km/u.s, u.dimensionless_angles()).value
        vdec = (self.pmdec * d).to(u.km/u.s, u.dimensionless_angles()).value
        return np.vstack((vra, vdec)).T * u.km/u.s

    def get_coord(self, lutz_kelker=True):
        """
        Return an `~astropy.coordinates.SkyCoord` object to represent
        all coordinates.
        """
        return coord.SkyCoord(ra=self.ra, dec=self.dec,
                              distance=self.get_distance(lutz_kelker=lutz_kelker))

    @property
    def parallax_snr(self):
        return self.parallax / self.parallax_error


class TGASStar(TGASData):
    def __init__(self, row, rv=None, rv_err=None):
        self._data = row
        self._cov = None # for caching
        self._Cinv = None # for caching

        # radial velocities
        if rv is not None:
            if rv_err is None:
                raise ValueError("If radial velocity is provided, you must also "
                                 "provide an error.")
            self._rv = rv.to(u.km/u.s).value
            self._rv_err = rv_err.to(u.km/u.s).value
        else:
            self._rv = 0. # need to do this so y can have float type
            self._rv_err = None

    def __len__(self):
        return 1

    def __getitem__(self, slc):
        object.__getitem__(self, slc)

    def __str__(self):
        infostr = '\n'.join([
            # 'index    = %i' %(i),
            'ra       = %s' % (self.ra),
            'dec      = %s' % (self.dec),
            'parallax = %s (snr = %.1f)' % (self.parallax, self.parallax_snr),
            'pmra     = %s (snr = %.1f)' % (self.pmra, self.pmra/self.pmra_error),
            'pmdec    = %s (snr = %.1f)' % (self.pmdec, self.pmdec/self.pmdec_error),
            'dist vra vdec = %s %s' % (self.get_distance(), self.get_vtan()),
        ])
        return infostr

    def get_cov(self):
        """
        The Gaia TGAS data table contains correlation coefficients and standard
        deviations for (ra, dec, parallax, pm_ra, pm_dec), but for most analysis
        we need covariance matrices. This converts the Gaia table into covariance
        matrix. If a radial velocity was specified on creation, this also contains
        the radial velocity variance. The base units are:
        [deg, deg, mas, mas/yr, mas/yr, km/s]
        """

        if self._cov is not None:
            return self._cov

        names = ['ra', 'dec', 'parallax', 'pmra', 'pmdec']

        C = np.zeros((6,6))

        # pre-load the diagonal
        for i,name in enumerate(names):
            full_name = "{}_error".format(name)
            C[i,i] = self._data[full_name]**2

        for i,name1 in enumerate(names):
            for j,name2 in enumerate(names):
                if j <= i:
                    continue
                full_name = "{}_{}_corr".format(name1, name2)
                C[i,j] = self._data[full_name] * np.sqrt(C[i,i]*C[j,j])
                C[j,i] = self._data[full_name] * np.sqrt(C[i,i]*C[j,j])

        if self._rv_err is not None:
            C[5,5] = self._rv_err**2

        self._cov = C
        return self._cov



def angdist_to(ra1, dec1, ra2, dec2):
    """Return angular distance in degrees"""
    ra1, dec1 = np.deg2rad(ra1), np.deg2rad(dec1)
    ra2, dec2 = np.deg2rad(ra2), np.deg2rad(dec2)
    return np.rad2deg(np.arccos(np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2)))


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



