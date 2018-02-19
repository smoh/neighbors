from sqlalchemy import case, cast, Column, types
from sqlalchemy.types import Integer, REAL, SmallInteger, String
from sqlalchemy.schema import ForeignKey
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.hybrid import hybrid_property

import astropy.units as u
import astropy.coordinates as coords

from .config import Base
from . import np_adapt

class TGASSource(Base):
    __tablename__ = 'tgas_source'

    id = Column(Integer, primary_key=True)
    # row_index = Column('row_index', Integer, nullable=False)

    # TGAS columns
    solution_id = Column('solution_id', Integer, nullable=False)
    source_id = Column('source_id', Integer, nullable=False, unique=True)
    tycho2_id = Column('tycho2_id', String)
    hip = Column('hip', Integer)
    random_index = Column('random_index', Integer, nullable=False)
    ref_epoch = Column('ref_epoch', REAL, nullable=False)
    ra = Column('ra', REAL, nullable=False)
    ra_error = Column('ra_error', REAL, nullable=False)
    dec = Column('dec', REAL, nullable=False)
    dec_error = Column('dec_error', REAL, nullable=False)
    parallax = Column('parallax', REAL, nullable=False)
    parallax_error = Column('parallax_error', REAL, nullable=False)
    pmra = Column('pmra', REAL, nullable=False)
    pmra_error = Column('pmra_error', REAL, nullable=False)
    pmdec = Column('pmdec', REAL, nullable=False)
    pmdec_error = Column('pmdec_error', REAL, nullable=False)
    ra_dec_corr = Column('ra_dec_corr', REAL, nullable=False)
    ra_parallax_corr = Column('ra_parallax_corr', REAL, nullable=False)
    ra_pmra_corr = Column('ra_pmra_corr', REAL, nullable=False)
    ra_pmdec_corr = Column('ra_pmdec_corr', REAL, nullable=False)
    dec_parallax_corr = Column('dec_parallax_corr', REAL, nullable=False)
    dec_pmra_corr = Column('dec_pmra_corr', REAL, nullable=False)
    dec_pmdec_corr = Column('dec_pmdec_corr', REAL, nullable=False)
    parallax_pmra_corr = Column('parallax_pmra_corr', REAL, nullable=False)
    parallax_pmdec_corr = Column('parallax_pmdec_corr', REAL, nullable=False)
    pmra_pmdec_corr = Column('pmra_pmdec_corr', REAL, nullable=False)
    phot_g_n_obs = Column('phot_g_n_obs', SmallInteger, nullable=False)
    phot_g_mean_flux = Column('phot_g_mean_flux', REAL, nullable=False)
    phot_g_mean_flux_error = Column('phot_g_mean_flux_error', REAL, nullable=False)
    phot_g_mean_mag = Column('phot_g_mean_mag', REAL, nullable=False)

    l = Column('l', REAL, nullable=False)
    b = Column('b', REAL, nullable=False)
    # ecl_lon = Column('ecl_lon', REAL, nullable=False)
    # ecl_lat = Column('ecl_lat', REAL, nullable=False)

    # photometric data
    j = Column('j', REAL)
    j_err = Column('j_err', REAL)
    h = Column('h', REAL)
    h_err = Column('h_err', REAL)
    ks = Column('ks', REAL)
    ks_err = Column('ks_err', REAL)
    #: SDSS u mag
    u = Column('u', REAL)
    #: SDSS u mag error
    u_err = Column('u_err', REAL)
    #: SDSS g mag
    g = Column('g', REAL)
    #: SDSS g mag error
    g_err = Column('g_err', REAL)
    #: SDSS r mag
    r = Column('r', REAL)
    #: SDSS r mag error
    r_err = Column('r_err', REAL)
    #: SDSS i mag
    i = Column('i', REAL)
    #: SDSS i mag error
    i_err = Column('i_err', REAL)
    #: SDSS z mag
    z = Column('z', REAL)
    #: SDSS z mag error
    z_err = Column('z_err', REAL)
    #: wise w1 mag (from profile-fitting)
    w1 = Column('w1', REAL)
    w1_err = Column('w1_err', REAL)
    w2 = Column('w2', REAL)
    w2_err = Column('w2_err', REAL)
    w3 = Column('w3', REAL)
    w3_err = Column('w3_err', REAL)
    w4 = Column('w4', REAL)
    w4_err = Column('w4_err', REAL)

    @hybrid_property
    def skycoord(self):
        return coords.SkyCoord(ra=self.ra*u.deg, dec=self.dec*u.deg,
                               distance=1000./self.parallax*u.pc)
    @hybrid_property
    def row_dict(self):
        row = dict()
        for k in self.__dict__:
            if not k.startswith('_'):
                v = getattr(self, k)

                if hasattr(v, 'value'):
                    row[k] = v.value

                elif not isinstance(v, float):
                    continue

                else:
                    row[k] = v

        return row

    @hybrid_property
    def any_id(self):
        return 'HIP {:d}'.format(self.hip) if self.tycho2_id == '' else 'TYC '+self.tycho2_id

    @any_id.expression
    def any_id(TGASSource):
        return case({False: 'TYC '+TGASSource.tycho2_id,
                     True: 'HIP '+cast(TGASSource.hip, String)},
                    TGASSource.tycho2_id=='')


class RAVE(Base):
    __tablename__ = 'rave'

    id = Column(types.Integer, primary_key=True)
    RAVE_OBS_ID = Column(types.String, nullable=False)
    RAVEID = Column(types.String, nullable=False)
    RAdeg = Column(types.Float)
    DEdeg = Column(types.Float)
    Glon = Column(types.Float)
    Glat = Column(types.Float)
    HRV = Column(types.Float)
    eHRV = Column(types.Float)
    teff = Column(types.Float)
    logg = Column(types.Float)

    # relationship
    tgas_source_id = Column(types.Integer)


class Star(Base):
    __tablename__ = 'star'

    id = Column(types.Integer, primary_key=True)
    source_id = Column(
        Integer, ForeignKey('tgas_source.source_id'),
        nullable=False)
    tgas = relationship("TGASSource", backref=backref('oh17star', uselist=False))

    group_id = Column(Integer, ForeignKey('comovinggroup.id'), nullable=False)
    group = relationship("Group", backref='members')

    # def add_neighbor(self, *nodes):
    #     for node in nodes:
    #         Edge(self, node)
    #     return self
    #
    # def add_neighbors(self, node):
    #     if node not in self.neighbors():
    #         self.neighbors(node)

    def neighbors(self):
        all_nodes = [x.lower_node for x in self.higher_edges]
        all_nodes.extend([x.higher_node for x in self.lower_edges])
        return all_nodes

    def __str__(self):
        return "Star id={:d} group_id={:d}\n".format(
            self.id, self.group_id)


class Pair(Base):
    __tablename__ = 'comovingpair'

    lower_id = Column(types.Integer,
                      ForeignKey('star.id'),
                      primary_key=True)

    higher_id = Column(types.Integer,
                       ForeignKey('star.id'),
                       primary_key=True)

    lower_node = relationship(Star,
                              primaryjoin=lower_id==Star.id,
                              backref='lower_edges')

    higher_node = relationship(Star,
                               primaryjoin=higher_id==Star.id,
                               backref='higher_edges')

    #: Marginalized likelihood ratio in natural log
    lnL1L2 = Column('lnL1L2', types.REAL, nullable=False)

    def __init__(self, n1, n2, lnL1L2):
        # make sure lower_node.id < higher_node.id
        if not(isinstance(n1, Star) & isinstance(n2, Star)):
            raise ValueError("n1 and n2 should be Star instances")
        if n1.id < n2.id:
            self.lower_node = n1
            self.higher_node = n2
        else:
            self.lower_node = n2
            self.higher_node = n1
        self.lnL1L2 = lnL1L2


class Group(Base):
    """Model for Oh17 groups"""
    __tablename__ = 'comovinggroup'
    #: Group index
    id = Column(types.Integer, primary_key=True)
    #: Group size
    size = Column(Integer, nullable=False)
    #: Mean R.A. in degrees
    mean_ra = Column(REAL, nullable=False)
    #: Mean Declination in degrees
    mean_dec = Column(REAL, nullable=False)
    #: Mean distance in pc
    mean_distance = Column(REAL, nullable=False)


# class MWSC(Base):
#     __table__ = 'mwsc'
#     id = Column(types.Integer, primary_key=True)
#     name = Column('name', String)
#     type = Column('type', String)
#     glon = Column('glon', REAL)
#     glat = Column('glat', REAL)
#     d = Column('distance', REAL)
#     EBV = Column('EBV', REAL)
#     logage = Column('logage', REAL)

class Gaia17Cl(Base):
    __tablename__ = 'gaia17cl'
    id = Column(Integer, primary_key=True)
    name = Column('name', String, nullable=False, unique=True)
    melotte = Column('melotte', Integer)
    nmembers = Column('nmembers', Integer)

    members = relationship(
        "Gaia17ClMember", backref='cluster')

    def __repr__(self):
        # return "name={:s} N={:d}".format(self.name, self.nmembers)
        return "name={:s}".format(self.name)


class Gaia17ClMember(Base):
    __tablename__ = 'gaia17cl_members'
    id = Column(Integer, primary_key=True)
    source_id = Column('source_id', Integer, nullable=False)
    ra = Column(REAL, nullable=False)
    dec = Column(REAL, nullable=False)
    Gmag = Column(REAL, nullable=False)
    hd = Column('hd', Integer)

    cluster_name = Column(String, ForeignKey('gaia17cl.name'))
    # source_id = Column(Integer, ForeignKey('tgas_source.source_id'), nullable=False)
    # tgas_source = relationship('TGASSource')

    def __repr__(self):
        return "{:s} id={:d} ra={:.3f} dec={:.3f}\n".format(
            self.cluster.name, self.id, self.ra, self.dec)
