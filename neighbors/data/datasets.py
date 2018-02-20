# -*- coding: utf-8 -*-
"""
module for loading local datasets quickly
"""
import os
import numpy as np
import pandas as pd

import fitsio
from astropy.table import Table
import astropy.coordinates as coord
import astropy.units as u

from ..utils import cached_property

__all__ = [
    'load_mwsc', 'load_oh17', 'load_malo13',
    'load_tgas', 'load_galah', 'load_gaia17cl',
    'MWSC'
]

# TODO: auto-download if data does not exist?
_DATADIR = os.path.dirname(os.path.dirname(__file__)[:-1])+'/data'


def load_mwsc():
    """Load Milky Way Star Cluster catalog in astropy Table"""
    dirname = _DATADIR + '/mwsc'
    table = Table.read(
        dirname+'/catalog.dat',
        readme=dirname+'/ReadMe',
        format='ascii.cds')
    return table


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


class MWSC(object):
    """class representing Milky Way Star Cluster catalog"""
    def __init__(self):
        self.table = load_mwsc()

    @cached_property
    def coords(self):
        return coord.Galactic(np.array(self.table['GLON'])*u.deg,
                              np.array(self.table['GLAT'])*u.deg,
                              np.array(self.table['d'])*u.pc)
