# -*- coding: utf-8 -*-
"""
module for generating fake data
"""
import numpy as np
import pandas as pd

from ..tools import sample_uniform_sphere, xyz_to_rradec, make_cov

__all__ = ['MockGroup']

class MockGroup(object):
    """class representing mock comovoing group"""
    def __init__(self, N, x_ref, v_ref, sigv, rmax):
        """

        N : int
            number of stars
        x_ref : tuple
            (x, y, z) of reference position
        v_ref : tuple
            mean velocity (vx, vy, vz) of the group
        sigv : float
            isotropic velocity dispersion
        rmax : float
            maximum physical size
        """

        self.N = N
        self.x_ref = x_ref
        self.v_ref = v_ref
        self.sigv = sigv
        self.rmax = rmax

        xr, yr, zr = x_ref
        vx, vy, vz = v_ref

        vi = np.random.multivariate_normal([vx, vy, vz], np.eye(3)*sigv**2, size=N)

        xi, yi, zi = sample_uniform_sphere(rmax, size=N).T
        xi, yi, zi = xi+xr, yi+yr, zi+zr

        ri, rai, deci = xyz_to_rradec(xi, yi, zi)

        vrai = -np.sin(np.deg2rad(rai)) * vi[:,0] + np.cos(np.deg2rad(rai)) * vi[:,1]
        vdeci = -np.sin(np.deg2rad(deci))*np.cos(np.deg2rad(rai)) * vi[:,0] \
                -np.sin(np.deg2rad(deci))*np.sin(np.deg2rad(rai)) * vi[:,1] \
                +np.cos(np.deg2rad(deci)) * vi[:,2]
        vri = np.cos(np.deg2rad(deci))*np.cos(np.deg2rad(rai)) * vi[:,0] \
                +np.cos(np.deg2rad(deci))*np.sin(np.deg2rad(rai)) * vi[:,1] \
                +np.sin(np.deg2rad(deci)) * vi[:,2]

        pmrai = vrai / (ri/1000.) / 4.74
        pmdeci = vdeci / (ri/1000.) / 4.74

        self.stars = pd.DataFrame.from_items((
            ('ra', rai),
            ('dec', deci),
            ('parallax', 1000./ri),
            ('pmra', pmrai),
            ('pmdec', pmdeci)
        ))

    def observe(self, error_table=None):
        """
        error_table : pd.DataFrame
            must have these columns
            ['parallax_error', 'pmra_error', 'pmdec_error',
            'parallax_pmra_corr', 'parallax_pmdec_corr', 'pmra_pmdec_corr']

        Returns pd.DataFrame of TGAS-like fake-data table
        with errors randomly attached from the error_table.
        """
        error_columns = [
            'parallax_error', 'pmra_error', 'pmdec_error',
            'parallax_pmra_corr', 'parallax_pmdec_corr', 'pmra_pmdec_corr']
        idx = np.random.randint(0, len(error_table), size=self.N)
        df = pd.merge(
            self.stars,
            error_table.loc[idx, error_columns].reset_index(drop=True),
            left_index=True, right_index=True)

        def noisify(row):
            cov = make_cov(row)
            return tuple(np.random.multivariate_normal(
                [row.parallax, row.pmra, row.pmdec],
                cov))
        df['parallax'], df['pmra'], df['pmdec'] = \
            zip(*df.apply(noisify, axis=1))

        return df
