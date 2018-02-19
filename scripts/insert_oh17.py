"""
Insert Oh+17 catalog into database with associated TGAS data incl. photometry
from 2MASS, SDSS, and WISE.
"""
import warnings
from sqlalchemy import create_engine
from astropy.table import Table

from neighbors import data
from neighbors.db.config import Session, Base
from neighbors.db import model
from neighbors.db.model import \
    Star, Pair, Group, TGASSource

engine = create_engine("sqlite:///data/neighbors.db")

# Turn on ForeignKey checking for sqlite
# See: https://stackoverflow.com/questions/2614984/sqlite-sqlalchemy-how-to-enforce-foreign-keys
def _fk_pragma_on_connect(dbapi_con, con_record):
    dbapi_con.execute('pragma foreign_keys=ON')

from sqlalchemy import event
event.listen(engine, 'connect', _fk_pragma_on_connect)

Session.configure(bind=engine)
Base.metadata.create_all(engine)
session = Session()


star, pair, group = data.load_oh17()
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    star_tgas = Table.read("data/oh17star-gaia.vot", format='votable')
    # convert bytes columns by hand
    for col in star_tgas.columns:
        if star_tgas[col].dtype == object:
            print("Converting bytes column {:s} to string".format(col))
            star_tgas[col] = [s.decode('utf-8') if type(s) == bytes else '' for s in star_tgas[col]]

stars, pairs, groups = [], [], []
tgas_sources = []
for row in group.itertuples():
    groups.append(Group(
        id=row.id,
        size=row.size,
        mean_ra=row.mean_ra,
        mean_dec=row.mean_dec,
        mean_distance=row.mean_distance
    ))
for i, row in enumerate(star.itertuples()):
    tgas_source = TGASSource(**dict(
        solution_id=star_tgas[i]['solution_id'],
        source_id=star_tgas[i]['source_id'],
        hip=star_tgas[i]['hip'],
        tycho2_id=star_tgas[i]['tycho2_id'],
        random_index=star_tgas[i]['random_index'],
        ref_epoch=star_tgas[i]['ref_epoch'],
        ra=star_tgas[i]['ra'],
        ra_error=star_tgas[i]['ra_error'],
        dec=star_tgas[i]['dec'],
        dec_error=star_tgas[i]['dec_error'],
        parallax=star_tgas[i]['parallax'],
        parallax_error=star_tgas[i]['parallax_error'],
        pmra=star_tgas[i]['pmra'],
        pmra_error=star_tgas[i]['pmra_error'],
        pmdec=star_tgas[i]['pmdec'],
        pmdec_error=star_tgas[i]['pmdec_error'],
        ra_dec_corr=star_tgas[i]['ra_dec_corr'],
        ra_parallax_corr=star_tgas[i]['ra_parallax_corr'],
        ra_pmra_corr=star_tgas[i]['ra_pmra_corr'],
        ra_pmdec_corr=star_tgas[i]['ra_pmdec_corr'],
        dec_parallax_corr=star_tgas[i]['dec_parallax_corr'],
        dec_pmra_corr=star_tgas[i]['dec_pmra_corr'],
        dec_pmdec_corr=star_tgas[i]['dec_pmdec_corr'],
        parallax_pmra_corr=star_tgas[i]['parallax_pmra_corr'],
        parallax_pmdec_corr=star_tgas[i]['parallax_pmdec_corr'],
        pmra_pmdec_corr=star_tgas[i]['pmra_pmdec_corr'],
        phot_g_n_obs=star_tgas[i]['phot_g_n_obs'],
        phot_g_mean_flux=star_tgas[i]['phot_g_mean_flux'],
        phot_g_mean_flux_error=star_tgas[i]['phot_g_mean_flux_error'],
        phot_g_mean_mag=star_tgas[i]['phot_g_mean_mag'],
        l=star_tgas[i]['l'],
        b=star_tgas[i]['b'],
        # 2mass
        j=star_tgas[i]['j_m'],
        j_err=star_tgas[i]['j_msigcom'],
        h=star_tgas[i]['h_m'],
        h_err=star_tgas[i]['h_msigcom'],
        ks=star_tgas[i]['ks_m'],
        ks_err=star_tgas[i]['ks_msigcom'],
        # sdss
        u=star_tgas[i]['u_mag'],
        u_err=star_tgas[i]['u_mag_error'],
        g=star_tgas[i]['g_mag'],
        g_err=star_tgas[i]['g_mag_error'],
        r=star_tgas[i]['r_mag'],
        r_err=star_tgas[i]['r_mag_error'],
        i=star_tgas[i]['i_mag'],
        i_err=star_tgas[i]['i_mag_error'],
        z=star_tgas[i]['z_mag'],
        z_err=star_tgas[i]['z_mag_error'],
        # wise
        w1=star_tgas[i]['w1mpro'],
        w1_err=star_tgas[i]['w1mpro_error'],
        w2=star_tgas[i]['w2mpro'],
        w2_err=star_tgas[i]['w2mpro_error'],
        w3=star_tgas[i]['w3mpro'],
        w3_err=star_tgas[i]['w3mpro_error'],
        w4=star_tgas[i]['w4mpro'],
        w4_err=star_tgas[i]['w4mpro_error']
    ))
    tgas_sources.append(tgas_source)
    s = Star(
        id=row.row_id,
        source_id=row.tgas_source_id,
        group=groups[row.group_id])
    stars.append(s)

session.add_all(stars)
session.add_all(tgas_sources)
session.commit()
for row in pair.itertuples():
    p = Pair(stars[row.star1], stars[row.star2], row._5)
    pairs.append(p)
session.add_all(pairs)
session.commit()
