import os
import numpy as np
from sqlalchemy import create_engine
from astropy.table import Table
from astroquery.gaia import Gaia

from neighbors import data
from neighbors.db.config import Session, Base
from neighbors.db.model import TGASSource, Gaia17Cl, Gaia17ClMember

# Gaia+2017 cluster catalog
if 0:
    gaia17_cl = Table.read(
        data._DATADIR + '/gaia17cl/table1.dat',
        readme=data._DATADIR + '/gaia17cl/ReadMe',
        format='ascii.cds')
    gaia17_mem = Table.read(
        data._DATADIR + '/gaia17cl/tabled.dat',
        readme=data._DATADIR + '/gaia17cl/ReadMe',
        format='ascii.cds')

    # Query GaiaArchive using member source ids
    # import tempfile
    # f = tempfile.NamedTemporaryFile()
    # gaia17_mem[['Source']].write(f, format='votable')
    # query = "SELECT gaiadr1.tgas_source.* FROM tap_upload.tmp LEFT JOIN gaiadr1.tgas_source ON gaiadr1.tgas_source.source_id = tap_upload.tmp.source"
    # j = Gaia.launch_job(query=query, upload_resource=f.name,
    #                 upload_table_name='tmp', verbose=True)
    # result = j.get_results()
    # f.close()
    # print((result['source_id'] == gaia17_mem['Source']).all())
    # print(result.colnames)

    clusters, members, tgas_sources = [], [], []
    # tgassource_columns = [
    #     str(c).split('.')[1] for c in TGASSource.__table__.columns]
    # tgassource_columns.remove('id')
    # for row in result:
    #     tgas_sources.append(TGASSource(**{
    #             col : row[col] for col in tgassource_columns}))
    #     session.add(tgas_sources[-1])
    #     try:
    #         session.commit()
    #         session.flush()
    #     except:
    #         print(row['source_id'])
    for row in gaia17_cl:
        cl = Gaia17Cl(
                name=row['Cluster'],
                melotte=row['Mel'],
                nmembers=row['Memb'])
        clusters.append(cl)
        for mem in gaia17_mem[gaia17_mem['Cluster'] == row['Cluster']]:
            members.append(
                Gaia17ClMember(
                    cluster=cl,
                    source_id=mem['Source'],
                    ra=mem['RAdeg'],
                    dec=mem['DEdeg'],
                    Gmag=mem['Gmag'],
                    hd=mem['HD']
                ))
    session.add_all(clusters)
    session.add_all(members)
    session.commit()

    cl = session.query(Gaia17Cl).first()
    print(cl)
    print(cl.members[0], cl.members[1])
