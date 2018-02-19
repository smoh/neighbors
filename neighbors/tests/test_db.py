import pytest
import numpy as np
from sqlalchemy import create_engine

from neighbors import data
from neighbors.db.model import Star, Group, Pair, TGASSource
from neighbors.db.config import Session, Base

@pytest.fixture(scope='session')
def db():
    engine = create_engine("sqlite:///:memory:")
    Base.metadata.create_all(engine)
    return engine

@pytest.fixture(scope='function')
def session(db):
    Session.configure(bind=db)
    session = Session()
    return session

tgassource_columns = [
    str(c).split('.')[1] for c in TGASSource.__table__.columns]

def test_oh17(session):
    star, pair, group = data.load_oh17()
    stars, groups = [], []
    for row in group.itertuples():
        groups.append(Group(id=row.id, size=row.size))
    for row in star.itertuples():
        s = Star(row_index=row.row_id, source_id=row.tgas_source_id,
                 group=groups[row.group_id])
        stars.append(s)
    session.add_all(stars)
    session.flush()
    for row in pair.itertuples():
        p = Pair(stars[row.star1], stars[row.star2], row._5)

    row = star.iloc[3264]
    assert [row.group_id] == session.query(Star.group_id).filter(Star.source_id==row.tgas_source_id).all()
    s1 = session.query(Star).filter(Star.source_id==row.tgas_source_id).all()[0]
    assert s1.group.size == len(s1.group.members)


# def test_db(session):
#     star, pair, group = data.load_oh17()
#     tgas = data.TGASData("/Users/semyeong/projects/gaia-wide-binaries/data/stacked_tgas.fits")
#     stars, tgassources = [], []
#     for s in star.iloc[:5].itertuples():
#         stars.append(
#             Star(
#                 row_index=s.row_id,
#                 source_id = s.tgas_source_id
#             ))

#         #Find entry in TGAS by source_id
#         tgas_idx = np.where(tgas.source_id == s.tgas_source_id)[0]
#         assert len(tgas_idx) == 1, \
#             "Panic: {:d} TGAS source found for star {:d}".format(
#                 len(tgas_idx), s.row_id)
#         tgas_idx = tgas_idx[0]
#         tgas_row = tgas._data[tgas_idx]
#         d = {}
#         for name in tgassource_columns:
#             if name in tgas_row.dtype.names:
#                 d[name] = tgas_row[name]
#         tgassources.append(
#             TGASSource(**d))
#     session.add_all(stars)
#     session.commit()

