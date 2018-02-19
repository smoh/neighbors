from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session
from sqlalchemy.ext.declarative import declarative_base

Session = scoped_session(sessionmaker(autoflush=True, autocommit=False))
Base = declarative_base()


def make_session():
    engine = create_engine("sqlite:////Users/semyeong/projects/neighbors/data/neighbors.db")
    Session.configure(bind=engine)
    session = Session()
    return session
