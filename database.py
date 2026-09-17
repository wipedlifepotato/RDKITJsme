from sqlalchemy import Column, String, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

Base = declarative_base()

class CachedName(Base):
    __tablename__ = "name_cache"
    smiles = Column(String, primary_key=True, index=True)
    name = Column(String, nullable=False)

engine = None
SessionLocal = None

def init_db(db_url: str):
    global engine, SessionLocal
    engine = create_engine(
        db_url, 
        connect_args={"check_same_thread": False} if "sqlite" in db_url else {},
        pool_pre_ping=True
    )
    Base.metadata.create_all(bind=engine)
    SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
