"""
This class generates the database and create connection
"""
import os

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

path = os.path.dirname(os.path.abspath(__file__))

db_path = os.path.join(path, 'database.db')

# Create Engine if doesn't exists
engine = create_engine(f'sqlite:///{db_path}')

# a sessionmaker(), also in the same scope as the engine
Session = sessionmaker(engine)




