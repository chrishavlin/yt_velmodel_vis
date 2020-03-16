'''
This script sets up your local filesystem database required to run the
example scripts.
'''
from yt_velmodel_vis import datamanager as dm
import argparse

parser = argparse.ArgumentParser(description='initialize the filesystem database')
parser.add_argument('-top_dir',type=str,help='the top level directory for database',required=True)
arg = parser.parse_args()
D=dm.initializeDB(top_level_dir=arg.top_dir)
