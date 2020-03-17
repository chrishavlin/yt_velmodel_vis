import os
try:
    import urllib.request as urlrequest
except ImportError:
    import urllib as urlrequest
import time
import pandas as pd
from yt_velmodel_vis import datafiles
from distutils.dir_util import copy_tree

class filesysDB(object):
    def __init__(self,top_level_dir=None):
        '''
        simple file system database manager

        Parameters
        ----------
        top_level_dir   the top level directory to search for data files. If
                        None, will check for environment variable YTVELMODELDIR
        '''
        if top_level_dir is None:
            top_level_dir = os.environ.get('YTVELMODELDIR')

        if top_level_dir is None:
            raise ValueError(("Could not intialize filesysDB: top_level_dir is None and "
                   "environment variable YTVELMODELDIR is not set."))

        self.db_path=top_level_dir
        self.buildFileDictionary()
        return

    def buildFileDictionary(self):
        self.FileDict = {}
        for root, subdirs, files in os.walk(self.db_path):
            fileList=files # files in this root
            for file in files:
                self.FileDict[file]=os.path.join(root,file)
        return

    def validateFile(self,fname):
        validFile=False
        if os.path.isfile(fname):
            validFile=fname
        else:
            if fname in self.FileDict.keys():
                validFile=self.FileDict[fname]

        return validFile

class initializeDB(filesysDB):
    """ initializes the filestyem database (cache)"""
    def __init__(self,**kwargs):
        lb=os.linesep
        print(lb+"Initializing filesystem data cache...")
        # find the top level directory to initialize at
        top_level_dir=kwargs.get('top_level_dir',None)
        if top_level_dir is None:
            top_level_dir = os.environ.get('YTVELMODELDIR')

        if top_level_dir is None:
            pathmsg=("YTVELMODELDIR environment variable is not set, please "
                "enter a path to install the filesystem cache (or hit enter"
                " to use default): ")
            top_level_dir=input(pathmsg)
            if len(top_level_dir)==0:
                top_level_dir=os.path.join(os.path.expanduser("~"),'.ytvelmodeldata')

        if top_level_dir is not None and os.path.isdir(top_level_dir) is False:
            try:
                os.mkdir(top_level_dir)
            except:
                raise ValueError('Could not create top level directory ' +top_level_dir)

        filesysDB.__init__(self,top_level_dir)

        self.buildDb() # builds the initial directory framework
        self.IRIS_pause_s=3
        self.fetchIRIS() # fetches iris data (will not re-copy)
        self.cpShp() # copies shapefile data

        print(lb+'Filesystem database initialization complete.')
        if os.environ.get('YTVELMODELDIR') is None:
            print(('Please set the environment variable YTVELMODELDIR'
                   ' to '+self.db_path+' . For unix systems with bash,'
                   ' add the following to your .bashrc or .bash_aliases file:'
                   +lb+lb+'    export YTVELMODELDIR='+self.db_path+lb
                   ))

    def buildDb(self):
        """ builds the top level subdirectories """
        # build the subdirectories
        print(os.linesep+"Building directory structure at "+self.db_path)
        for db_dir in ['IRIS_models','shapedata','IRIS_refModels']:
            newdir=os.path.join(self.db_path,db_dir)
            if os.path.isdir(newdir) is True:
                print('    '+db_dir+' already exists')
            else:
                try:
                    os.mkdir(newdir)
                except:
                    msg=('Could not create subdirectories in '+self.db_path+'. Check '
                        'write permissions.')
                    raise ValueError(msg)

    def fetchIRIS(self):
        ''' fetches a number of models from IRIS '''

        # IRIS earth models
        url_base='https://ds.iris.edu/files/products/emc/emc-files/'
        models_to_fetch=['NA07_percent.nc','NWUS11-S_percent.nc',
            'DNA10-S_percent.nc','US.2016.nc']

        # IRIS reference models
        # reference models are csv files with no header. Following dictionaries
        # set the columns for each reference model and include multiplication
        # factors to normalize units between models.
        Rmodels={
            'AK135F/AK135F_AVG.csv':
                {'cols':['depth_km','rho_kgm3','Vp_kms','Vs_kms','Qkappa','Qmu'],
                 'factors':[1.,1000.,1.,1.,1.,1.]},
            'PREM500/PREM500.csv':{
                 'cols':['radius_km','rho_kgm3','Vpv_kms','Vsv_kms','Qkappa','Qmu',
                                   'Vph_kms','Vsh_kms','eta'],
                 'factors':[0.001,1.,0.001,0.001,1.,1.,0.001,0.001,1.]
                                   }
                }
        url_base_R='http://ds.iris.edu/files/products/emc/data/'

        # fetch the earth models
        print(os.linesep+"Fetching models from IRIS")
        for model in models_to_fetch:
            print("    fetching "+model+" from IRIS...")
            full_url=url_base+model
            local_path=os.path.join(self.db_path,'IRIS_models',model)
            if os.path.isfile(local_path) is False:
                try:
                    urlrequest.urlretrieve(full_url, local_path)
                except:
                    msg=('    Could not fetch '+full_url)
                    print(msg)
                time.sleep(self.IRIS_pause_s)
            else:
                print('    '+model+' already downloaded.')


        # fetch and normalize the reference models
        for model in Rmodels.keys():
            full_url=url_base_R+model
            print("    fetching "+model+" from IRIS...")
            model_short=model.split('/')[-1]
            local_path=os.path.join(self.db_path,'IRIS_refModels',model_short)
            if os.path.isfile(local_path) is False:
                try:
                    urlrequest.urlretrieve(full_url, local_path)
                except:
                    msg=('    Could not fetch '+full_url)
                    print(msg)

                df=pd.read_csv(local_path,header=None,names=Rmodels[model]['cols'])

                colfact=dict(zip(Rmodels[model]['cols'],Rmodels[model]['factors']))
                for col in Rmodels[model]['cols']:
                    df[col]=df[col]*colfact[col]

                if 'depth_km' not in Rmodels[model]['cols'] and 'radius_km' in Rmodels[model]['cols']:
                    df['depth_km']=6371.-df['radius_km']

                df.to_csv(local_path,index=False)
                time.sleep(self.IRIS_pause_s)
            else:
                print('    '+model+' already downloaded.')
        return

    def cpShp(self):
        ''' copies shapefile data from package to local db cache '''
        dataDir=datafiles.__path__._path[0]
        shapeDest=os.path.join(self.db_path,'shapedata')
        print(os.linesep+'Copying packaged shapefile data')
        for shpDir in os.listdir(dataDir):
            fullSrc=os.path.join(dataDir,shpDir)
            fullDest=os.path.join(shapeDest,shpDir)
            if os.path.isdir(fullSrc):
                print('    copying '+shpDir)
                copy_tree(fullSrc,fullDest,verbose=0)

        return
