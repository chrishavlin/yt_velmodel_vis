import os
try:
    import urllib.request as urlrequest
except ImportError:
    import urllib as urlrequest
import time
import pandas as pd
from zipfile import ZipFile
from yt_velmodel_vis import datafiles
from distutils.dir_util import copy_tree

global default_db_path
default_db_path=os.path.join(os.path.expanduser("~"),'.ytvelmodeldata')

class filesysDB(object):
    def __init__(self,top_level_dir=None):
        """
        simple file system database manager

        Parameters
        ----------
        top_level_dir   the top level directory to search for data files. If
                        None, will check for environment variable YTVELMODELDIR
        """
        if top_level_dir is None:
            top_level_dir = os.environ.get('YTVELMODELDIR')

        if top_level_dir is None:
            top_level_dir=default_db_path

        self.db_path=top_level_dir
        self.buildFileDictionary()
        return

    def buildFileDictionary(self):
        """builds dictionaries of available files

        Returns
        -------
        None


        """
        self.FileDict = {}
        self.FilesByDir = {
            'IRIS_models':[],
            'IRIS_refModels':[],
            'shapedata':[]
        }
        for root, subdirs, files in os.walk(self.db_path):
            fileList=files # files in this root
            for file in files:
                full_fi=os.path.join(root,file)
                self.FileDict[file]=full_fi
                for dirname in self.FilesByDir.keys():
                    if dirname in full_fi:
                        self.FilesByDir[dirname].append(file)
        return

    def validateFile(self,fname):
        """ checks if file exists

        Returns
        -------
        str or bool
            returns the filename if it exists, False if not
        """
        validFile=False
        if os.path.isfile(fname):
            validFile=fname
        else:
            if fname in self.FileDict.keys():
                validFile=self.FileDict[fname]

        return validFile

class initializeDB(filesysDB):
    """initializes the filestyem database (cache)

    builds directory structure and fetches file. Will not fetch files that are
    already within the directories.

    Parameters
    ----------
    **kwargs : keyword dict with following possible parameters

    top_level_dir : str
        the top level directory of filesystem db cache. If not provided, will
        check if the YTVELMODELDIR environment variable is set

    build : bool
        if True (default), will build the db and fetch sample files.

    """
    def __init__(self,**kwargs):
        lb=os.linesep
        print(lb+"Initializing filesystem data cache...")
        # find the top level directory to initialize at
        top_level_dir=kwargs.get('top_level_dir',None)
        if top_level_dir is None:
            top_level_dir = os.environ.get('YTVELMODELDIR')

        if top_level_dir is None:
            top_level_dir=default_db_path

        if top_level_dir is not None and os.path.isdir(top_level_dir) is False:
            try:
                os.mkdir(top_level_dir)
            except:
                raise ValueError('Could not create top level directory ' +top_level_dir)

        filesysDB.__init__(self,top_level_dir)
        self.url_pause=3
        if kwargs.get('build',True):
            self.buildDb() # builds the initial directory framework
            self.fetchIRIS() # fetches iris data (will not re-copy)
            self.cpShp() # copies shapefile data
            self.fetchNatEarth() # fetches natural earth shapefiles

            print(lb+'Filesystem database initialization complete.')
            if os.environ.get('YTVELMODELDIR') is None:
                print(('Please set the environment variable YTVELMODELDIR'
                       ' to '+self.db_path+' . For unix systems with bash,'
                       ' add the following to your .bashrc or .bash_aliases file:'
                       +lb+lb+'    export YTVELMODELDIR='+self.db_path+lb
                       ))
        return

    def buildDb(self):
        """ builds the top level subdirectories """
        # build the subdirectories
        print(os.linesep+"Building directory structure at "+self.db_path)
        for db_dir in ['IRIS_models','shapedata','IRIS_refModels','interpolated_models','output']:
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
        """fetches some netcdf files from IRIS, saves in local db cache

        Returns
        -------
        None
        """

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
                time.sleep(self.url_pause)
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
                time.sleep(self.url_pause)
            else:
                print('    '+model+' already downloaded.')
        return

    def cpShp(self):
        """copies shapefile data from package to local db cache

        Returns
        -------
        None
        """
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

    def fetchNatEarth(self):
        """attempts to fetch and unpack useful shapefiles from
        https://www.naturalearthdata.com

        Returns
        -------
        None

        """
        print(os.linesep+"Fetching remote shapefiles")
        # build list of URLs to fetch
        urlbase='https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/'

        resolutions=['10m','50m','110m']
        FiList=[
            ['cultural/ne_','_admin_0_countries.zip'],
            ['cultural/ne_','_admin_1_states_provinces.zip'],
            ['physical/ne_','_coastline.zip'],
        ]

        full_urls=[]
        for res in resolutions:
            for fi in FiList:
                thisfi=res+'/'+fi[0]+res+fi[1]
                full_urls.append(urlbase+thisfi)

        # try to fetch each one
        shapeDest=os.path.join(self.db_path,'shapedata')
        for fi in full_urls:
            zipname=fi.split('/')[-1]
            folderName=zipname.split('.')[0]
            destDir=os.path.join(shapeDest,folderName)
            if os.path.isdir(destDir) is False:
                # unpacked folder doesnt exist
                zip_path=os.path.join(shapeDest,zipname)
                if os.path.isfile(zip_path) is False:
                    # the zip file doesnt exist, go fetch it
                    print('    fetching ' +zipname)
                    try:
                        urlrequest.urlretrieve(fi, zip_path)
                    except:
                        msg=('    WARNING: Could not fetch '+fi)
                        print(msg)
                    time.sleep(self.url_pause)
                else:
                    print('    '+zipname+' already downloaded')

                removeTheZip=True
                try:
                    # unpack the zip
                    print('    unpacking ' +zipname)
                    with ZipFile(zip_path, 'r') as zip_ref:
                        zip_ref.extractall(destDir)
                except:
                    print('    WARNING: Could not unpack '+zip_path)
                    removeTheZip=False

                if removeTheZip:
                    os.remove(zip_path)


        return
