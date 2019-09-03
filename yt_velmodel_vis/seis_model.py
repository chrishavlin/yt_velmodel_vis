'''
seis_model

classes for initial loads and processing of seismic models for use with yt
'''
from netCDF4 import Dataset
import os

class netcdf(object):
    '''
    netcdf

    class for working with 3D models from iris
    '''
    def __init__(self,fname=None,load_ds=True):
        self.fname=fname
        self.model_dir=os.environ['YTVELMODELDIR']
        if fname is not None:
            self.validateFile(fname)
            self.load()

        return

    def validateFile(self,fname):
        '''
        validateFile

        checks that model file exists, searches model directory if not found
        initially.

        Parameters
        ----------
        fname : the filename provided (required)

        Returns
        -------
        validFile : boolean for existence. 
        '''
        validFile=False
        if os.path.isfile(fname):
            validFile=True
        else:
            if self.model_dir is not None:
                print("file not found, searching model directory, "+self.model_dir)
                for root, dirs, files in os.walk(self.model_dir):
                    if fname in files:
                        self.fname=os.path.join(root,fname)
                        validFile=True
                if validFile:
                    print("found matching model file:"+self.fname)
                else:
                    print("could not find model file "+fname)

        return validFile


    def load(self,fname=None):
        '''
        load

        loads netcdf file, stores useful attributes

        Parameters
        ----------
        fname if None (default), will check attributes for fname

        Returns
        -------
        Loads netcdf dataset into self.data. Other useful attributes are:
        variables: list of variables, list of str
        coords: list of coordinates, list of str
        ndim: number of dimensions, integer
        varshape: the shape of the variables, tuple of len 3
        '''
        if fname is None:
            try:
                fname=self.fname
            except:
                print("fname required for seis_model.netcdf.load()")

        self.data=Dataset(fname)
        self.variables=list(self.data.variables.keys())
        self.coords=list(self.data.dimensions.keys())
        self.ndim=len(self.coords)
        if len(self.variables)>0:
            v=self.variables[0]
            self.varshape=self.data[v].shape

        return
