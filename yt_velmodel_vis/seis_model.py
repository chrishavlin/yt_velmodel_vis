'''
seis_model

classes for initial loads and processing of seismic models for use with yt
'''
from netCDF4 import Dataset
import os
import numpy as np
import yt

def sphere2cart(phi,theta,radius):
    '''
    sphere2cart

    Parameters
    ----------
    phi: angle from north in radians (0 = north pole)
    theta: latitudinal angle in radians
    radius: radius, in any units

    Returns
    -------
    (x,y,z) : tuple of cartesian x,y,z in units of radus
    '''
    x=radius * np.sin(phi) * np.sin(theta)
    y=radius * np.sin(phi) * np.cos(theta)
    z=radius * np.cos(phi)
    return (x,y,z)

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

        self.setUsefulInfo()
        self.coordTransform()

        return

    def setUsefulInfo(self):
        '''
        setUsefulInfo()

        copies netcdf characteristics to attributes for easy access

        '''

        self.coords={}
        for key, value in self.data.dimensions.items():
            self.coords[key]=value.size
        self.ndim=len(self.coords)

        self.varinfo={}
        for key,value in self.data.variables.items():
            if key not in list(self.coords.keys()):
                self.varinfo[key]={'shape':value.shape}

        return

    def coordTransform(self,method='geo2spherical'):
        '''
        coordTransform

        some coordinate transforms
        '''

        if method=='geo2spherical':
            # convert lat/lon/depth to spherical coordinates:
            # phi = azimuth (longitudinal), 0 <= phi <= 2 pi
            # theta = angle from pole (inclination), 0 <= theta <= pi
            #Note that in yt, "spherical" means that it is ordered $r$, $\theta$, $\phi$,
            # where $\theta$ is the declination from the azimuth (running from $0$ to $\pi$
            # and $\phi$ is the angle around the zenith (running from $0$ to $2\pi$).
            lon=np.array(self.data.variables['longitude'])
            lon[lon<0.]=lon[lon<0]+360.

            if 'depth' in self.coords.keys():
                Radius = 6371 -np.array(self.data.variables['depth'])
            elif 'radius' in self.coords.keys():
                Radius = np.array(self.data.variables['radius'])
            else:
                print("depth/radius coordinate not found!")

            self.spherical={
                'theta':(90. -np.array(self.data.variables['latitude']))*np.pi/180.,
                'phi':lon*np.pi/180.,
                'R':Radius / 6371.
            }

            self.spherical['bounds']={}
            for crd in ['phi','theta','R']:
                self.spherical['bounds'][crd]=[self.spherical[crd].min(),
                                               self.spherical[crd].max()]

            self.spherical['bbox']=np.array([
                self.spherical['bounds']['R'],
                self.spherical['bounds']['theta'],
                self.spherical['bounds']['phi']
            ])

        elif method=='sphere2cart':
            # create 3d meshgrid, preserving the variable dimension ordering
            var=list(self.varinfo.keys())[0]
            ordr=self.data.variables[var].dimensions

            newdims=np.meshgrid(np.array(self.data.variables[ordr[0]]),
                                np.array(self.data.variables[ordr[1]]),
                                np.array(self.data.variables[ordr[2]]),
                                indexing='ij')
            crds={}
            for idim in range(0,3):
                crds[ordr[idim]]=newdims[idim]
            crds['lat']=crds['latitude']
            crds['lon']=crds['longitude']

            # assume lat/lon points lie on perfect sphere
            crds['depth'] =6371. - crds['depth'] # depth is now Radius
            crds['lat'] = ( 90. - crds['lat'] ) * np.pi / 180. # lat is now deg from North
            crds['lon'][crds['lon']<0.0]=crds['lon'][crds['lon']<0.0]+360.
            crds['lon'] = crds['lon'] * np.pi / 180.

            X,Y,Z=sphere2cart(crds['lat'],crds['lon'],crds['depth']*1000.)
            self.cart={'X':X,'Y':Y,'Z':Z}

            self.cart['bounds']={}
            for crd in ['X','Y','Z']:
                self.cart['bounds'][crd]=[self.cart[crd].min(),
                                               self.cart[crd].max()]

            self.cart['bbox']=np.array([
                self.cart['bounds']['X'],
                self.cart['bounds']['Y'],
                self.cart['bounds']['Z']
            ])
        else:
            print(method + " is not a defined method")

        return

    def generateDepthSlices(self,dfield,depths=[],ds=None,latlonRnge=None,
        outDir=None,bbox=None,SavePrefix='',data_min=None,data_max=None,fillval=np.nan):
        '''generateDepthSlices()
        generates depth slices naively using uniform data loader
        
        '''
        if latlonRnge is None:
            bnds={
            'longitude':[float(self.data.geospatial_lon_min),
                        float(self.data.geospatial_lon_max)],
            'latitude':[float(self.data.geospatial_lat_min),
                        float(self.data.geospatial_lat_max)],
            'depth':[float(self.data.geospatial_vertical_min),
                    float(self.data.geospatial_vertical_max)]
            }

            for bnd in [0,1]:
                if bnds['latitude'][bnd] > 180:
                   bnds['latitude'][bnd]=bnds['latitude'][bnd]-360
        else:
            bnds=latlonRnge

        clalon=(np.mean(bnds['latitude']),np.mean(bnds['longitude']))
        data={}
        data[dfield]=self.data.variables[dfield][:]
        data[dfield] = np.transpose(data[dfield], (2, 1, 0))
        data[dfield][data[dfield]==self.data[dfield].missing_value]=fillval

        if data_min is not None:
            data[dfield][data[dfield]<data_min]=fillval

        if data_max is not None:
            data[dfield][data[dfield]>data_max]=fillval

        if ds is None:
            sc_mult=1.0 # scale multiplier
            bbox = np.array([bnds['longitude'],bnds['latitude'],bnds['depth']])
            ds = yt.load_uniform_grid(data,data[dfield].shape,sc_mult,bbox=bbox,nprocs=1,periodicity=(False,False,False))

        clon=np.mean(bnds['longitude'])
        clat=np.mean(bnds['latitude'])
        for d_slice in depths:
            cnt=[clon,clat,d_slice]
            slc=yt.SlicePlot(ds,'z',dfield,center=cnt)
            slc.set_log(dfield, False)
            slc.set_ylabel("lat - center lat")
            slc.set_xlabel("lon - center lon")
            slc.set_cmap(dfield, 'gist_heat')
            slc.annotate_title("Depth = "+str(int(d_slice))+' km, C='+str(clalon))
            if outDir is not None:
               slc.save(os.path.join(outDir,SavePrefix+'_slice_'+str(int(d_slice))))

        return
