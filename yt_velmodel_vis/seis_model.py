"""
seis_model

classes for initial loads and processing of seismic models for use with yt
"""
from netCDF4 import Dataset
import h5py
import os
import numpy as np
import yt
from scipy import spatial
from scipy.interpolate import interp1d
from . import datamanager as dm
from scipy.spatial import Delaunay
import pandas as pd

def sphere2cart(phi,theta,radius):
    """
    seis_model.sphere2cart(phi,theta,radius)

    transformation from yt spherical coordinates to cartesian

    Parameters
    ----------
    phi : ndarray or scalar float/ing
        angle from north in radians (0 = north pole)
    theta : ndarray or scalar float/ing
        longitudinal angle in radians
    radius : ndarray or scalar float/ing
        radius in any units

    all arrays must be the same size (or 2 of 3 can be scalars)

    Returns
    -------
    (x,y,z) : tuple of cartesian x,y,z in same units as radius
    """
    x=radius * np.sin(phi) * np.sin(theta)
    y=radius * np.sin(phi) * np.cos(theta)
    z=radius * np.cos(phi)
    return (x,y,z)

def cart2sphere(x,y,z,geo=True):
    """
    seis_model.cart2sphere(x,y,z,geo=True)

    transformation from cartesian to spherical coordinates

    Parameters
    ----------
    x, y, z   cartesian coordinate arrays
    geo       boolean, if True then latitude is 0 at equator, otherwise 0 at
              the north pole.

    all arrays must be the same size (or 2 of 3 can be scalars)

    Returns
    -------
    (R,lat,lon) : tuple of cartesian radius, lat and lon (lat,lon in degrees)

    """
    
    xy = x**2 + y**2
    R = np.sqrt(xy + z**2)
    lat = np.arctan2(np.sqrt(xy), z)*180./np.pi
    lon = np.arctan2(y, x)*180./np.pi
    if geo:
        lat = lat - 90. # equator is at 0, +90 is N pole

    return (R,lat,lon)

class netcdf(object):
    """
    netcdf

    class for working with 3D models from iris
    """
    def __init__(self,fname=None,interpDict={}):

        self.db=dm.filesysDB()
        fullfilename=self.db.validateFile(fname)
        if fullfilename is not False:
            self.fname=fullfilename
            self.load()
        else:
            print('Could not find file '+fname)

        interpfield=interpDict.get('field',None)
        if interpfield is not None:
            self.loadInterpolated(**interpDict)

        return


    def load(self,fname=None):
        """
        load

        loads netcdf file, stores useful attributes

        Parameters
        ----------
        fname : str
            filename to load, if None (default), will check attributes for fname

        Returns
        -------
        Loads netcdf dataset into self.data. Other useful attributes are:
            variables: list of variables, list of str
            coords: list of coordinates, list of str
            ndim: number of dimensions, integer
            varshape: the shape of the variables, tuple of len 3
        """
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
        """ copies netcdf characteristics to attributes for easy access """

        self.coords={}
        for key, value in self.data.dimensions.items():
            self.coords[key]=value.size
        self.ndim=len(self.coords)

        self.varinfo={}
        for key,value in self.data.variables.items():
            if key not in list(self.coords.keys()):
                self.varinfo[key]={'shape':value.shape}

        self.fields=list(self.varinfo.keys())
        return

    def coordTransform(self,method='geo2spherical'):
        """
        seis_model.netcdf.coordTransform()

        some coordinate transformations

        Parameters
        ----------
        method : str
            the method to use, must be:

            'geo2spherical': converts lat/lon/radius to spherical coordinates
                             conventions used by yt (i.e., N. pole is 0, not +90)
            'sphere2cart': conerts spherical coordinates to cartesian
                           coordinates.

        Returns
        ------
        Sets attributes depending on method.
        For geo2spherical, sets:
            netcdf.spherical :
                a dictionary with coordinate arrays and bounding box info
        For sphere2cart, sets:
            netcdf.cart :
                a dictionary with coordinate arrays and bounding box info
        """

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

    def moveToUnstructured(self,fields=[],coordfields=[]):
        """ moves regular geo-spherical mesh to unstructured mesh in x,y,z

        (x,y,z)=(0,0,0) is planet center

        data should be depth, latitude, longitude coords

        (not worth it, use interp2cartesian())
        """

        uData={}

        # calculate x,y,z from regular geo-spherical data
        self.coordTransform('sphere2cart') # builds self.cart={'X':X,'Y':Y,'Z':Z}
        print("coord transformed!")

        # initialize data value store
        if len(fields)>0:
            data1d={} # data arrays in 1d
            uData['data']={} # data values on vertices for each element
            for fi in fields:
                print("copying "+fi)
                uData['data']['connect1',fi]=[]
                data1d[fi] = self.data.variables[fi][:].data.ravel() # 1d data array

        # add spherical coordinate fields
        if len(coordfields)>0:
            if 'radius' in coordfields:
                uData['data']['connect1','radius']=[]#np.sqrt(X**2+Y**2+Z**2)
            elif 'depth' in coordfields:
                uData['data']['connect1','depth']=[]#6371000. - np.sqrt(X**2+Y**2+Z**2)

        # build the array of all vertices. cartesian is already 3d
        print("building X,Y,Z")
        X=self.cart['X'].ravel()
        Y=self.cart['Y'].ravel()
        Z=self.cart['Z'].ravel()

        uData['vertices']=[[X[i],Y[i],Z[i]] for i in range(0,len(X))]
        uData['vertices']=np.array(uData['vertices'])
        uData['vertices']=np.float64(uData['vertices'])

        # build the connectivity array for uniform spherical volume elements
        # X,Y,Z=sphere2cart(crds['lat'],crds['lon'],crds['depth']*1000.)
        Nx = len(self.cart['X'][0][0]) # verts in x
        Ny = len(self.cart['X'][0]) # verts in y
        Nz = len(self.cart['X']) # verts in z
        N = Nx*Nz*Ny # number of verts
        print("total verts [x,y,z,total,len()]: ")
        print([Nx,Ny,Nz,N,len(X)])
        Neles=(Nx-1)*(Ny-1)*(Nz-1) # number of elements
        print("total elements: "+str(Neles))
        M = 8 # vertices per element
        uData['cnnct']=[] # connnectivity array (wihch vertex for each element)
        for ele in range(0,Neles):

            verts = [
                ele, ele + 1,
                ele + Nx, ele + Nx + 1,
                ele + Nx * Ny, ele + Ny*Nx + 1,
                ele + Nx * Ny + Nx, ele + Ny*Nx + Nx + 1
            ]

            uData['cnnct'].append(verts) # indeces of vertices

            if len(fields)>0:
                for fi in fields:
                    uData['data']['connect1',fi].append(data1d[fi][verts])

            xverts=X[verts]
            yverts=Y[verts]
            zverts=Z[verts]
            if 'radius' in coordfields:
                uData['data']['connect1','radius'].append(np.sqrt(xverts**2+yverts**2+zverts**2))
            elif 'depth' in coordfields:
                uData['data']['connect1','depth'].append(6371000.-np.sqrt(xverts**2+yverts**2+zverts**2))


        # convert final to np arrays
        uData['cnnct']=np.array(uData['cnnct'])
        for fi in fields + coordfields:
            uData['data']['connect1',fi]=np.array(uData['data']['connect1',fi])


        setattr(self,'unstructured',uData)
        return

    def checkDomain(self,xdata,ydata,zdata):
        """ checks if cartesian coordinates fall within spherical domain

        Builds a convex hull of the original spherical domain and checks if the
        supplied coordinates fall within the domain.

        Generally not worth the overhead.

        Parameters
        -----------
        xdata,ydata,zdata:  ndarrays
            1-D arrays of length N where N is number of points

        Returns
        -------
        in_hull: ndarray
            1-D boolean array of length N, True if in domain.
        """

        # build hull framework from lat/lon/depth extents
        hull_coords=[]

        # pull out the lat/lon/depth arrays
        lat=np.array(self.data.variables['latitude'])
        lon=np.array(self.data.variables['longitude'])

        try:
            radius=np.array(self.data.variables['depth'])
            radius=6371. - radius
        except:
            radius=np.array(self.data.variables['radius'])

        # units and angle restrictions
        radius=radius*1000
        lat=( 90. - lat ) * np.pi / 180.
        lon[lon<0]=lon[lon<0]+360.
        lon = lon * np.pi / 180.

        # min/max radius "planes"
        lats,lons=np.meshgrid(lat,lon)
        for rad in [radius.min(),radius.max()]:
            x,y,z=sphere2cart(lats,lons,np.full(lats.shape,rad))
            hull_coords.append(np.column_stack((x.ravel(),y.ravel(),z.ravel())))

        # min/max lat "planes"
        lons,rads=np.meshgrid(lon,radius)
        for lat0 in [lat.min(),lat.max()]:
            x,y,z=sphere2cart(np.full(lons.shape,lat0),lons,rads)
            hull_coords.append(np.column_stack((x.ravel(),y.ravel(),z.ravel())))

        # min/max lon "planes"
        lats,rads=np.meshgrid(lat,radius)
        for lon0 in [lon.min(),lon.max()]:
            x,y,z=sphere2cart(lats,np.full(lats.shape,lon0),rads)
            hull_coords.append(np.column_stack((x.ravel(),y.ravel(),z.ravel())))

        hull_coords=np.concatenate(hull_coords)

        return in_convexHull(np.column_stack((xdata,ydata,zdata)),hull_coords)

    def interp2cartesian(self,fields=[],res=[10,10,10], input_units='km',
                        max_dist=100,store_trees=False,interpChunk=500000):
        """
        seis_model.netcdf4.interp2cartesian()

        moves geo-spherical data (radius/depth, lat, lon) to earth-centered
        cartesian coordinates using a kdtree

        Parameters
        ----------
        fields : list of str
            the fields to interpolate to the grid (default [])
        res: list of int
            list of resolution in x,y,z (default [10,10,10])
        input_units : str
            the units of res and max_dist (default 'km'). The final coord system
            will be in these units.
        max_dist : int or float
            the max distance away for nearest neighbor search (default 100)
        store_trees : boolean
            if True, will store the kdtree(s) generated (default False)
        interpChunk : int
            the chunk size for querying the kdtree (default 500000)

        Returns
        -------
        adds an 'interp' attribute to netcdf:
        netcdf.interp['grid'] : dict of ndarrays
            {'x':X,'y':Y,'z':Z}, dict of x,y,z 1D arrays (not meshgrid)
        netcdf.interp['data'][field] : ndarray
            3D np array of shape (Nx,Ny,Nz) for each field in fields
        """

        self.interp={}

        if len(fields)==0:
            print("no fields provided, building grid only")

        if getattr(self,'cart',None) is None:
            # get cartesian coordinates of the data (no interpolation here)
            self.coordTransform('sphere2cart') # builds self.cart={'X':X,'Y':Y,'Z':Z}

        # do interpolation in meters, build res dict
        res = np.array(res)
        if input_units=='km':
            res=res * 1000. # coord transform are in m
            max_dist=max_dist*1000.
        res = dict(zip(['X','Y','Z'],list(res)))

        # build dimensional info
        dimInfo={}
        for dim in ['X','Y','Z']:
            dimInfo[dim]={
                1:self.cart['bounds'][dim][1],
                0:self.cart['bounds'][dim][0],
                'delta':abs(self.cart['bounds'][dim][1]-self.cart['bounds'][dim][0])
            }
            dimInfo[dim]['N']=int(dimInfo[dim]['delta']/res[dim])

        # store bounding box here too
        self.interp['bbox']=np.array([
            self.cart['bounds']['X'],
            self.cart['bounds']['Y'],
            self.cart['bounds']['Z']
        ])

        # build the 1d spatial vectors
        X=np.linspace(dimInfo['X'][0],dimInfo['X'][1],dimInfo['X']['N'])
        Y=np.linspace(dimInfo['Y'][0],dimInfo['Y'][1],dimInfo['Y']['N'])
        Z=np.linspace(dimInfo['Z'][0],dimInfo['Z'][1],dimInfo['Z']['N'])
        self.interp['grid']={'x':X,'y':Y,'z':Z}

        # interpolate data fields onto the grid, ignore fill values
        if len(fields)>0:
            # self.cart['X']={'X':X,'Y':Y,'Z':Z} The cartesian coords of data are here
            xdata=self.cart['X'].ravel()
            ydata=self.cart['Y'].ravel()
            zdata=self.cart['Z'].ravel()

            self.interp['data']={}
            trees={}
            for fi in fields:
                # build kdtree of non-null data for each variable
                # (might not be OK to assume same location for non-null in each field)

                if fi == 'dvs':
                    data=self.check_dvs()
                    data=data.ravel()
                    fillval=np.nan
                else:
                    fillval=getattr(self.data.variables[fi],'missing_value',np.nan)
                    data=self.data.variables[fi][:].data.ravel()

                x_fi=xdata[data!=fillval]
                y_fi=ydata[data!=fillval]
                z_fi=zdata[data!=fillval]
                data=data[data!=fillval]
                xyz=np.column_stack((x_fi,y_fi,z_fi))
                print("building kd tree for "+fi)
                trees[fi]={'tree':spatial.cKDTree(xyz),'data':data}
                print("    kd tree built")

                # initialize interpolated field
                self.interp['data'][fi]=np.full((dimInfo['X']['N'],dimInfo['Y']['N'],dimInfo['Z']['N']),np.nan)

            # build the cartesian grid
            print('building cartesian grid')
            xdata,ydata,zdata=np.meshgrid(X,Y,Z,indexing='ij')
            orig_shape=xdata.shape
            xdata=xdata.ravel(order='C')
            ydata=ydata.ravel(order='C')
            zdata=zdata.ravel(order='C')

            # limit to the points only in domain (not worth the overhead)
            #     full_linear_indxs=np.array(range(0,len(xdata)))
            #     in_domain=self.checkDomain(xdata,ydata,zdata)
            #     full_linear_indxs=full_linear_indxs[in_domain]
            #     xdata=xdata[in_domain]
            #     ydata=ydata[in_domain]
            #     zdata=zdata[in_domain]

            # query the tree at each new grid point and weight nearest neighbors
            # by inverse distance. proceed in chunks.
            N_grid=len(xdata)
            print("querying kdtree on interpolated grid")
            chunk=interpChunk
            N_chunks=int(N_grid/chunk)+1
            print("breaking into "+str(N_chunks)+' chunks')
            for i_chunk in range(0,N_chunks):
                print("   processing chunk "+str(i_chunk+1)+" of "+str(N_chunks))
                i_0=i_chunk*chunk
                i_1=i_0+chunk
                if i_1>N_grid:
                    i_1=N_grid
                pts=np.column_stack((xdata[i_0:i_1],ydata[i_0:i_1],zdata[i_0:i_1]))
                indxs=np.array(range(i_0,i_1)) # the linear indeces of this chunk
                for fi in fields:
                    (dists,tree_indxs)=trees[fi]['tree'].query(pts,k=8,distance_upper_bound=max_dist)

                    # remove points with all infs (no NN's within max_dist)
                    m=np.all(~np.isinf(dists),axis=1)
                    tree_indxs=tree_indxs[m]
                    indxs=indxs[m]
                    dists=dists[m]

                    # IDW with array manipulation
                    # Build weighting matrix
                    wts=1/dists
                    wts= wts / np.sum(wts,axis=1)[:,np.newaxis] # shape (N,8)
                    vals=trees[fi]['data'][tree_indxs] # shape (N,8)
                    vals=vals * wts
                    vals=np.sum(vals,axis=1) # shape (N,)

                    # store in proper indeces
                    full_indxs=np.unravel_index(indxs,orig_shape,order='C')
                    self.interp['data'][fi][full_indxs]=vals

            if store_trees:
                self.interp['trees']=trees

        # adjust final grid units
        if input_units=='km':
            self.interp['bbox']=self.interp['bbox']/1000.
            for dim in ['x','y','z']:
                self.interp['grid'][dim]=self.interp['grid'][dim]/1000.



        return

    def interpFilename(self,field,res,max_dist):
        """
        seis_model.netcdf4.interpFilename(field,res,max_dist)

        returns the filename for an interpolated file

        Parameters
        ----------
        field : str
            the data field (e.g., 'dvs')
        res : list of int or float
            resolution list (e.g., [10,10,10])
        mx_dist : int or float
            max distance for NN search (e.g., 100)

        Returns
        ------
        fname : str
            the full path filename of the interpolated file. Does not check
            that the file exists
        """

        db_dir=os.path.join(self.db.db_path,'interpolated_models')

        fullfname=getattr(self,'fname','')
        model=os.path.splitext(os.path.basename(fullfname))[0]
        fname='_'.join([model,field,str(max_dist),str(res[0]),str(res[1]),str(res[2])])
        fname=fname+'.h5py'
        return os.path.join(db_dir,fname)

    def parseInterpFilename(self,fname):
        """
        seis_model.netcdf4.parseInterpFilename(fname)

        parses the filename of an interpolated file and returns the interpolation
        settings

        Parameters
        -----------
        fname : str
            the filename to parse

        Returns
        -------
        interp_settings : dict
            dict with interpolation settings

        """
        basename=os.path.splitext(os.path.basename(fname))[0]
        settings=basename.split('_')
        setfields=['model','field','max_distance','res_x','res_y','res_z']
        return dict(zip(setfields,settings))

    def loadInterpolated(self,field='dvs',**kwargs):
        """
        seis_model.netcdf.loadInterpolated()

        loads interpolated model file data. Will generate interpolated data if
        it does not exist (unless kwarg generate=False)

        Parameters
        ----------
        field : str
            the field to load/interpolate (default 'dvs')
        **kwargs : dict
            can include:
                - the kwargs for for interp2cart: res, input_units, max_dist
                - kwargs['generate'] : boolean
                    if True, will build the interpolation and save it if not
                    found.

        Returns
        ------
        adds an 'interp' attribute to netcdf:
        netcdf.interp['grid'] : dict of ndarrays
            {'x':X,'y':Y,'z':Z} dict of x,y,z 1D arrays (not meshgrid)
        netcdf.interp['data'][field] : ndarray
            3D np array of shape (Nx,Ny,Nz) for each field in fields

        """

        res=kwargs.get('res',[10,10,10])
        input_units=kwargs.get('input_units','km')
        max_dist=kwargs.get('max_dist',100)
        chunk=kwargs.get('interpChunk',500000)

        fname=self.interpFilename(field,res,max_dist)
        if os.path.isfile(fname):
            # load it
            if ~hasattr(self,'interp'):
                self.interp={'data':{},'grid':{}}

            with h5py.File(fname, 'r') as hf:
                self.interp['data'][field] = hf[field][:]
                for dim in ['x','y','z']:
                    self.interp['grid'][dim] = hf[dim][:]

            if ~hasattr(self,'cart'):
                self.cart={}

            self.cart['bbox']=[]
            for dim in ['x','y','z']:
                self.cart['bbox'].append([self.interp['grid'][dim].min(),
                                          self.interp['grid'][dim].max()])
            self.cart['bbox']=np.array(self.cart['bbox'])

        elif kwargs.get('generate',True):
            # generate it
            self.interp2cartesian([field],res,input_units,max_dist,interpChunk=chunk)

            # save it
            with h5py.File(fname, 'w') as hf:
                hf.create_dataset(field,  data=self.interp['data'][field],compression='gzip')
                for dim in ['x','y','z']:
                    hf.create_dataset(dim,  data=self.interp['grid'][dim],compression='gzip')

        return

    def check_dvs(self,ref_model='ref_model'):
        """
        seis_model.netcdf.check_dvs()

        checks for shear wave velocity perturbation, calcuates it if not there.

        Parameters
        -----------
        ref_type : str
            the reference model type to use in calculating perturbations:

            'h_average'     ref model is horizontal average of current model at
                            every depth
            'ref_model':    will use a 1D reference model. If the data has 'vsh'
                            and 'vsv' fields, will compare voigt average to
                            AK135F_AVG. If not, will compare 'vsv' or 'vsh' to
                            PREM. (default)

        Returns
        -------
        dvs : ndarray
            the shear wave velocity perturbation,(v_model - v_ref)/v_ref

        """

        if 'dvs' in self.fields:
            return self.returnFilled_np('dvs',np.nan)
        else:
            print(("seis_model.netcdf.check_dvs: no dvs field, attempting to "
                    "calculate it."))
            return self.calc_dvs()

    def calc_dvs(self,ref_model='1d_ref'):
        """ calculates shear wave velocity perturbations from a reference model

        Parameters
        ----------
        ref_model : str
            ref_model type, 'h_average' or '1d_ref' (the default is 'h_average')

            'h_average'     ref model is horizontal average of current model at
                            every depth
            '1d_ref':       will use a 1D reference model. If the data has 'vsh'
                            and 'vsv' fields, will compare voigt average to
                            AK135F_AVG. If not, will compare 'vsv' or 'vsh' to
                            PREM.

        Returns
        -------
        ndarray
            velocity perturbation on model grid.

        """

        v=None
        if 'vsv' in self.fields and 'vsh' in self.fields:
            vsv=self.returnFilled_np('vsv',np.nan)
            vsh=self.returnFilled_np('vsh',np.nan)
            v=simpleVoigt(vsv,vsh)
            ref_type='voigt'
        elif 'vsv' in self.fields:
            v=self.returnFilled_np('vsv',np.nan)
            ref_type='vsv'
        elif 'vsh' in self.fields:
            v=self.returnFilled_np('vsh',np.nan)
            ref_type='vsh'

        if v is not None:
            if ref_model == 'h_average':
                v_ave=np.nanmean(v,axis=1)
                v_ave=np.nanmean(v_ave,axis=1)[:,np.newaxis,np.newaxis]
            elif ref_model=='1d_ref':
                # interpolate a reference model at the depths of the present
                # model to get average velocity vs depth
                if ref_type=='voigt':
                    ref_m=refModel('AK135F_AVG')
                    v_func=ref_m.interp('Vs_kms')
                elif ref_type=='vsv':
                    ref_m=refModel('PREM500')
                    v_func=ref_m.interp('Vsv_kms')
                elif ref_type=='vsh':
                    ref_m=refModel('PREM500')
                    v_func=ref_m.interp('Vsh_kms')
                model_depths=np.array(self.data.variables['depth'][:])
                v_ave=v_func(model_depths)[:,np.newaxis,np.newaxis]

            v_ave[v_ave<1e-8]=np.min(v_ave[v_ave>1e-8]) # interp close to surf can give 0
            return (v - v_ave)/(v_ave)
        else:
            return np.array([])

    def returnFilled_np(self,field,fillval=None):
        """
        seis_model.netcdf.returnFilled_np(field,fillval=np.nan)

        returns a copy of variable data as a numpy array

        Parameters
        ----------
        field :
            the field to return
        fillval :
            optional fill value to use. If None, will return the default fill
            value of the netcdf variable

        Returns
        -------
        npRA : ndarray
            array of desired field

        """
        npRA=np.array(self.data.variables[field][:])
        if fillval is not None:
            npRA[self.data.variables[field][:].mask]=fillval

        return npRA

class refModel(object):
    """
    seis_model.refModel() IRIS refence model reader for plaintext csv files.

    m = seis_model.refModel(ref_model)

    Parameters
    ----------
    ref_model  the reference model: can be a file path or just a filename if
               the model resides in the yt_velmodel_vis filesystem database
    """
    def __init__(self,ref_model):

        self.db=dm.filesysDB()

        if os.path.isfile(ref_model):
            filename=ref_model
        else:
            filename=self.db.validateFile(ref_model)
            if filename is False:
                filename=self.db.validateFile(ref_model+'.csv')
                if filename is False:
                   raise ValueError('Could not find reference model '+fname)

        self.model_file=filename
        self.df=None
        return

    def return_df(self,column_list=None,**kwargs):
        """
        seis_model.refModel.return_df(column_list=None,**kwargs)

        returns reference model as pandas dataframe

        Parameters
        -----------
        column_list     list of columns to return. If None, returns all
        kwargs          keyword dict, passed to pandas.read_csv()
        """

        df=self.df
        if df is None:
            df=pd.read_csv(self.model_file,**kwargs)

        if column_list is None:
            column_list=df.columns
        return df[column_list]

    def interp(self,column,**kwargs):
        """
        seis_model.refModel.interp(column)

        returns a scipy.interpolate.interp1d function for the depth dependence
        of specified column

        Parameters
        -----------
        column : str
            the column to build interpolation function for
        kwargs : dict
            keyword dict, passed to pandas.read_csv()

        Returns
        ------
        interp :
            the interpolation function

        Example:
        --------
        To calculate Vp at 100 and 200 km:
            Vp_vs_depth = seis_model.refModel.interp('Vp_kms')
            Vp=Vp_vs_depth(np.array([100,200]))
        """
        df=self.df
        if df is None:
            df=pd.read_csv(self.model_file,**kwargs)

        # deal with discontinuities: pull columns, ensure depth is increasing
        col_v=df[column].to_numpy()
        depths=df['depth_km'].to_numpy()
        sorted_i=np.argsort(depths)
        depths=depths[sorted_i]
        col_v=col_v[sorted_i]

        # offset disc depths by a small number
        # disc_vals=[]
        eps_off=1e-8
        for i_z in range(0,len(depths)):
            if depths[i_z]==depths[i_z-1]:
                # disc_vals.append(depths[i_z]) # not used... might be useful
                depths[i_z]=depths[i_z]+eps_off

        # build and return the interpolation function
        interp_f=interp1d(depths,col_v)
        return interp_f

def simpleVoigt(vsh,vsv):
    """
    seis_model.simpleVoigt(vsh,vsv)
    voigt average of horizontal and vertical shear wave velocities

        v = 0.5 * (vsh + vsv)

    Parameters
    -----------
    vsh     horizontal shear wave velocity, array or scalar
    vsv     vertical shear wave velocity, array or scalar

    if vsh and vsv are both arrays, must be the same size

    Returns
    -------
    v       voigt average

    """
    return (vsh+vsv)*0.5


def in_convexHull(testPoints, hullPoints):
    """
    seis_model.in_convexHull()

    checks if testPoints fall within convex hull defined by points in hullPoints

    Parameters
    ----------
    testPoints  array of coordinates to check: N x n_dim, N is number of points
                to check, n_dim is number of dimensions
    hullPoints  coordinates defining the convex hull: M x n_dim, M is number of
                points on the hull, n_dim is number of dimensions.

    Returns
    -------
    in_hull    boolean array of length N

    """

    hull = Delaunay(hullPoints) # builds Delaunay triangulation
    in_hull = hull.find_simplex(testPoints)>=0 # find_simplex returns -1 if not in hull
    return in_hull
