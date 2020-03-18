'''
seis_model

classes for initial loads and processing of seismic models for use with yt
'''
from netCDF4 import Dataset
import h5py
import os
import numpy as np
import yt
from scipy import spatial
from . import datamanager as dm
from scipy.spatial import Delaunay

def sphere2cart(phi,theta,radius):
    '''
    sphere2cart

    Parameters
    ----------
    phi: angle from north in radians (0 = north pole)
    theta: longitudinal angle in radians
    radius: radius, in any units

    Returns
    -------
    (x,y,z) : tuple of cartesian x,y,z in units of radus
    '''
    x=radius * np.sin(phi) * np.sin(theta)
    y=radius * np.sin(phi) * np.cos(theta)
    z=radius * np.cos(phi)
    return (x,y,z)

def cart2sphere(x,y,z,geo=True):
    xy = x**2 + y**2
    R = np.sqrt(xy + z**2)
    lat = np.arctan2(np.sqrt(xy), z)*180./np.pi
    lon = np.arctan2(y, x)*180./np.pi
    if geo:
        lat = lat - 90. # equator is at 0, +90 is N pole

    return (R,lat,lon)

def in_convexHull(testPoints, hullPoints):
    """
    in_convexHull()

    checks if testPoints fall within convex hull defined by points in hullPoints
    Test if points in `p` are in `hull`

    Parameters:
    ----------
    testPoints  array of coordinates to check: N x n_dim, N is number of points
                to check, n_dim is number of dimensions
    hullPoints  coordinates defining the convex hull: M x n_dim, M is number of
                points on the hull, n_dim is number of dimensions.

    Output:
    -------
    in_hull    boolean array of length N

    """

    hull = Delaunay(hullPoints) # builds Delaunay triangulation
    in_hull = hull.find_simplex(testPoints)>=0 # find_simplex returns -1 if not in hull
    return in_hull

class netcdf(object):
    '''
    netcdf

    class for working with 3D models from iris
    '''
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

    # def addBoundaryHighlights(self,data_field,bound_vals,data=None,newFill=None):
    #     '''
    #     sets the values of the data on the boundaries to specified values and
    #     returns the modified data array. When volume rendering, add TF gaussians
    #     at these boundary values to highlight the domain extents.
    #
    #     Parameters
    #     ----------
    #     data_field   the data field to set and return
    #     bound_vals   dict, the boundary values to set. e.g.,:
    #             bound_vals={'x':[9999.,9998],'y':[9997.,9996],'z':[9995.,9994]}
    #     data    (optional) the data array. If None, will pull data_field from
    #             loaded dataset
    #
    #     '''
    #
    #     try:
    #         dims = list(self.data.variables[data_field].dimensions)
    #         if data is None:
    #             data = self.data.variables[data_field][:]
    #     except:
    #         raise ValueError('netcdf object does not have data loaded')
    #
    #     if newFill is not None:
    #         data.fill_value=newFill
    #
    #     for dim in dims:
    #         print(dim)
    #         if dim in bound_vals.keys():
    #             vals=bound_vals[dim]
    #             print(vals)
    #             this_dim=self.data.variables[dim][:]
    #             data[this_dim==this_dim.min()]=vals[0]
    #             data[this_dim==this_dim.max()]=vals[1]
    #     return data

    def moveToUnstructured(self,fields=[],coordfields=[]):
        '''
        moves regular geo-spherical mesh to unstructured mesh in x,y,z

        (x,y,z)=(0,0,0) is planet center

        data should be depth, latitude, longitude coords
        '''

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
        '''
        checkDomain()

        checks if cartesian coordinates fall within the original domain

        Parameters:
        -----------
        xdata,ydata,zdata:   1-D arrays of length N where N is number of points

        Output:
        -------
        in_hull:             1-D boolean of length N, True if in domain.
        '''

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
        '''
        seis_model.netcdf4.interp2cartesian()
        
        moves geo-spherical data (radius/depth, lat, lon) to earth-centered
            cartesian coordinates using a kdtree

        Parameters:
        ----------
        fields          the fields to interpolate to the grid
        res             list of resolution in x,y,z
        input_units     the units of res (final coord system will be in these units)
        max_dist        the max distance away for nearest neighbor search
        store_trees     if True, will store the kdtree(s) generated
        interpChunk     the chunk size for querying the kdtree.

        Output:
        ------
        adds an 'interp' attribute to netcdf:
        netcdf.interp['grid']={'x':X,'y':Y,'z':Z}   dict of x,y,z 1D arrays (not meshgrid)
        netcdf.interp['data'][field]             3D np array of shape (Nx,Ny,Nz)
                                                 for each field in fields
        '''

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

    def interpFilename(self,field,res=[10,10,10], input_units='km',max_dist=100):
        '''
        interpFilename()
        returns the filename for an interpolated file

        Parameters
        ----------
        **kwargs    the kwargs for for interp2cartesian()

        Output
        ------
        fname       the full path filename of the interpolated file, whether
                    or not it exists
        '''

        db_dir=os.path.join(self.db.db_path,'interpolated_models')

        fullfname=getattr(self,'fname','')
        model=os.path.splitext(os.path.basename(fullfname))[0]
        fname='_'.join([model,field,str(max_dist),str(res[0]),str(res[1]),str(res[2])])
        fname=fname+'.h5py'
        return os.path.join(db_dir,fname)

    def parseInterpFilename(self,fname):
        '''
        parseInterpFilename()
        parses the filename of an interpolated file and returns the interpolation
        settings

        Parameters:
        -----------
        fname  the filename to parse

        Output:
        -------
        interp_settings  dict with interpolation settings

        '''
        basename=os.path.splitext(os.path.basename(fname))[0]
        settings=basename.split('_')
        setfields=['model','field','max_distance','res_x','res_y','res_z']
        return dict(zip(setfields,settings))

    def loadInterpolated(self,field='dvs',**kwargs):
        '''
        loadInterpolated()
        loads interpolated model file data

        Parameters
        ----------
        field   the field to load/interpolate, string
        **kwargs
            the kwargs for for interp2cart: res, input_units, max_dist
            generate    if True, will build the interpolation and save it if not found.
        '''

        res=kwargs.get('res',[10,10,10])
        input_units=kwargs.get('input_units','km')
        max_dist=kwargs.get('max_dist',100)
        chunk=kwargs.get('interpChunk',500000)

        fname=self.interpFilename(field,res=res, input_units=input_units,max_dist=max_dist)
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
