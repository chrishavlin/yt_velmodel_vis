'''
shapeplotter module: support for plotting shapefile data on yt scenes
'''

from . import seis_model as sm
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import yt
import os
import geopandas as gpd
from . import datamanager as dm

class shapeTrace(object):
    '''
    individual trace of a shape.

    shapeTrace(lat,lon,radius) with lat,lon in degrees.

    lat, lon can be scalar int or float, lists or arrays of the same size and
    shape, radius can be single value or list/array of the same shape
    '''
    def __init__(self,lat,lon,radius=6371.):

        coord_input={'lat':lat,'lon':lon,'radius':radius}
        self.coords=self.processCoordInputs(coord_input)
        self.projectCartesian()

    def processCoordInputs(self,coord_input):
        '''
        initialization function, checks coordinate input after initial
        processing.

        Parameters:
        ----------
        coord_input     dictionary of 'lat', 'lon', 'radius' arrays

        Output:
        ------
        coord_input     modified dictionary of 'lat', 'lon', 'radius' arrays
                        Projects radius into a np array of same size as lat/lon
                        arrays if input radius is a single value.

        Raises errors for failed criteria
        '''

        # make everything an array
        for coord, coord_vals in coord_input.items():
            if type(coord_vals)==int or type(coord_vals)==float:
                coord_input[coord]=np.array([coord_vals])
            if type(coord_vals)==list:
                coord_input[coord]=np.array([coord_vals])

        # error checks (everything should be np array by now)
        if coord_input['lat'].shape!= coord_input['lon'].shape:
            msg=('You must supply lat, lon values of equal length.'
            ' lat, lon values can be scalars, lists or arrays but must be the '
            ' length (and shape if array/list).')
            raise ValueError(msg)

        # check radius
        if coord_input['lat'].shape != coord_input['radius'].shape:
            if coord_input['radius'].shape==(1,):
                R=coord_input['radius'][0]
                coord_input['radius']= np.full(coord_input['lat'].shape, R)
            else:
                msg=('Radius must be a single value or have the same length and '
                ' shape as lat, lon ')
                raise ValueError(msg)

        # lon in (0,360)
        coord_input['lon'][coord_input['lon']<0]=coord_input['lon'][coord_input['lon']<0]+360.

        return coord_input

    def projectCartesian(self):
        '''
        projects coordinates to x,y,z
        '''
        phi=( 90.-self.coords['lat'])*np.pi/180.
        theta=self.coords['lon']*np.pi/180
        (x,y,z)=sm.sphere2cart(phi,theta,self.coords['radius'])
        self.projection={'x':x,'y':y,'z':z}
        self.bbox={'x':[x.min(),x.max()],'y':[y.min(),y.max()],'z':[z.min(),z.max()]}
        return

    def mplot(self,fig=None,ax=None,c=None,ptype='scatter'):
        '''
        adds shape to matplotlib 3d plot

        Parameters:
        ----------
        fig (optional) matplotlib figure handle. will create a new figure if it
            does not exist.

        ax (optional)  the 3d axis handle to plot on, will create a new axis if
            creating a new figure, will use the first axis if fig is supplied
            without ax

        '''

        if fig is None:
            fig=plt.figure()

        if ax is None and len(fig.get_axes())==0:
            ax = fig.add_subplot(111, projection='3d')
        elif ax is None:
            axs=fig.get_axes()
            ax=axs[0]

        if ptype=='scatter':
            ax.scatter(self.projection['x'],self.projection['y'],self.projection['z'],c=c)
        elif ptype=='plot':
            ax.plot(self.projection['x'],self.projection['y'],self.projection['z'],c=c)

        return fig

    def buildYtSource(self,src_type='LineSource',RGBa=[1.,1.,1.,0.05],pt_size=1):
        '''
        builds the LineSource or PointSource to add to a yt view

        Parameters:
        ----------
        src_type    'LineSource' (default) or 'PointSource', corresponding to the
                    yt source types
        RGBa        4 element list-like for RGBa values to use
        '''

        clrs=[]
        x=self.projection['x'].ravel()
        y=self.projection['y'].ravel()
        z=self.projection['z'].ravel()

        if src_type=='LineSource':
            segments=[]
            for i in range(0,len(x)-1):
                segment=[[x[i],y[i],z[i]],
                        [x[i+1],y[i+1],z[i+1]]]
                segments.append(segment)
                clrs.append(RGBa)
            clrs=np.array(clrs)
            segments=np.array(segments)
            OutSource=yt.visualization.volume_rendering.api.LineSource(segments,clrs)
        elif src_type=='PointSource':
            verts=np.stack((x,y,z),axis=1)
            for vert in verts:
                clrs.append(RGBa)
            clrs=np.array(clrs)
            pt_sizes=np.full(x.shape,pt_size)
            OutSource=yt.visualization.volume_rendering.api.PointSource(verts,clrs,radii=pt_sizes)
        return OutSource

def addShapeToScene(sc,lats,lons,rads,src_type='LineSource',RGBa=[1.,1.,1.,0.005],pt_size=3):
    '''
    appends a shapeTrace to the current scene, wrapps PointSource and LineSource

    parameters
    ----------
    sc              a yt scene instance or an empty list
    lats, lons      latitude and longitude, can be scalar int or float, lists or
                    arrays of the same size and shape
    rads            radius, single value scalar or list/array same shape as
                    lats, lons
    src_type        (optional) either 'PointSource' or 'LineSource', default
                    is 'LineSource'
    RGBa            (optional) The RGBa value to use for all points or line
                    segements, default is [1.,1.,1.,0.005]
    pt_size         (optional) The pixel size of point data, default is 3

    Output
    ------
    sc    scene with shapeTrace added
    '''

    shp=shapeTrace(lats,lons,rads)
    if type(sc)==list:
        sc.append(shp)
    else:
        sc.add_source(shp.buildYtSource(src_type=src_type,RGBa=RGBa,pt_size=pt_size))
    return sc

class sphericalChunk(object):
    '''

    class for adding annotations to a spherical chunk in cartesian coordinates

    SC=sphericalChunk(lat_range,lon_range,radius_range)

    parameters
    ----------
    lat_range       list or tuple of latitude range, [min_lat, max_lat]
    lon_range       list or tuple of longitude range, [min_lon, max_lon]
    radius_range    list or tuple of radius range, [min_radius, max_radius]

    '''
    def __init__(self,lat_range,lon_range,radius_range):
        self.lat_range=lat_range
        self.lon_range=lon_range
        self.radius_range=radius_range
        return

    def domainExtent(self,sc,RGBa=[1.,1.,1.,0.005],n_latlon=100,n_rad=25):
        '''
        adds domain boundary for spherical grid interpolated to cartesian grid

        parameters
        ----------
        sc              a yt scene instance or an empty list
        RGBa            (optional) The RGBa value to use for all points or line
                        segements, default is [1.,1.,1.,0.005]
        n_latlon        (optional) number of points to use to create line segments
                        for lat/lon segment at fixed radius, default is 100
        n_rad           (optional) number of points for variable radius segments,
                        default is 25

        Output
        ------
        sc              modified yt scene
        '''

        # extents of this chunk
        lat_range=self.lat_range
        lon_range=self.lon_range
        radius_range=self.radius_range

        # constant radius, variable lat/lon boundaries
        lats=np.linspace(lat_range[0],lat_range[1],n_latlon)
        lons=np.linspace(lon_range[0],lon_range[1],n_latlon)
        for this_rad in radius_range:
            sc=addShapeToScene(sc,lats,np.full(lats.shape,lon_range[0]),this_rad,RGBa=RGBa)
            sc=addShapeToScene(sc,lats,np.full(lats.shape,lon_range[1]),this_rad,RGBa=RGBa)
            sc=addShapeToScene(sc,np.full(lons.shape,lat_range[0]),lons,this_rad,RGBa=RGBa)
            sc=addShapeToScene(sc,np.full(lons.shape,lat_range[1]),lons,this_rad,RGBa=RGBa)

        # boundary lines of constant lat/lon and variable radius
        rads=np.linspace(radius_range[0],radius_range[1],n_rad)
        rshp=rads.shape
        for lat in lat_range:
            for lon in lon_range:
                lats=np.full(rshp,lat)
                lons=np.full(rshp,lon)
                sc=addShapeToScene(sc,lats,lons,rads,RGBa=RGBa)

        return sc

    def latlonGrid(self,sc,n_lat=10,n_lon=10,radius=None,n_lat2=50,n_lon2=50,RGBa=[1.,1.,1.,0.005]):
        '''
        latlonGrid(n_lat=10,n_lon=10)

        adds a lat/lon grid at fixed radius (max radius by default)

        Parameters
        ----------
        sc      the yt scene to add to
        n_lat   (optional) number of latitudinal lines to add, default 10
        n_lon   (optional) number of longitudinal lines to add, default 10
        radius  (optional) the radius to add the grid at, default is None,
                which will pull the max radius from radius_range
        n_lat2  (optional) for a given longitude, number of lat points to use
                for line segments, default is 50
        n_lon2  (optional) for a given latitude, number of lon points to use
                for line segments, default is 50
        RGBa    (optional) The RGBa for lat/lon grid, default [1.,1.,1.,0.005]

        Output
        ------
        sc      the modified yt scene
        '''

        if radius is None:
            radius = max(self.radius_range)

        # fixed lat, vary longitude
        lat_pts=np.linspace(self.lat_range[0],self.lat_range[1],n_lat)
        lon_pts=np.linspace(self.lon_range[0],self.lon_range[1],n_lon2)
        for lat in lat_pts:
            lats=np.full(lon_pts.shape,lat)
            sc=addShapeToScene(sc,lats,lon_pts,radius,RGBa=RGBa)

        # fixed longitude, vary latitude
        lon_pts=np.linspace(self.lon_range[0],self.lon_range[1],n_lon)
        lat_pts=np.linspace(self.lat_range[0],self.lat_range[1],n_lat2)
        for lon in lon_pts:
            lons=np.full(lat_pts.shape,lon)
            sc=addShapeToScene(sc,lat_pts,lons,radius,RGBa=RGBa)

        return sc

    def wholeSphereReference(self,sc,RGBa=[1.,1.,1.,0.001],radius=None):
        '''
        adds whole sphere reference
        '''

        if radius is None:
            radius = max(self.radius_range)

        # radius lines from 0 to radius
        rads=np.linspace(0,radius,10)
        rshp=rads.shape
        for lat in self.lat_range:
            for lon in self.lon_range:
                lats=np.full(rshp,lat)
                lons=np.full(rshp,lon)
                sc=addShapeToScene(sc,lats,lons,rads,RGBa=RGBa)



        return sc

def availableShapeFiles():
    ''' returns a dictionary of available shapefiles '''

    # those included in package [filename,short name, category, description,source]

    keynames=['file','short_name','category','description','source']
    included=[
        ['cb_2018_us_state_20m.shp','us_states','political_boundaries',
         'US state boundaries, 20m resolution','https://census.gov'],
        ['GLB_VOLC.shp','global_volcanos','tectonic',
         'global volcanic fields with eruptions in last 10k years',
         'https://earthworks.stanford.edu/catalog/harvard-glb-volc'],
        ['ridge.shp','ridge','tectonic',
         'divergent plate boundaries',
         'http://www-udc.ig.utexas.edu/external/plates/data.htm'],
        ['transform.shp','transform','tectonic',
         'transform plate boundaries',
         'http://www-udc.ig.utexas.edu/external/plates/data.htm'],
        ['trench.shp','trench','tectonic',
         'convergent plate boundaries',
         'http://www-udc.ig.utexas.edu/external/plates/data.htm']
    ]

    # add on the natural earth files
    src='https://www.naturalearthdata.com/'
    res_s={'10m':'h','50m':'m','110m':'l'}
    for res in ['10m','50m','110m']:
        a=res_s[res]
        included.append(
            ['ne_'+res+'_admin_0_countries.shp','countries_'+a,'political_boundaries',src]
        )
        included.append(
            ['ne_'+res+'_admin_1_states_provinces.shp','states_'+a,'political_boundaries',src]
        )
        included.append(
            ['ne_'+res+'_coastline.shp','coast_'+a,'tectonic',src]
        )

    db=dm.filesysDB()

    shapeDict={'available':[],'details':{}}
    fullfiles=[]
    for shp in included:
        if db.validateFile(shp[0]):
            new_row=dict(zip(keynames,shp))
            shapeDict['available'].append(shp[1])
            shapeDict['details'][shp[1]]=new_row
            fullfiles.append(shp[0].split('.')[0])

    # look for other shapefiles in db
    for fi in db.FilesByDir['shapedata']:
        shrtnm=fi.split('.')[0]
        file_ext=fi.split('.')[-1]
        if shrtnm not in fullfiles and file_ext in ['shp']:
            fullfiles.append(fi)
            this_fi=[fi,shrtnm,'','','']
            shapeDict['available'].append(shrtnm)
            shapeDict['details'][shrtnm]=dict(zip(keynames,this_fi))

    return shapeDict

class shapedata(object):
    '''
    parses shapefiles using geopandas to construct yt line and point sources
    from points, lines and polygons in shapefiles.

    shp=shapedata(filename,buildTraces=True,bbox=None)

    Parameters
    ----------
    filename    the full path filename of the shapefile OR the short_name from
                availableShapeFiles()
    buildTraces if True, will build the traces on instantiating class
    bbox        bounding box to use when reading in shapefile, four element list
                [lon_min,lat_min,lon_max,lat_max]

    '''
    def __init__(self,filename,buildTraces=True,bbox=None,radius=6371.):

        self.db=dm.filesysDB()
        if os.path.isfile(filename):
            self.filename=filename
        else:
            shpDict=availableShapeFiles()
            if filename in shpDict['available']:
                # given filename is a short name
                filename=shpDict['details'][filename]['file'] # now file name
                filename=self.db.validateFile(filename) # now full path
                if filename is False:
                    raise ValueError(filename + ' does not exist.')
                else:
                    self.filename=filename
            else:
                raise ValueError(filename + ' does not exist.')

        self.radius=radius
        self.Traces=[]
        if buildTraces:
            self.Traces=self.buildTraces(bbox=bbox)

        return

    def buildTraces(self,traces=[],bbox=None,sc=None,include_points=True,
                    include_lines=True,include_polygons=True,
                    RGBa=[1.,1.,1.,0.05],pt_size=3):
        '''
        loads a shapefile and builds the yt traces.

        shapedata.buildTraces(traces=[],bbox=None,sc=None,include_points=True,
                              include_lines=True,include_polygons=True,
                              RGBa=[1.,1.,1.,0.05],pt_size=3):

        Parameters
        ----------
        traces            list of yt point or line sources
        bbox              bounding box for shapefile read
        sc                the yt scence to add traces to
        include_points    boolean, include point data from shapefile?
        include_lines     boolean, include line data from shapefile?
        include_polygons  boolean, include polygon data from shapefile?
        RGBa              RGBa list or tuple for yt
        pt_size           pixel size for point data

        Output
        ------
        if sc is provided, will return:

        sc      the modified yt scene, if

        otherwise, returns

        traces  list of yt line and point sources


        kwargs passed to shapeTrace.buildYtSource
        '''

        R0=self.radius
        def traversePoints(df,traces=[]):
            '''
            traverses shapefile points, appends to traces
            '''
            # print("traversing point data")
            # handle Points
            pt_df=df[df.geometry.type=='Point']['geometry']
            if len(pt_df)>0:
                # print(len(pt_df))
                pts=shapeTrace(pt_df.y.to_numpy(),pt_df.x.to_numpy(),R0)
                # print('appending with point size '+str(pt_size)+' and RGBa')
                # print(RGBa)
                traces.append(pts.buildYtSource('PointSource',RGBa,pt_size))

            # handle MultiPoints
            # print("traversing multipoint")
            pt_df=df[df.geometry.type=='MultiPoint'].geometry.tolist()
            lons=[]
            lats=[]
            for multipt in pt_df:
                for pt in multipt:
                    lons.append(pt.x)
                    lats.append(pt.y)
            # print('assembled lat lons for multipoint')
            # print(len(lons))
            # print(len(lats))
            if len(lons)>0 and len(lats)>0:
                pts=shapeTrace(lats,lons,R0)
                traces.append(pts.buildYtSource('PointSource',RGBa,pt_size))
            return traces

        def traverseLines(df,traces=[]):
            '''
            traverses shapefile lines, appends line segments to traces
            '''
            # Lines
            pt_df=df[df.geometry.type=='LineString'].geometry.tolist()
            for ln in pt_df:
                pts=shapeTrace(np.array(ln.xy[1]),np.array(ln.xy[0]),R0)
                traces.append(pts.buildYtSource('LineSource',RGBa))

            # MultiLines
            pt_df=df[df.geometry.type=='MultiLine'].geometry.tolist()
            for lns in pt_df:
                for ln in lns:
                    pts=shapeTrace(np.array(ln.xy[1]),np.array(ln.xy[0]),R0)
                    traces.append(pts.buildYtSource('LineSource',RGBa))
            return traces

        def traversePoly(poly,traces=[]):
            '''
            appends traces for each line segment in a polygon
            '''
            # print("    traversing single polygon")
            if poly.boundary.type=='LineString':
                # print('     points in this polygon boundary:')
                # print(len(poly.boundary.xy[1]))
                pts=shapeTrace(np.array(poly.boundary.xy[1]),np.array(poly.boundary.xy[0]),R0)
                traces.append(pts.buildYtSource('LineSource',RGBa))
            else:
                # print("    looping over multiline")
                for ln in poly.boundary:
                    # print(len(ln.xy[1]))
                    pts=shapeTrace(np.array(ln.xy[1]),np.array(ln.xy[0]),R0)
                    traces.append(pts.buildYtSource('LineSource',RGBa))
            return traces

        def traversePolygons(df,traces=[]):
            '''
            traverses shapefile polygons, appends line segments to traces
            '''
            # Polygons
            pt_df=df[df.geometry.type=='Polygon'].geometry.tolist()
            # print("traversing "+str(len(pt_df))+' polygons')
            for poly in pt_df:
                traces=traversePoly(poly,traces)

            # Multi-Polygons
            pt_df=df[df.geometry.type=='MultiPolygon'].geometry.tolist()
            # print("traversing "+str(len(pt_df))+' multi polygons')
            for multipoly in pt_df:
                for poly in multipoly:
                    traces=traversePoly(poly,traces)
            return traces

        df=gpd.read_file(self.filename,bbox=bbox)
        if include_points:
            traces=traversePoints(df,traces)
        if include_lines:
            traces=traverseLines(df,traces)
        if include_polygons:
            traces=traversePolygons(df,traces)

        if sc is not None:
            return self.addToScene(sc,traces)
        else:
            return traces


    def addToScene(self,sc,traces=[]):
        '''

        shapedata.addToScene(sc,traces)

        adds the yt line and point sources to the yt scene

        Parameters
        ----------
        sc      the yt scene to modify
        traces  list of yt line and point sources

        Output
        ------
        sc      the modified yt scene
        '''

        if len(traces)==0:
            traces=self.Traces

        for Trc in traces:
            sc.add_source(Trc)

        return sc
