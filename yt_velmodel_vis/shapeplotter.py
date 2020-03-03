'''
shapeplotter module: support for plotting shapefile data on yt scenes
'''

from . import seis_model as sm
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import yt

class shapedata(object):
    '''
    shapefile support
    '''
    def __init__(self,filename):
        ncdf=sm.netcdf(filename,load_ds=False) # use for file path validation
        self.file=ncdf.validateFile()

        return

class shapeTrace(object):
    '''
    individual trace of a shape.

    shapeTrace(lat,lon,radius) with lat,lon in degrees.

    lat, lon can be scalar int or float, lists or arrays of the same size and
    shape, radius can be single value or the same
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
            print(clrs.shape)
            print(segments.shape)
            OutSource=yt.visualization.volume_rendering.api.LineSource(segments,clrs)
        elif src_type=='PointSource':
            verts=np.stack((x,y,z),axis=1)
            for vert in verts:
                clrs.append(RGBa)
            clrs=np.array(clrs)
            pt_sizes=np.full(x.shape,pt_size)
            OutSource=yt.visualization.volume_rendering.api.PointSource(verts,clrs,radii=pt_sizes)
        return OutSource
