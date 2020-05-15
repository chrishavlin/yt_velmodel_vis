"""
test script for new transfer function methods and plots

"""
import yt
from yt_velmodel_vis import seis_model as sm
from yt_velmodel_vis import shapeplotter as sp
from yt_velmodel_vis import transferfunctions as TFs
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.special import erf,erfc
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import cm

# choose model and plot settings
fname='NWUS11-S_percent.nc' # model selection
settings={
    'sigma_clip':1.5,
    'res_factor':1.5,
}
# settings['interp']={'field':'dvs','max_dist':50000,'res':[10000,10000,10000],
#           'input_units':'m','interpChunk':int(1e7)}
settings['interp']={'field':'dvs','max_dist':75000,'res':[10000,10000,10000],
          'input_units':'m','interpChunk':int(1e7)}
out_dir='./output/NWUS11_slab_vs_YS' # output directory for figures

# load the model
model=sm.netcdf(fname,settings['interp'])
shortname=fname.split('.')[0]
datafld=settings['interp']['field']
data={datafld:model.interp['data'][datafld]} # dict container for yt scene

X1d=model.interp['grid']['x']
Y1d=model.interp['grid']['y']
Z1d=model.interp['grid']['z']
X,Y,Z=np.meshgrid(X1d,Y1d,Z1d,indexing='ij')
nanMask=np.isnan(data['dvs'])
data['dvs'][nanMask]=99999.
sc_mult=1.0 # scale multiplier
bbox = model.cart['bbox']
ds = yt.load_uniform_grid(data,data['dvs'].shape,sc_mult,bbox=bbox,nprocs=1,
    periodicity=(True,True,True),unit_system="mks")

data['dvs'][nanMask]=-99999.
sc_mult=1.0 # scale multiplier
bbox = model.cart['bbox']
dsPos = yt.load_uniform_grid(data,data['dvs'].shape,sc_mult,bbox=bbox,nprocs=1,
    periodicity=(True,True,True),unit_system="mks")

center=np.array(ds.domain_center)
dm={'left':np.array(ds.domain_left_edge),
    'right':np.array(ds.domain_right_edge)}

dm['left']=dm['left']+200.*1000
dm['right']=dm['right']-200*1000
boxregion=ds.region(ds.domain_center,dm['left'],dm['right'])

# Draw the domain boundary and useful grids
lat_rnge=[np.min(model.data.variables['latitude']),np.max(model.data.variables['latitude'])]
lon_rnge=[np.min(model.data.variables['longitude']),np.max(model.data.variables['longitude'])]
min_dep=0.
max_dep=1200.
R=6371.
r_rnge=[(R-max_dep)*1000.,(R-min_dep)*1000.]
Chunk=sp.sphericalChunk(lat_rnge,lon_rnge,r_rnge)
Traces=[]
Traces=Chunk.domainExtent(Traces,RGBa=[1.,1.,1.,0.002],n_latlon=100,n_rad=50)
Traces=Chunk.latlonGrid(Traces,RGBa=[1.,1.,1.,0.005])
Traces=Chunk.latlonGrid(Traces,RGBa=[1.,1.,1.,0.002],radius=(R-410.)*1000.)
Traces=Chunk.latlonGrid(Traces,RGBa=[1.,1.,1.,0.002],radius=(R-max_dep)*1000.)
Traces=Chunk.wholeSphereReference(Traces,RGBa=[1.,1.,1.,0.002])


first_val=True
minval=-4
maxval=4
def scale_val(val):
    return (val-minval)/(maxval-minval)

frame=0
isovals=np.linspace(-4,4,20)
for isoval in isovals[isovals>0]:
    print("processing isoval "+str(isoval))

    if isoval>0:
        surface = dsPos.surface(boxregion, "dvs", isoval)
    else:
        surface = ds.surface(boxregion, "dvs", isoval)
    p3dc = Poly3DCollection(surface.triangles, linewidth=0.0)
    p3dc.set_facecolor(cm.seismic_r(scale_val(isoval)))
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.add_collection(p3dc)

    # add the reference grids
    for Tr in Traces:
        ax.plot(Tr.projection['x'],Tr.projection['y'],Tr.projection['z'],'k',alpha=0.3)

    ax.set_axis_off()
    if first_val:
        max_extent = (surface.vertices.max(axis=1) - surface.vertices.min(axis=1)).max()
        centers = (surface.vertices.max(axis=1) + surface.vertices.min(axis=1)) / 2
        bounds = np.zeros([3,2])
        bounds[:,0] = centers[:] - max_extent/2
        bounds[:,1] = centers[:] + max_extent/2
        print(bounds)
    ax.auto_scale_xyz(bounds[0,:], bounds[1,:], bounds[2,:])
    ax.view_init(elev=4.675324675324646,azim=-68.70967741935488)
    first_val=False
    frameid=str(frame).zfill(3)
    plt.savefig('output/iso_dvs_'+frameid+'_'+str(isoval)+'.png')
    plt.close()
    frame=frame+1
