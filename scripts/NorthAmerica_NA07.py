"""
example script: North America continent scale using model NA07

shear wave perturbations

will output figures in ./output/NA07_percent/

figures include a volume rendering and the transfer function used in volume
rendering.
"""
import yt
from yt_velmodel_vis import seis_model as sm
from yt_velmodel_vis import shapeplotter as sp
import numpy as np
import os
import matplotlib.pyplot as plt


# choose model and plot settings
fname='NA07_percent.nc' # model selection
settings = {
    'interp':{'field':'dvs','max_dist':50000,'res':[10000,10000,10000],
              'interpChunk':int(1e7)},
    'tfs':[
          {'tftype':'step','values':(-5, -.5, [1.0, 0., 0., .5])},
          {'tftype':'step','values':(0.0001, 0.5, [0.5, 1., .5, .6])},
          {'tftype':'step','values':(0.5, 1, [0., 1., .5, .7])},
          {'tftype':'step','values':(1., 1.5, [0., 0., .9, .8])}
        ],
    'tf_bounds':[-10.,10.],
    'res_factor':1.,
    'sigma_clip':1
}
out_dir='./output/'+fname.split('.')[0] # output directory for figures

# interpDict=dict(field='dvs',res=[10000,10000,10000], input_units='m',max_dist=50000,interpChunk=int(1e6))
# load the model
model=sm.netcdf(fname,settings['interp'])
shortname=fname.split('.')[0]
datafld=settings['interp']['field']
data={datafld:model.interp['data'][datafld]} # dict container for yt scene

# add the transfer functions
bnds=settings.get('tf_bounds', [-5,5])
tf = yt.ColorTransferFunction((bnds[0],bnds[1]))
for TFtoadd in settings['tfs']:
    v=TFtoadd['values']
    if TFtoadd['tftype']=='step':
        tf.add_step(v[0],v[1],v[2])
    elif TFtoadd['tftype']=='gaussian':
        tf.add_gaussian(v[0],v[1],v[2])

# plot the TF with a histogram
x = np.linspace(bnds[0],bnds[1],tf.nbins)
y = tf.funcs[3].y
w = np.append(x[1:]-x[:-1], x[-1]-x[-2])
colors = np.array([tf.funcs[0].y, tf.funcs[1].y, tf.funcs[2].y,
                   tf.funcs[3].y]).T
fig = plt.figure(figsize=[6, 3])
ax = fig.add_axes([0.2, 0.2, 0.75, 0.75])
d_hist=ax.hist(data[datafld][~np.isnan(data[datafld])].ravel(),bins=100,density=True,log=False)
ax.bar(x, tf.funcs[3].y, w, edgecolor=[0.0, 0.0, 0.0, 0.0],
       log=False, color=colors, bottom=[0])

if os.path.exists(out_dir) is False:
    os.mkdir(out_dir)

plt.savefig(os.path.join(out_dir,'transfer_function.png'))

# load the data as a uniform grid, create the 3d scene
sc_mult=1.0 # scale multiplier
bbox = model.cart['bbox'] # list-like [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
ds = yt.load_uniform_grid(data,data[datafld].shape,sc_mult,bbox=bbox,nprocs=1,
                        periodicity=(True,True,True),unit_system="mks")

sc = yt.create_scene(ds,datafld)

# Draw the domain boundary and useful grids
lat_rnge=[np.min(model.data.variables['latitude']),np.max(model.data.variables['latitude'])]
lon_rnge=[np.min(model.data.variables['longitude']),np.max(model.data.variables['longitude'])]
min_dep=0.
max_dep=1200.
R=6371.
r_rnge=[(R-max_dep)*1000.,(R-min_dep)*1000.]
Chunk=sp.sphericalChunk(lat_rnge,lon_rnge,r_rnge)
sc=Chunk.domainExtent(sc,RGBa=[1.,1.,1.,0.002],n_latlon=100,n_rad=50)
sc=Chunk.latlonGrid(sc,RGBa=[1.,1.,1.,0.005])
sc=Chunk.latlonGrid(sc,RGBa=[1.,1.,1.,0.002],radius=(R-max_dep)*1000.)
sc=Chunk.wholeSphereReference(sc,RGBa=[1.,1.,1.,0.002])

# Add shapfeil data
print('adding volcanic fields')
shp_bbox=[lon_rnge[0],lat_rnge[0],lon_rnge[1],lat_rnge[1]]
volcs=sp.shapedata('global_volcanos',radius=R*1000.,buildTraces=False)
sc=volcs.buildTraces(RGBa=[0.,0.8,0.,0.05],bbox=shp_bbox,sc=sc)

print("adding state bounds")
continents=sp.shapedata('us_states',bbox=shp_bbox,radius=R*1000.)
sc=continents.addToScene(sc)

print('adding plate boundaries')
clrs={
    'transform':[0.8,0.,0.8,0.05],
    'ridge':[0.,0.,0.8,0.05],
    'trench':[0.8,0.,0.,0.05],
}
for bound in ['transform','ridge','trench']:
    tect=sp.shapedata(bound,radius=R*1000.,buildTraces=False)
    sc=tect.buildTraces(RGBa=clrs[bound],sc=sc,bbox=shp_bbox)

# some camera settings
pos=sc.camera.position
camera_height=.5*6371.*1000# + 1000*1000.
cam_xyz=sm.sphere2cart(90.-0.,np.mean(lon_rnge),camera_height)

sc.camera.set_position(cam_xyz,north_vector=np.array([0.0, 0.0, 1.0]))
source = sc.sources['source_00']
source.tfh.set_log(False)

res=sc.camera.get_resolution()
res_factor=settings.get('res_factor', 1)
new_res=(int(res[0]*res_factor),int(res[1]*res_factor))
sc.camera.set_resolution(new_res)

zoom_factor=.7 # < 1 zooms in
init_width=sc.camera.width
sc.camera.width = (init_width * zoom_factor)
# sc.camera.rotate(-5*np.pi/180)

source.set_transfer_function(tf) # apply the TF and render it


nm='vol_render.png'
sc.save(os.path.join(out_dir,nm),sigma_clip=settings.get('sigma_clip', 2))
