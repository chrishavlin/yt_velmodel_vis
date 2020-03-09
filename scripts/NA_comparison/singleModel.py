import yt
from yt_velmodel_vis import seis_model as sm
from yt_velmodel_vis import shapeplotter as sp
import numpy as np
import os
import matplotlib.pyplot as plt

model_settings={
    'NWUS11-S_percent.nc': {
        'field':'dvs',
        'interp':{'max_dist':50000,'res':[10000,10000,10000]},
        'view':{'rotate_0':0.,'north_vec':'center_vec'}
        },
    'NA07_percent.nc' :{
        'field':'dvs',
        'interp':{'max_dist':100000,'res':[50000,50000,50000]},
        'view':{
            'rotate_0':0.,'north_vec':'north_vec'}
        }
}

# '3D2018-08Sv-radius.nc'
# 3D2018-08Sv-radius.nc
# csem-north-america-2019.12.01.nc
# DNA09_percent.nc                  Porter_Liu_Holt_2015.nc
# DNA10-S_percent.nc                SEMum-NA14_kmps.nc
# MC35.csv
# NA04_percent.nc                   Shen_Ritzwoller_2016.nc
# NA07_percent.nc


fname='NA07_percent.nc'
out_dir='./output/'+fname.split('.')[0] # for figures
derived_data_dir='./data' # for derived data (interpolations )
model=sm.netcdf(fname)

shortname=fname.split('.')[0]
settings=model_settings[fname]
datafld=settings['field']

# interpolate to the cartesian grid (or load it if it exists)
mx=settings['interp']['max_dist']
res=settings['interp']['res']
interp_file='_'.join([shortname,datafld,str(mx),str(res[0])])+'.npy'
dvs_n=os.path.join(derived_data_dir,interp_file)
if os.path.isfile(dvs_n):
    model.coordTransform('sphere2cart')
    print("Loading cartesian-interpolated data from file")
    dvs = np.load(dvs_n)
else:
    model.interp2cartesian(fields=[datafld],res=res, input_units='m',max_dist=mx)
    dvs=model.interp['data'][datafld]
    np.save(dvs_n,dvs)
# dvs[dvs<0]=0
data={}
data[datafld]=dvs

# set some gaussians for the TF
bnds=[-10.,10.]
tf = yt.ColorTransferFunction((bnds[0],bnds[1]))
# [center location, peak width, (red, green, blue, alpha)]
# TF_gausians=[[-2,1,(1.,0.0,0.0,.8)]]
# for gau in TF_gausians:
#     tf.add_gaussian(gau[0],gau[1],gau[2])
tf.add_step(-5, -0.001, [1.0, 0., 0., .5])
tf.add_step(0.0001, 0.5, [0.5, 1., .5, .6])
tf.add_step(0.5, 1, [0., 1., .5, .7])
tf.add_step(1., 1.5, [0., 0., .9, .8])

# tf.add_step(-1,-0.0001, [1., 0.,0.,.3])
# tf.add_step(1,2, [0., 1., 1., .8])
# tf.add_step(2,4, [0., 0., 1., .9])
# tf.add_step(-1.5, -0, [1.0, 1., 0., .1])
# tf.add_gaussian(4,3,(1.,0,0,0.8))
# tf.add_gaussian(-2,1,(1.,0,0,0.8))
# tf.add_gaussian(2,1.5,(1.,0,0,0.8))


# plot the TF with a histogram
x = np.linspace(bnds[0],bnds[1],tf.nbins)
y = tf.funcs[3].y
w = np.append(x[1:]-x[:-1], x[-1]-x[-2])
colors = np.array([tf.funcs[0].y, tf.funcs[1].y, tf.funcs[2].y,
                   tf.funcs[3].y]).T
fig = plt.figure(figsize=[6, 3])
ax = fig.add_axes([0.2, 0.2, 0.75, 0.75])
d_hist=ax.hist(data['dvs'][~np.isnan(dvs)].ravel(),bins=100,density=True,log=False)
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

# Draw the domain boundary

# draw true boundary extents
lat_rnge=[np.min(model.data.variables['latitude']),np.max(model.data.variables['latitude'])]
lon_rnge=[np.min(model.data.variables['longitude']),np.max(model.data.variables['longitude'])]
min_dep=0.
max_dep=1200.
R=6371.
r_rnge=[(R-max_dep)*1000.,(R-min_dep)*1000.]
Chunk=sp.sphericalChunk(lat_rnge,lon_rnge,r_rnge)
sc=Chunk.domainExtent(sc,RGBa=[1.,1.,1.,0.001],n_latlon=100,n_rad=50)
sc=Chunk.latlonGrid(sc,RGBa=[1.,1.,1.,0.001])
# sc=Chunk.latlonGrid(sc,RGBa=[1.,1.,1.,0.001],radius=(R-410.)*1000.)
sc=Chunk.latlonGrid(sc,RGBa=[1.,1.,1.,0.001],radius=(R-max_dep)*1000.)
sc=Chunk.wholeSphereReference(sc,RGBa=[1.,1.,1.,0.001])

YS_lat=np.array([44.429764])
YS_lon=np.array([-110.584663])
YS_rads=np.array([6371.*1e3])
sc=sp.addShapeToScene(sc,YS_lat,YS_lon,YS_rads,'PointSource',[1.,1.,1.,0.005],6)

# set north through domain center
# bbox = model.cart['bbox'] # list-like [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
if settings['view']['north_vec']=='center_vec':
    Rmax=6371*1000.
    center_vec=np.array([np.mean(bbox[0])/Rmax,np.mean(bbox[1])/Rmax,np.mean(bbox[2])/Rmax])
elif settings['view']['north_vec']=='north_vec':
    center_vec=np.array([0.0, 0.0, 1.0])

# some camera settings
pos=sc.camera.position
camera_height=.5*6371.*1000# + 1000*1000.
cam_xyz=sm.sphere2cart(90.-0.,np.mean(lon_rnge),camera_height)
sc.camera.set_position(cam_xyz,north_vector=center_vec)
source = sc.sources['source_00']
source.tfh.set_log(False)


res=sc.camera.get_resolution()
res_factor=1
new_res=(res[0]*res_factor,res[1]*res_factor)
sc.camera.set_resolution(new_res)
# apply the TF and render it
source.set_transfer_function(tf)

initRotation=settings['view']['rotate_0']*np.pi/180
sc.camera.rotate(initRotation)

nm='vol_render.png'
sc.save(os.path.join(out_dir,nm),sigma_clip=0.5)
