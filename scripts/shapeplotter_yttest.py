import yt
from yt_velmodel_vis import seis_model as sm
from yt_velmodel_vis import shapeplotter as sp
import numpy as np
import os
import matplotlib.pyplot as plt

fname='NWUS11-S_percent.nc'
out_dir='./output'
model=sm.netcdf(fname)

# interpolate the cartesian (or load it if it exists)
mx=50000
dvs_n=os.path.join(out_dir,fname.split('.')[0]+'_dvs_interp_'+str(mx)+'.npy')
if os.path.isfile(dvs_n):
    model.coordTransform('sphere2cart')
    print("Loading cartesian-interpolated data from file")
    dvs = np.load(dvs_n)
else:
    model.interp2cartesian(fields=['dvs'],res=[10000,10000,10000], input_units='m',max_dist=mx)
    dvs=model.interp['data']['dvs']
    np.save(dvs_n,dvs)
data={}
data['dvs']=dvs

# set some gaussians for the TF
bnds=[-6,10.]
tf = yt.ColorTransferFunction((bnds[0],bnds[1]))
# [center location, peak width, (red, green, blue, alpha)]
TF_gausians=[[-2,1,(1.,0.0,0.0,0.8)],
             [2,1,(0.,0.0,1.0,0.8)]]
for gau in TF_gausians:
    tf.add_gaussian(gau[0],gau[1],gau[2])

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
plt.savefig(os.path.join(out_dir,'shapeplotter_test_tf.png'))

# load the data as a uniform grid, create the 3d scene
sc_mult=1.0 # scale multiplier
bbox = model.cart['bbox'] # list-like [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
print(bbox)
ds = yt.load_uniform_grid(data,data['dvs'].shape,sc_mult,bbox=bbox,nprocs=1,
                        periodicity=(True,True,True),unit_system="mks")

sc = yt.create_scene(ds,'dvs')

# Draw the domain boundary
# sc.annotate_domain(ds, color=[1, 1, 1, 0.001])

# draw true boundary extents
lat_rnge=[np.min(model.data.variables['latitude']),np.max(model.data.variables['latitude'])]
lon_rnge=[np.min(model.data.variables['longitude']),np.max(model.data.variables['longitude'])]
min_dep=0.
max_dep=1200.
radval=6371*1000.
boundTrs=[]
lats=np.linspace(lat_rnge[0],lat_rnge[1],100)
lons=np.linspace(lon_rnge[0],lon_rnge[1],100)


for depth in [min_dep,max_dep]:
    this_rad=(6371-depth)*1000.
    # this_sh=sp.shapeTrace(lats,np.full(lats.shape,lon_rnge[0]),this_rad)
    boundTrs.append(sp.shapeTrace(lats,np.full(lats.shape,lon_rnge[0]),this_rad))
    boundTrs.append(sp.shapeTrace(lats,np.full(lats.shape,lon_rnge[1]),this_rad))
    boundTrs.append(sp.shapeTrace(np.full(lons.shape,lat_rnge[0]),lons,this_rad))
    boundTrs.append(sp.shapeTrace(np.full(lons.shape,lat_rnge[1]),lons,this_rad))

# surface lat grid
lons=np.linspace(lon_rnge[0],lon_rnge[1],10)
this_rad=6371*1000
for lon in lons:
    boundTrs.append(sp.shapeTrace(lats,np.full(lats.shape,lon),this_rad))


rads=(6371-np.linspace(min_dep,max_dep,50))*1000
boundTrs.append(sp.shapeTrace(np.full(rads.shape,lat_rnge[0]),np.full(rads.shape,lon_rnge[0]),rads))
boundTrs.append(sp.shapeTrace(np.full(rads.shape,lat_rnge[0]),np.full(rads.shape,lon_rnge[1]),rads))
boundTrs.append(sp.shapeTrace(np.full(rads.shape,lat_rnge[1]),np.full(rads.shape,lon_rnge[0]),rads))
boundTrs.append(sp.shapeTrace(np.full(rads.shape,lat_rnge[1]),np.full(rads.shape,lon_rnge[1]),rads))
#
for Tr in boundTrs:
    ln_src=Tr.buildYtSource(src_type='LineSource',RGBa=[1.,1.,1.,0.001])
    sc.add_source(ln_src)


YS_lat=np.array([44.429764])
YS_lon=np.array([-110.584663])
YS_rads=np.array([6371.*1e3])
YS=sp.shapeTrace(YS_lat,YS_lon,YS_rads)
sc.add_source(YS.buildYtSource(src_type='PointSource',RGBa=[1.,1.,1.,0.005],pt_size=3))
# some camera settings

pos=sc.camera.position
print(pos)
sc.camera.set_position(pos,north_vector=np.array([0.0, 0.0, 1.0]))
source = sc.sources['source_00']
source.tfh.set_log(False)


res=sc.camera.get_resolution()
print(res)
res_factor=2
new_res=(res[0]*res_factor,res[1]*res_factor)
sc.camera.set_resolution(new_res)
# apply the TF and render it
source.set_transfer_function(tf)
current_angle=0
angle_step=5;
for frame in range(0,10):
    sc.camera.rotate(current_angle)
    current_angle=current_angle+angle_step*np.pi/180
    nm='shapeplotter/shapeplotter_'+str(frame)+'.png'
    sc.save(os.path.join(out_dir,nm),sigma_clip=3)
