import yt
from yt_velmodel_vis import seis_model as sm
from yt_velmodel_vis import datamanager as dm
from yt_velmodel_vis import shapeplotter as sp
import numpy as np
import os
import matplotlib.pyplot as plt
import geopandas as gpd

fname='US.2016.nc'
out_dir='./output'

model=sm.netcdf(fname)

# interpolate the cartesian (or load it if it exists)
mx=50000
dvs_n=os.path.join(out_dir,fname.split('.')[0]+'_vsv_interp_'+str(mx)+'.npy')
if os.path.isfile(dvs_n):
    model.coordTransform('sphere2cart')
    print("Loading cartesian-interpolated data from file")
    dvs = np.load(dvs_n)
else:
    model.interp2cartesian(fields=['vsv'],res=[10000,10000,10000], input_units='m',max_dist=mx)
    dvs=model.interp['data']['vsv']
    np.save(dvs_n,dvs)
data={}
data['vsv']=dvs

# set some gaussians for the TF
bnds=[2,6.]
tf = yt.ColorTransferFunction((bnds[0],bnds[1]))
# [center location, peak width, (red, green, blue, alpha)]
# TF_gausians=[[-2,1,(1.,0.0,0.0,.8)]]
# for gau in TF_gausians:
#     tf.add_gaussian(gau[0],gau[1],gau[2])
tf.add_step(4.5, 5, [0.8, 0., 0., .8])
tf.add_step(4, 4.49, [0., 0.8, 0., .8])
tf.add_step(2, 3.99, [0., 0., 0.8, .8])
# tf.add_step(-2, -1, [1.0, 0.5, 0., .8])
# tf.add_step(-1.5, -0, [1.0, 1., 0., .1])

# load the data as a uniform grid, create the 3d scene
sc_mult=1.0 # scale multiplier
bbox = model.cart['bbox'] # list-like [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
ds = yt.load_uniform_grid(data,data['vsv'].shape,sc_mult,bbox=bbox,nprocs=1,
                        periodicity=(True,True,True),unit_system="mks")

sc = yt.create_scene(ds,'vsv')

# plot the TF with a histogram
x = np.linspace(bnds[0],bnds[1],tf.nbins)
y = tf.funcs[3].y
w = np.append(x[1:]-x[:-1], x[-1]-x[-2])
colors = np.array([tf.funcs[0].y, tf.funcs[1].y, tf.funcs[2].y,
                   tf.funcs[3].y]).T
fig = plt.figure(figsize=[6, 3])
ax = fig.add_axes([0.2, 0.2, 0.75, 0.75])
d_hist=ax.hist(data['vsv'][~np.isnan(dvs)].ravel(),bins=100,density=True,log=False)
ax.bar(x, tf.funcs[3].y, w, edgecolor=[0.0, 0.0, 0.0, 0.0],
       log=False, color=colors, bottom=[0])
plt.savefig(os.path.join(out_dir,'CP_shenritz_tf.png'))


# draw true boundary extents
lat_rnge=[np.min(model.data.variables['latitude']),np.max(model.data.variables['latitude'])]
lon_rnge=[np.min(model.data.variables['longitude']),np.max(model.data.variables['longitude'])]
print(lat_rnge)
print(lon_rnge)
min_dep=0.
max_dep=1200.
R=6371.
r_rnge=[(R-max_dep)*1000.,(R-min_dep)*1000.]
Chunk=sp.sphericalChunk(lat_rnge,lon_rnge,r_rnge)
sc=Chunk.domainExtent(sc,RGBa=[1.,1.,1.,0.005],n_latlon=100,n_rad=50)
sc=Chunk.latlonGrid(sc,RGBa=[1.,1.,1.,0.005])
sc=Chunk.latlonGrid(sc,RGBa=[1.,1.,1.,0.001],radius=(R-410.)*1000.)
sc=Chunk.latlonGrid(sc,RGBa=[1.,1.,1.,0.001],radius=(R-max_dep)*1000.)
sc=Chunk.wholeSphereReference(sc,RGBa=[1.,1.,1.,0.001])

YS_lat=np.array([44.429764])
YS_lon=np.array([-110.584663])
YS_rads=np.array([6371.*1e3])
sc=sp.addShapeToScene(sc,YS_lat,YS_lon,YS_rads,'PointSource',[1.,.9,.9,0.005],6)


## add data from a shapefile
print('adding volcanic fields')
shp_bbox=[lon_rnge[0]-360,lat_rnge[0],lon_rnge[1]-360,lat_rnge[1]]
print(shp_bbox)
volcs=sp.shapedata('global_volcanos',radius=R*1000.,buildTraces=False)
sc=volcs.buildTraces(RGBa=[0.,0.8,0.,0.05],bbox=shp_bbox,sc=sc)

print("adding state bounds")
continents=sp.shapedata('us_states',bbox=shp_bbox,radius=R*1000.)
sc=continents.addToScene(sc)


print("adding tectonic boundaries")
clrs={
    'transform':[0.8,0.,0.8,0.05],
    'ridge':[0.,0.,0.8,0.05],
    'trench':[0.8,0.,0.,0.05],
}
for bound in ['transform','ridge','trench']:
    tect=sp.shapedata(bound,radius=R*1000.,buildTraces=False)
    sc=tect.buildTraces(RGBa=clrs[bound],sc=sc,bbox=shp_bbox)


Rmax=6371*1000.
center_vec=np.array([np.mean(bbox[0])/Rmax,np.mean(bbox[1])/Rmax,np.mean(bbox[2])/Rmax])

pos=sc.camera.position
print(pos)
sc.camera.set_position(pos,north_vector=center_vec)
source = sc.sources['source_00']
source.tfh.set_log(False)


res=sc.camera.get_resolution()
res_factor=1
new_res=(res[0]*res_factor,res[1]*res_factor)
sc.camera.set_resolution(new_res)
# apply the TF and render it
source.set_transfer_function(tf)

nm='CP_shrenritz_vol.png'
print("saving "+nm)
sc.save(os.path.join(out_dir,nm),sigma_clip=0.5)
