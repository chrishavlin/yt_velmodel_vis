import yt
from yt_velmodel_vis import seis_model as sm
from yt_velmodel_vis import datamanager as dm
from yt_velmodel_vis import shapeplotter as sp
import numpy as np
import os
import matplotlib.pyplot as plt
import geopandas as gpd

fname='NWUS11-S_percent.nc'
fname='NA07_percent.nc'
out_dir='./output'

# load model, include interpolation
interpDict=dict(field='dvs',res=[10000,10000,10000], input_units='m',max_dist=50000,interpChunk=int(1e6))
model=sm.netcdf(fname,interpDict)

# interpolate the cartesian (or load it if it exists)
data={}
data['dvs']=model.interp['data']['dvs']


# load the data as a uniform grid, create the 3d scene
sc_mult=1.0 # scale multiplier
bbox = model.cart['bbox'] # list-like [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
ds = yt.load_uniform_grid(data,data['dvs'].shape,sc_mult,bbox=bbox,nprocs=1,
                        periodicity=(True,True,True),unit_system="mks")
sc = yt.create_scene(ds,'dvs')


# add TF
bnds=[-6,10.]
tf = yt.ColorTransferFunction((bnds[0],bnds[1]))
tf.add_step(-3, -2, [1.0, 0., 0., 1.0])
tf.add_step(-2, -1, [1.0, 0.5, 0., .8])

# Draw boundary extents, helpful grids
lat_rnge=[np.min(model.data.variables['latitude']),np.max(model.data.variables['latitude'])]
lon_rnge=[np.min(model.data.variables['longitude']),np.max(model.data.variables['longitude'])]
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

# Add data from a shapefile
shp_bbox=[lon_rnge[0],lat_rnge[0],lon_rnge[1],lat_rnge[1]]
volcs=sp.shapedata('global_volcanos',radius=R*1000.,buildTraces=False)
sc=volcs.buildTraces(RGBa=[0.,0.8,0.,0.05],bbox=shp_bbox,sc=sc)
continents=sp.shapedata('us_states',bbox=shp_bbox,radius=R*1000.)
sc=continents.addToScene(sc)
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
sc.camera.set_position(pos,north_vector=center_vec)
source = sc.sources['source_00']
source.tfh.set_log(False)

res=sc.camera.get_resolution()
res_factor=1
new_res=(res[0]*res_factor,res[1]*res_factor)
sc.camera.set_resolution(new_res)
source.set_transfer_function(tf)

nm='interpolate_g5py_test_vol.png'
print("saving "+nm)
sc.save(os.path.join(out_dir,nm),sigma_clip=0.5)
