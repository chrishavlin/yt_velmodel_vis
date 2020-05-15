"""
test script for new transfer function methods and plots

"""
import yt
from yt_velmodel_vis import seis_model as sm, shapeplotter as sp
from yt_velmodel_vis import transferfunctions as TFs, animator as ani
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.special import erf,erfc


# choose model and plot settings
fname='NWUS11-S_percent.nc' # model selection
settings={
    'sigma_clip':1.5,
    'res_factor':1.5,
}

settings['interp']={'field':'dvs','max_dist':75000,'res':[10000,10000,10000],
          'input_units':'m','interpChunk':int(1e7)}
out_dir=os.path.join(os.environ.get('YTVELMODELDIR'),'output','NWUS11-S_animate') # output directory for figures
if os.path.isdir(out_dir) is False:
    os.mkdir(out_dir)
# load the model
model=sm.netcdf(fname,settings['interp'])
shortname=fname.split('.')[0]
datafld=settings['interp']['field']
data={datafld:model.interp['data'][datafld]} # dict container for yt scene

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
sc=Chunk.latlonGrid(sc,RGBa=[1.,1.,1.,0.009])
sc=Chunk.latlonGrid(sc,RGBa=[1.,1.,1.,0.006],radius=(R-410.)*1000.)
sc=Chunk.latlonGrid(sc,RGBa=[1.,1.,1.,0.002],radius=(R-max_dep)*1000.)
sc=Chunk.wholeSphereReference(sc,RGBa=[1.,1.,1.,0.002])

# Add shapefile data
print('adding shapefile data to scene')
shp_bbox=[lon_rnge[0],lat_rnge[0],lon_rnge[1],lat_rnge[1]]

for shpfi in ['us_states']:
    thisshp=sp.shapedata(shpfi,bbox=shp_bbox,radius=R*1000.)
    sc=thisshp.addToScene(sc)

clrs={
    'transform':[0.8,0.,0.8,0.05],
    'ridge':[0.,0.,0.8,0.05],
    'trench':[0.8,0.,0.,0.05],
    'global_volcanos':[0.,0.8,0.,0.05]
}
for bound in ['transform','ridge','trench','global_volcanos']:
    tect=sp.shapedata(bound,radius=R*1000.,buildTraces=False)
    sc=tect.buildTraces(RGBa=clrs[bound],sc=sc,bbox=shp_bbox)



# build and set the transfer function
tfOb = TFs.dv(data[datafld].ravel(),bounds=[-4,4])

# segment 1, slow anomalies
bnds=[-1.3,-.3]
TFseg=TFs.TFsegment(tfOb,bounds=bnds,cmap='OrRd_r')
alpha_o=0.95
Dalpha=-0.85
alpha=alpha_o + Dalpha/(bnds[1]-bnds[0]) * (TFseg.dvbins_c-bnds[0])
tfOb.addTFsegment(alpha,TFseg)

# segment 2, fast anomalies
bnds=[.1,.35]
TFseg=TFs.TFsegment(tfOb,bounds=bnds,cmap='winter_r')
alpha_o=.6
Dalpha=.4
alpha=alpha_o+ Dalpha/(bnds[1]-bnds[0]) * (TFseg.dvbins_c-bnds[0])
tfOb.addTFsegment(alpha,TFseg)

# plot the TF and save it
f=plt.figure()
ax=plt.axes()
ax=tfOb.addHist(ax=ax,density=True,color=(0.,0.,0.,1.))
ax=tfOb.addTFtoPlot(ax=ax)
ax.set_xlabel('$\mathregular{dV_s}$')
plt.savefig(os.path.join(out_dir,'transfer_function.png'))

source = sc.sources['source_00']
source.tfh.set_log(False)
# source.set_transfer_function(tf) # apply the TF and render it
source.set_transfer_function(tfOb.tf)

# adjust initial view

# initial zoom
zoom_factor=.7 # < 1 zooms in
init_width=sc.camera.width
sc.camera.width = (init_width * zoom_factor)


# use vector passing through center as north
x_c=np.mean(bbox[0])
y_c=np.mean(bbox[1])
z_c=np.mean(bbox[2])
center_vec=np.array([x_c,y_c,z_c])
center_vec=center_vec / np.linalg.norm(center_vec)

# use initial position, adjust the north_vector
pos=sc.camera.position
sc.camera.set_position(pos,north_vector=center_vec)

# adjust resolution
res=sc.camera.get_resolution()
res_factor=settings.get('res_factor', 1)
new_res=(int(res[0]*res_factor),int(res[1]*res_factor))
sc.camera.set_resolution(new_res)


# set up frame directory
frameDir=os.path.join(out_dir,'frames')
if os.path.isdir(frameDir) is False:
    os.mkdir(frameDir)

frame_sigdigs=4
frameN=0
init_pos=[i for i in pos]

# adjust viewing angle
for fracval in np.linspace(1,1.3,15):
    pos[0]=init_pos[0]*fracval
    sc.camera.set_position(pos,north_vector=center_vec)

    frameid=str(frameN).zfill(frame_sigdigs)
    nm='vol_render_'+frameid+'.png'
    sc.save(os.path.join(frameDir,nm),sigma_clip=settings.get('sigma_clip', 2))

    frameN=frameN+1

# zoom in
for zoomval in np.linspace(zoom_factor,0.4,5):
    sc.camera.width = (init_width * zoomval)
    frameid=str(frameN).zfill(frame_sigdigs)
    nm='vol_render_'+frameid+'.png'
    sc.save(os.path.join(frameDir,nm),sigma_clip=settings.get('sigma_clip', 2))
    frameN=frameN+1

# rotate 360 degrees
for deg in np.linspace(0,360,361):
    sc.camera.rotate(np.pi/180., rot_vector=center_vec)
    frameid=str(frameN).zfill(frame_sigdigs)
    nm='vol_render_'+frameid+'.png'
    sc.save(os.path.join(frameDir,nm),sigma_clip=settings.get('sigma_clip', 2))
    frameN=frameN+1

# build a gif (will be lower quality than individual frames)
an=ani.animator(frameDir)
an.writeGif(os.path.join(out_dir,'NWUS11-S-rotation.gif'),sec_per_frame=.1)

# to build a higher quality video, use a tool like ffmpeg.
# For ffmpeg, cd into directory with the frames then, e.g., :
# ffmpeg -framerate 10 -pattern_type glob -i '*.png' -c:v libx264 -r 30 -pix_fmt yuv420p NWUS11-S-rotation.mp4
