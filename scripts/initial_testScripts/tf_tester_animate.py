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


# choose model and plot settings
fname='NWUS11-S_percent.nc' # model selection
settings={
    'interp':{'field':'dvs','max_dist':50000,'res':[10000,10000,10000],
              'input_units':'m','interpChunk':int(1e7)},
    'sigma_clip':.4,
    'res_factor':3,
}
out_dir='./output/tf_tester/animate' # output directory for figures

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

def set_the_scene():
    """ creates the yt scene, except for transfer function and rendering """
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
    sc=Chunk.latlonGrid(sc,RGBa=[1.,1.,1.,0.002],radius=(R-410.)*1000.)
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


    # some camera settings
    pos=sc.camera.position
    Rmax=6371*1000.
    center_vec=np.array([np.mean(bbox[0])/Rmax,np.mean(bbox[1])/Rmax,np.mean(bbox[2])/Rmax])
    sc.camera.set_position(pos,north_vector=center_vec)
    source = sc.sources['source_00']
    source.tfh.set_log(False)

    res=sc.camera.get_resolution()
    res_factor=settings.get('res_factor', 1)
    new_res=(int(res[0]*res_factor),int(res[1]*res_factor))
    sc.camera.set_resolution(new_res)

    zoom_factor=0.7 # < 1 zooms in
    init_width=sc.camera.width
    sc.camera.width = (init_width * zoom_factor)
    sc.camera.rotate(0*np.pi/180)

    return (sc,source)

sc,source=set_the_scene()

def plotSaveTf(tfOb,savename):
    f=plt.figure()
    ax=plt.axes()
    ax=tfOb.addHist(ax=ax,density=True,color=(0.,0.,0.,1.))
    ax=tfOb.addTFtoPlot(ax=ax)
    ax.set_xlabel('dvs')
    ax.set_title(savename)
    if os.path.exists(out_dir) is False:
        os.mkdir(out_dir)
    plt.savefig(os.path.join(out_dir,savename+'.png'))
    return

def TF_test_6(centerval,dvwid=1,cmap='summer_r'):
    tfOb = TFs.dv(data[datafld].ravel(),bounds=[-4,4])
    bnds=[centerval-dvwid/2,centerval+dvwid/2]
    TFseg=TFs.TFsegment(tfOb,bounds=bnds,cmap=cmap)
    dv=TFseg.dvbins_c # dv in this segment
    alpha=np.full(TFseg.dvbins_c.shape,0.6)
    tfOb.addTFsegment(alpha,TFseg)

    return tfOb

frame=0
for centerval in np.linspace(-2,-.51,20):
    framestr=str(frame).zfill(3)
    tfOb=TF_test_6(centerval,0.5,'OrRd')
    plotSaveTf(tfOb,'hot_TF'+framestr)
    source.set_transfer_function(tfOb.tf) # apply the TF and render it
    nm='hot_vol_'+framestr+'.png'
    sc.save(os.path.join(out_dir,nm),sigma_clip=settings.get('sigma_clip', 2))
    frame=frame+1


frame=0
for centerval in np.linspace(1.,.251,20):
    framestr=str(frame).zfill(3)
    tfOb=TF_test_6(centerval,0.25,'cool')
    plotSaveTf(tfOb,'cold_TF'+framestr)
    source.set_transfer_function(tfOb.tf) # apply the TF and render it
    nm='cold_vol_'+framestr+'.png'
    sc.save(os.path.join(out_dir,nm),sigma_clip=settings.get('sigma_clip', 2))
    frame=frame+1
