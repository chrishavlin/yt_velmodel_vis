'''
loadsNWUS11-S_percent.nc and plots some slices of dVs at different depths,
saves to provided directory.
'''
import yt
from yt_velmodel_vis import seis_model as sm
import numpy as np
import os

if __name__=='__main__':
    fname='NWUS11-S_percent.nc'
    out_dir='./output'
    model=sm.netcdf(fname)
    model.coordTransform('sphere2cart')
    data={}
    data['dvs']=model.data.variables['dvs'][:]
    data['dvs'][data['dvs']>900.00000]=np.nan
    data['dvs'][data['dvs']>=0]=np.nan
    data['dvs']=data['dvs']*-1

    # yt spherical expects R, theta, phi (depth, ~lat,~lon)
    sc_mult=1.0 # scale multiplier
    bbox = model.cart['bbox']
    ds = yt.load_uniform_grid(data,data['dvs'].shape,sc_mult,bbox=bbox,nprocs=1,
                            periodicity=(False,False,False),unit_system="mks")
    sc = yt.create_scene(ds,'dvs')
    sc.camera.set_width(ds.quan(40*1000., 'm'))
    # sc.camera.rotate(np.pi/4)
    # sc.camera.roll(np.pi/6)
    # cam = sc.add_camera()
    # cam.resolution=1000
    tf = yt.ColorTransferFunction((-2, 2))
    tf.add_layers(10, w=0.01)
    source = sc.sources['source_00']
    source.set_transfer_function(tf)
    sc.save('volumetest.png',sigma_clip=6)
