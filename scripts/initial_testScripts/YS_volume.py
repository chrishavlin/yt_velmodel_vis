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


    center_lat=43.
    lat_range=center_lat+np.array([-5,5])
    center_lon=-117.
    lon_range=center_lon+np.array([-5,10])
    dep_range=np.array([100,1000])

    phi=(90-lat_range) * np.pi / 180.
    theta=(360+lon_range)*np.pi/180.
    R=(6371. - dep_range)*1000.
    x,y,z=sm.sphere2cart(phi,theta,R)


    lats=np.array(model.data.variables['latitude'])
    lons=np.array(model.data.variables['longitude'])
    deps=np.array(model.data.variables['depth'])

    m2=(lats>=lat_range[0])&(lats<=lat_range[1])
    m3=(lons>=lon_range[0])&(lons<=lon_range[1])
    m1=(deps>=dep_range[0])&(deps<=dep_range[1])
    data['dvs']=data['dvs'][m1,:,:]
    data['dvs']=data['dvs'][:,m2,:]
    data['dvs']=data['dvs'][:,:,m3]
    data['dvs']=np.flipud(data['dvs'])
    # data['dvs']=np.transpose(data['dvs'],axes=[0,1,2])
    bbox=np.array([x,y,z])

    # yt spherical expects R, theta, phi (depth, ~lat,~lon)
    sc_mult=1.0 # scale multiplier
    # bbox = model.cart['bbox']
    ds = yt.load_uniform_grid(data,data['dvs'].shape,sc_mult,bbox=bbox,nprocs=1,
                            periodicity=(False,False,False),unit_system="mks")

    sc = yt.create_scene(ds,'dvs')
    sc.camera.set_width(ds.quan(15*1000., 'm'))
    # sc.camera.rotate(10*np.pi/180.)
    # sc.camera.roll(-5*np.pi/180)
    # sc.camera.yaw(-5*np.pi/180)
    # cam = sc.add_camera()
    # cam.resolution=1000
    tf = yt.ColorTransferFunction((0,3))
    tf.add_layers(10, w=0.01)


    # Set the bounds of the transfer function
    source = sc.sources['source_00']
    source.set_transfer_function(tf)
    source.tfh.set_log(False)
    sc.save(os.path.join(out_dir,'YS_volume.png'),sigma_clip=.25)
