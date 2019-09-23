'''
loadsNWUS11-S_percent.nc and plots some slices of dVs at different depths,
saves to provided directory.
'''
import yt
from yt_velmodel_vis import seis_model as sm
import numpy as np
import os
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt

if __name__=='__main__' or gogo:
    fname='NWUS11-S_percent.nc'
    out_dir='./output'
    model=sm.netcdf(fname)
    model.coordTransform('sphere2cart')
    data={}

    data['dvs']=model.data.variables['dvs'][:]
    # data['dvs'][data['dvs']>900.00000]=np.nan
    # data['dvs'][data['dvs']>=0]=np.nan
    #
    # data['dvs'][data['dvs']==np.nan]=0
    data['dvs'][data['dvs']>900.00000]=0
    data['dvs'][data['dvs']>=0]=0
    data['dvs']=data['dvs']*-1
    data['dvs'] = np.transpose(data['dvs'], (1, 2, 0))
    print(data['dvs'].max())

    # yt spherical expects R, theta, phi (depth, ~lat,~lon)
    sc_mult=1.0 # scale multiplier
    bbox = model.cart['bbox']
    print(bbox)
    ds = yt.load_uniform_grid(data,data['dvs'].shape,sc_mult,bbox=bbox,nprocs=1,
                            periodicity=(False,False,False),unit_system="mks")
    sc = yt.create_scene(ds,'dvs')
    # sc.annotate_axes()
    # cam = sc.add_camera()

    pos=sc.camera.position
    # pos[2]=np.mean(bbox[2])
    # pos[]
    sc.camera.set_position(pos,north_vector=np.array([0.0, 0.0, 1.0]))
    # cam.position = pos
    # cam.focus = np.array([np.mean(bbox[0]),np.mean(bbox[1]),np.mean(bbox[2])])
    # cam.north_vector = np.array([0.0, 0.0, 1.0])

    # sc.camera.zoom(.3)
    # wid=np.array([bbox[0][1]-bbox[0][0],bbox[1][1]-bbox[1][0],bbox[2][1]-bbox[2][0]])
    # wid=wid*0.1
    # wid_q=[ds.quan(wid[i],'m') for i in range(0,3)]
    # sc.camera.set_width(wid_q)#ds.quan(50*1000., 'm'))
    # sc.camera.north_vector = np.array([0.0, 0.0, 1.0])
    sc.camera.rotate(180*np.pi/180.)
    # sc.camera.roll(90.*np.pi/180.)
    # sc.camera.pitch(np.pi/4)
    # sc.camera.yaw(45.*np.pi/180)
    # cam = sc.add_camera()
    # cam.resolution=1000
    tf = yt.ColorTransferFunction((0, 5))
    tf.add_layers(4, w=0.01)

    source = sc.sources['source_00']
    source.tfh.set_log(False)
    source.set_transfer_function(tf)
    sc.save(os.path.join(out_dir,'WUS_volume.png'),sigma_clip=4)

    bbox = model.cart['bbox']
    ds = yt.load_uniform_grid(data,data['dvs'].shape,sc_mult,bbox=bbox,nprocs=1,
        periodicity=(True,True,True),unit_system="mks")


    center=np.array(ds.domain_center)
    dm={'left':np.array(ds.domain_left_edge),
        'right':np.array(ds.domain_right_edge)}


    dm['left']=dm['left']+200.*1000
    dm['right']=dm['right']-200*1000
    boxregion=ds.region(ds.domain_center,dm['left'],dm['right'])

    isoval=1
    surface = ds.surface(boxregion, "dvs", isoval)
    p3dc = Poly3DCollection(surface.triangles, linewidth=0.0)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.add_collection(p3dc)
    max_extent = (surface.vertices.max(axis=1) - surface.vertices.min(axis=1)).max()
    centers = (surface.vertices.max(axis=1) + surface.vertices.min(axis=1)) / 2
    bounds = np.zeros([3,2])
    bounds[:,0] = centers[:] - max_extent/2
    bounds[:,1] = centers[:] + max_extent/2
    print(bounds)
    ax.auto_scale_xyz(bounds[0,:], bounds[1,:], bounds[2,:])

    plt.savefig('output/WUS_isosurf_'+str(isoval)+'.png')
    # plt.close('all')
